
from collections import defaultdict
from pathlib import Path
from glob import glob
from itertools import islice
import numpy as np
import time
from multiprocessing import shared_memory
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from Bio import Align


##############################################################################
# Settings

# Additional chain scoring threshold for cutoff
SCORE_THRESH = 40

# How similar to be to ignore it
ALIGNMENT_SCORE_CUTOFF = None

SEQUENCE_FILE = '/home/projects/msu_nsf_pangenomics/pgrp/data/arabidopsis/atacseq/tair10.fa'
# SEQUENCE_FILE = '/home/kyu/Test_scripts/test_full_genome.fa'

# Hits
HITS_FILE = '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/post_fimo/BLAST_hits/hits_DAPv1_scorev2l08.tsv'
# HITS_FILE = '/home/kyu/Test_scripts/test_hits.tsv'

chrs_to_use = {'Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5'}

# OUT_FILE = '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/hits_ArabidopsisPBM_with_align.tsv'
OUT_FILE = '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/post_fimo/BLAST_hits/hits_DAPv1_scorev2l08_aligned.tsv'


##############################################################################
# [nt, nt...]
global_seq_array = None
global_seq_shm = None
global_seq_name = 'seq'

BATCH_SIZE = 10000

#############################################################################

def parse_sequence() :
    print("Reading Sequences")

    sequence_dict = defaultdict(list)
    sequence_sizes = defaultdict(int)
    with open(SEQUENCE_FILE, 'r') as seq_file :
        curr_chr = ""
        to_use = True
        for line in seq_file :
            if '>' in line :
                curr_chr = line.rstrip().split('>')[1]
                if curr_chr in chrs_to_use :
                    to_use = True
                else :
                    to_use = False
            elif to_use :
                sequence_dict[curr_chr].append(line.rstrip())
                sequence_sizes[curr_chr] += len(line.rstrip())
    
    total_size = 0
    for single_size in sequence_sizes.values() :
        total_size += single_size


    try: 
        shm = shared_memory.SharedMemory(name='seq')
        shm.unlink()  # This removes the shared memory segment
        shm.close()
    except:
        _ = 5

    print("Creating Shared Array")

    # Create the shared memory
    global global_seq_array, global_seq_shm
    global_seq_shm = shared_memory.SharedMemory(create=True, size=total_size * 1, name = global_seq_name)  # 1 byte each
    global_seq_array = np.ndarray(total_size, dtype='S1', buffer=global_seq_shm.buf)

    # {Chr: (start, end)}
    seq_loc_dict = {}
    curr_start = 0
    for chr, size in sequence_sizes.items() :
        full_seq = ''.join(sequence_dict[chr])
        global_seq_array[curr_start:curr_start + size] = np.frombuffer(bytes(full_seq, encoding='ascii'), dtype='S1')
        seq_loc_dict[chr] = (curr_start, curr_start + size)
        curr_start += size

    return seq_loc_dict, total_size
    


# Final cleanup step for shared memory
def cleanup_global_arrays():
    global_seq_shm.close()
    global_seq_shm.unlink()

# Read in inputs, make queries list, and list of lines to save
def parse_input(seq_loc_dict) :
    print("Reading Hit File")
    # To ensure no repeats
    lines_dict = set()

    # (unbatched) inputs: [(line_idx, qstart, qend, rstart, rend, isReverse)]
    queries = []
    lines_hits = []
    curr_line_idx = 0
    with open(HITS_FILE, 'r') as hits :
        next(hits)
        next(hits)
        for line in hits :
            if line not in lines_dict :
                lines_dict.add(line)
                line_arr = line.rstrip().split('\t')

                if SCORE_THRESH == None or float(line_arr[1]) >= SCORE_THRESH :
                    qchr = line_arr[0].split('_')[0]
                    qstart = int(line_arr[0].split('_')[1].split('to')[0])

                    # Account for reversed versions
                    if 'rev' not in line_arr[0] :
                        qend = int(line_arr[0].rstrip().split('to')[1])
                        if qchr in chrs_to_use and line_arr[2] in chrs_to_use :
                            queries.append((curr_line_idx, 
                                            qstart + seq_loc_dict[qchr][0], 
                                            qend + seq_loc_dict[qchr][0], 
                                            int(line_arr[3]) + seq_loc_dict[line_arr[2]][0], 
                                            int(line_arr[4])  + seq_loc_dict[line_arr[2]][0], False))

                            lines_hits.append(line.rstrip())
                            curr_line_idx += 1
                    else :
                        qend = int(line_arr[0].split('to')[1].split('_')[0])
                        if qchr in chrs_to_use and line_arr[2] in chrs_to_use :
                            queries.append((curr_line_idx, 
                                            qstart + seq_loc_dict[qchr][0], 
                                            qend + seq_loc_dict[qchr][0], 
                                            int(line_arr[3]) + seq_loc_dict[line_arr[2]][0], 
                                            int(line_arr[4])  + seq_loc_dict[line_arr[2]][0], True))

                            lines_hits.append(line.rstrip())
                            curr_line_idx += 1

    batched_queries = list(chunkify(queries, BATCH_SIZE))

    return batched_queries, lines_hits


def chunkify(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def batched_worker(args) :
    batch, seq_arr_size = args

    #[(line_idx, identity)]
    results = []
    # Open shared memory
    seq_shm = shared_memory.SharedMemory(name='seq')
    shared_motif_array = np.ndarray(seq_arr_size, dtype='S1', buffer=seq_shm.buf)

    aligner = Align.PairwiseAligner(mode = 'global', scoring = "blastn", query_end_open_gap_score = 0, query_end_extend_gap_score = 0)
    for query in batch :
        if query[5] :
            alignment = aligner.align(shared_motif_array[query[1]: query[2] + 1].astype('S').tobytes().decode('ascii')[::-1], 
                                  shared_motif_array[query[3]: query[4] + 1].astype('S').tobytes().decode('ascii'))
        else :
            alignment = aligner.align(shared_motif_array[query[1]: query[2] + 1].astype('S').tobytes().decode('ascii'), 
                                  shared_motif_array[query[3]: query[4] + 1].astype('S').tobytes().decode('ascii'))
        best = alignment[0]
        # print(best)
        # print(f"Score: { best.score}")

        results.append((query[0], best.score /  (query[4] + 1 - query[3])))
    return results

def align_all_semi_global(batched_queries, seq_arr_size, lines_hits) :
    print("Aligning")

    parallelized_inputs = [(batch, seq_arr_size) for batch in batched_queries]
    with open(OUT_FILE, 'w') as out :
        with open(HITS_FILE, 'r') as hits :
            line = hits.readline().rstrip()
            out.write(f"{line}\tScore_thresh_align: {SCORE_THRESH}\tAlign_score_max: {ALIGNMENT_SCORE_CUTOFF}\n")

            line = hits.readline().rstrip()
            out.write(f"{line}\tAlign_score\n")


        with Pool(cpu_count()) as pool :
            try :
                for lines in tqdm(pool.imap_unordered(batched_worker, parallelized_inputs), total=len(parallelized_inputs)):
                    for line in lines :
                        if ALIGNMENT_SCORE_CUTOFF != None :
                            if line[1] <= ALIGNMENT_SCORE_CUTOFF :
                                out.write(f"{lines_hits[line[0]]}\t{line[1]}\n")
                        else :

                            out.write(f"{lines_hits[line[0]]}\t{line[1]}\n")

            finally :
                pool.close()
                pool.join()

seq_loc_dict, seq_arr_size = parse_sequence()
batched_queries, lines_hits = parse_input(seq_loc_dict)
align_all_semi_global(batched_queries, seq_arr_size, lines_hits)
cleanup_global_arrays()