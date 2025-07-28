
from collections import defaultdict
from pathlib import Path
from glob import glob
from itertools import islice
import numpy as np
import time
from multiprocessing import shared_memory
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

##############################################################################
# Settings

# Motif sequences for chromosomes
INPUT_FOLDER = '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/post_fimo/BLAST_inputs/v2ArabidopsisPBM_2014021_fimo'

# ACR motif sequences
ACR_FILE = '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/post_fimo/BLAST_inputs/v2acr_PBM_motifs.txt'

OUTFILE = '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/post_fimo/BLAST_hits/hits_PBM_scorev2l08.tsv'

# Initial kmer size to set as a seed
KMER_SIZE = 4

# ACR set to consider. If all, set to None
ACR_SET = '/home/mwarr/Data/seta_half.txt'
# ACR_SET = None

# Ignore BREAK motifs
IGNORE_BREAK = True

# How much larger the sequence window is to the ACR window
WINDOW_SIZE_MULTIPLIER = 2

# Scoring
MATCH = 5
MISMATCH = -1
GAP = -1

# Score cutoff
SCORE_THRESH = 0

# If none, then score is set to the MATCH score. Otherwise, provide file
CUSTOM_SCORE = '/home/mwarr/Data/motif_relevance/scores_PBM_lambda08.tsv'
# CUSOM_SCORE = None

# If not none, sets custom score's average to this number
SCORE_CENTER = 5


###########################################################################
# Global memory variables

# [motif, motif...]
global_motif_array = None
global_motif_shm = None
global_motif_name = 'motif'

# [motif, motif...]
global_acr_array = None
global_acr_shm = None
global_acr_name = 'acr'

BATCH_SIZE = 2000

# Fills the global motif array and kmer array, returns the location dict
def parse_input(input_folder) :
    # {(KMER) : [(chromosome, index), (chromosome, index)]}
    kmer_dict = defaultdict(list)

    # {(chromosome, index) : location}
    location_dict = {}

    # {chromosome: [motif, motif, ..., motif]}
    motif_lists = defaultdict(list)

    # {(chromosome, index): end}
    end_dict = {}

    # Remove if still exists
    try:
        shm = shared_memory.SharedMemory(name='motif')
        shm.unlink()  # This removes the shared memory segment
        shm.close()
    except:
        _ = 5
    try: 
        shm = shared_memory.SharedMemory(name='acr')
        shm.unlink()  # This removes the shared memory segment
        shm.close()
    except:
        _ = 5


    print("Reading Input Sequence")


    # Create the motif lists and location_dict
    for file_path in Path(input_folder).iterdir():
        chromosome = str(file_path).rstrip().split('/')[-1]
        with open(file_path, 'r') as motif_file :

            motif_count = 0
            for line in motif_file :
                if IGNORE_BREAK :
                    if 'BREAK\n' not in line:
                        # should be [motif, location]
                        line_arr = line.rstrip().split('\t')

                        # Fill locaiton dict
                        location_dict[(chromosome, motif_count)] = int(line_arr[1])

                        # fill in motif_list
                        motif_lists[chromosome].append(line_arr[0])

                        end_dict[(chromosome, motif_count)] = int(line_arr[2])

                        motif_count += 1

    # Convert to ints
    motif_to_int = {}
    int_to_motif = {}
    curr_motif_idx = 0

    # Make them numpy instead, and use ints
    for chr, motif_list in motif_lists.items() :
        np_int_motif_list = np.zeros(len(motif_list))
        for i, motif in enumerate(motif_list) :
            # Found a new motif
            if motif not in motif_to_int :
                motif_to_int[motif] = curr_motif_idx
                int_to_motif[curr_motif_idx] = motif
                curr_motif_idx += 1
            
            np_int_motif_list[i] = motif_to_int[motif]
        motif_lists[chr] = np_int_motif_list

                        
    for chr, motif_list in motif_lists.items() :
        iterators = [islice(motif_list, i, None) for i in range(KMER_SIZE)]
        for idx, kmer in enumerate(zip(*iterators)) :
            kmer_dict[kmer].append((chr, idx))

    score_dict = {}
    if CUSTOM_SCORE != None :
        with open(CUSTOM_SCORE, 'r') as cs :
           next(cs)
           for line in cs :
               score_dict[motif_to_int[line.split('\t')[0]]] = float(line.rstrip().split('\t')[1])

        # Adjust the score
        if SCORE_CENTER != None :
            current_mean = sum(score_dict.values()) / len(score_dict)
            scale = SCORE_CENTER / current_mean
            score_dict = {k: v * scale for k, v in score_dict.items()}

    else :
        for int_motif in int_to_motif.keys() :
            score_dict[int_motif] = MATCH


    return kmer_dict, location_dict, motif_lists, motif_to_int, int_to_motif, end_dict, score_dict

# Turn into global array
# Returns mappings from motif to int, int to motif, and a place where all the motifs are
def init_global_motif(motif_lists) :
    print("Initializing Global Motif")

    # First calculate the size
    motif_array_size = 0
    for motif_list in motif_lists.values() :
        motif_array_size += len(motif_list)

    global global_motif_array, global_motif_shm
    global_motif_shm = shared_memory.SharedMemory(create=True, size=motif_array_size * 4, name = global_motif_name)  # 4 bytes each
    global_motif_array = np.ndarray(motif_array_size, dtype=np.int32, buffer=global_motif_shm.buf)

    # Maps motif location to shared array location
    motif_array_locs = {}

    # Tracks motif number
    curr_array_idx = 0
    for chr, motif_list in motif_lists.items() :
        motif_array_locs[chr] = (curr_array_idx, curr_array_idx + len(motif_list))
        for motif in motif_list :
            
            global_motif_array[curr_array_idx] = motif
            curr_array_idx += 1

    return motif_array_size, motif_array_locs
    

# Final cleanup step for shared memory
def cleanup_global_arrays():
    global_motif_shm.close()
    global_motif_shm.unlink()

    global_acr_shm.close()
    global_acr_shm.unlink()


def batched_worker(args) :
    batch, motif_array_size, acr_array_size, score_dict, MISMATCH, GAP, SCORE_THRESH = args

    results = []

    # Open shared memory
    motif_shm = shared_memory.SharedMemory(name='motif')
    shared_motif_array = np.ndarray(motif_array_size, dtype=np.int32, buffer=motif_shm.buf)

    acr_shm = shared_memory.SharedMemory(name = 'acr')
    shared_acr_array = np.ndarray(acr_array_size, dtype=np.int32, buffer=acr_shm.buf)

    for motif_arr_start, motif_arr_stop, acr_arr_start, acr_arr_stop, acr, chr in batch :
        local_score, local_start, local_stop = local_align(shared_motif_array[motif_arr_start : motif_arr_stop], shared_acr_array[acr_arr_start : acr_arr_stop], score_dict)

        if local_score > SCORE_THRESH :
            results.append((acr, local_score, chr, motif_arr_start, motif_arr_stop, local_start, local_stop))

    return results

def chunkify(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def find_matches(motif_lists, acr_file, outfile, kmer_dict, motif_array_locs, motif_array_size, location_dict, motif_to_int, end_dict, score_dict) :
    print("Parsing ACR File")

    # If restricting by set
    if ACR_SET != None :
        acrs_to_use = set()
        with open(ACR_SET, 'r') as acr_set :
            for line in acr_set :
                acrs_to_use.add(line.rstrip())

    # {ACR: [motif, motif...]}
    acr_dict = defaultdict(list)
    with open(acr_file, 'r') as acrf :
        curr_acr = ""
        for line in acrf :
            # End of one
            if 'ACR: ' in line:
                # RESET
                curr_acr = line.rstrip().split('ACR: ')[1]

            else :
                if ACR_SET != None :
                    if curr_acr in acrs_to_use :
                        # Check if motif is found in the dict, if not, safe to just make it -1
                        if line.rstrip() in motif_to_int :
                            acr_dict[curr_acr].append(motif_to_int[line.rstrip()])
                        else :
                            acr_dict[curr_acr].append(-1)

                else :
                    if line.rstrip() in motif_to_int :
                        acr_dict[curr_acr].append(motif_to_int[line.rstrip()])
                    else :
                        acr_dict[curr_acr].append(-1)

    # Get reverse:
    for acr in list(acr_dict.keys()) :
        acr_dict[f"{acr}_rev"] = acr_dict[acr][::-1]


    # Put into global array
    # First calculate the size
    acr_array_size = 0
    for motif_list in acr_dict.values() :
        acr_array_size += len(motif_list)

    global global_acr_array, global_acr_shm
    global_acr_shm = shared_memory.SharedMemory(create=True, size=acr_array_size * 4, name = global_acr_name)  # 4 bytes each
    global_acr_array = np.ndarray(acr_array_size, dtype=np.int32, buffer=global_acr_shm.buf)

    # fill the array
    acr_array_locs = {}
    # Tracks array start
    curr_array_idx = 0
    for acr, motif_list in acr_dict.items() :
        acr_array_locs[acr] = (curr_array_idx, curr_array_idx + len(motif_list))
        for motif in motif_list :

            global_acr_array[curr_array_idx] = motif
            curr_array_idx += 1

    print("Finding Matches")

            
    # Actually match the seeds and extend
    with open(outfile, 'w') as out :
        # Print out info
        out.write(f"kmer_length: {KMER_SIZE}\twindow_size: {WINDOW_SIZE_MULTIPLIER}\tscore_thresh: {SCORE_THRESH}\tmatch: {MATCH}\tmismatch: {MISMATCH}\tgap: {GAP}\tignore_break: {IGNORE_BREAK}\tacr_set: {ACR_SET}\n")
        out.write("ACR\tScore\tChr\tStart\tStop\n")  

        # Create parallelized inputs
        # [(motif_arr_start, motif_arr_stop, acr_arr_start, acr_arr_stop, acr, chr)]
        parallelized_inputs = []

        for acr, acr_motif_list in acr_dict.items() :
            iterators = [islice(acr_motif_list, i, None) for i in range(KMER_SIZE)]
            len_acr_motif = len(acr_motif_list)
            for acr_idx, kmer in enumerate(zip(*iterators)) :

                if kmer in kmer_dict :
                    for chr, sequence_idx in kmer_dict[kmer] :
                        motif_arr_start = max(motif_array_locs[chr][0], motif_array_locs[chr][0] + sequence_idx - (acr_idx * WINDOW_SIZE_MULTIPLIER))
                        motif_arr_stop = min(motif_array_locs[chr][1], motif_array_locs[chr][0] + sequence_idx + KMER_SIZE + ((len_acr_motif - (acr_idx + KMER_SIZE)) * WINDOW_SIZE_MULTIPLIER))
                        
                        # print(f"{chr}, {sequence_idx} | {acr}, {acr_idx} : {kmer} | {motif_arr_start} to {motif_arr_stop}")

                        parallelized_inputs.append((motif_arr_start, motif_arr_stop, acr_array_locs[acr][0], acr_array_locs[acr][1], acr, chr))

        print(f"Total matches: {len(parallelized_inputs)}")
        batches = list(chunkify(parallelized_inputs, BATCH_SIZE))
        final_inputs =  [(batch, motif_array_size, acr_array_size, score_dict, MISMATCH, GAP, SCORE_THRESH) for batch in batches]
        with Pool(cpu_count()) as pool :
            try :
                for lines in tqdm(pool.imap_unordered(batched_worker, final_inputs), total=len(final_inputs)):
                    for line in lines :
                        acr, local_score, chr, start, stop, local_start, local_stop = line
                        out.write(f"{acr}\t{local_score}\t{chr}\t{location_dict[(chr, start + local_start - motif_array_locs[chr][0])]}\t{end_dict[(chr, start + local_stop - motif_array_locs[chr][0])]}\n")
            finally :
                pool.close()
                pool.join()



# SW alignment
def local_align(list_seq, list_acr, score_dict) :
    l1 = len(list_seq)
    l2 = len(list_acr)
    score_matrix = np.zeros((l1 + 1, l2 + 1), dtype=float)
    
    max_score = 0
    max_pos = (0, 0)

    # Fill the scoring matrix
    for i in range(1, l1 + 1):
        for j in range(1, l2 + 1):
            if list_seq[i - 1] == list_acr[j - 1]:
                score = score_dict[list_seq[i - 1]]
            else:
                score = MISMATCH

            diag = score_matrix[i - 1][j - 1] + score
            up   = score_matrix[i - 1][j] + GAP
            left = score_matrix[i][j - 1] + GAP

            score_matrix[i][j] = max(0, diag, up, left)

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    # Traceback to find start position
    i, j = max_pos
    end_i, end_j = i, j

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        current = score_matrix[i][j]
        if list_seq[i - 1] == list_acr[j - 1]:
            score = score_dict[list_seq[i - 1]]
        else:
            score = MISMATCH

        if current == score_matrix[i - 1][j - 1] + score:
            i -= 1
            j -= 1
        elif current == score_matrix[i - 1][j] + GAP:
            i -= 1
        else:
            j -= 1

    start_i, start_j = i, j

    return max_score, start_i, end_i - 1

##################################################################################
# Drivers

def run_all(input_folder, acr_file, outfile) :
    kmer_dict, location_dict, motif_lists, motif_to_int, int_to_motif, end_dict, score_dict = parse_input(input_folder)
    motif_array_size, motif_array_locs = init_global_motif(motif_lists)
    find_matches(motif_lists, acr_file, outfile, kmer_dict, motif_array_locs, motif_array_size, location_dict, motif_to_int, end_dict, score_dict)
    cleanup_global_arrays()

################################################################################


run_all(INPUT_FOLDER, ACR_FILE, OUTFILE)
