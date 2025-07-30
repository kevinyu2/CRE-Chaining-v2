# Parallelizes chaining but NOT finding anchors

from pathlib import Path
from glob import glob
import sys
from chaining import chain_driver_np, chain_local_driver_np
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
from collections import defaultdict
from itertools import combinations, product
import os
from multiprocessing import Manager
from multiprocessing import shared_memory
import numpy as np
import time

'''
The one genome pipeline, parallelized, allowing varied scores and starting from ACR motif lists

This version finds chains between two sets

Uses preprocessed data (motifs for each genome for each ACR) and the motif pair data 
to output the anchors for each pair of ACRs across all genomes
'''

##########################################################
# Settings

# Set to 'full' or 'only_acr'; if the acr motifs were found from the full fimo, set to full
motif_mode = 'full'

# 'local' or 'global'
chain_mode = 'local'

# if motif mode is 'full', provide the ACR motif list. Otherwise, provide the FIMO folder
ACR_MOTIF_SET_A = '/home/mwarr/Data/One_Genome/exp3_randomforest/seta_rep_dapv1_motifs_preclust.txt'

ACR_MOTIF_SET_B = '/home/mwarr/Data/One_Genome/exp3_randomforest/setb_rand_dapv1_motifs_preclust.txt'

MATCH = 5
MISMATCH = -1
GAP = -1

# Set to None if not weighted
CUSTOM_SCORE = '/home/mwarr/Data/motif_relevance/scores_DAPv1_clustered.tsv'

# Adjust scores (good for local)
SCORE_CENTER = 5

OUTPUT = "/home/mwarr/Data/One_Genome/exp3_randomforest/Chaining_local_rep.tsv"


###########################################################
# Global memory variables

global_anchor_array = None
global_anchor_shm = None
global_anchor_name = f"anchors_{int(time.time())}"
BATCH_SIZE = 10000


###########################################################
# Helpers for parallelization

def chunkify(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

#######################################################
# Motif dict creation

# Helpers
def get_initial_motif_loc_only_acr(motif_loc_dict, data, prefix) :
    acr_set = set()
    pattern = os.path.join(data, 'fimo_out_*/', 'fimo.tsv')
    fimo_files = glob(pattern)

    for fimo_file in fimo_files:
        
        # Fill in the dict
        with open(fimo_file, "r") as f:
            # Ignore first line
            next(f)
            for line in f: 
                line_arr = line.rstrip().split('\t')
                # The file ends tsv info early
                if len(line_arr) < 5: 
                    break

                acr = prefix + line_arr[2]
                acr_set.add(acr)

                if int(line_arr[3]) not in motif_loc_dict[line_arr[0].split('-')[1]][acr] :
                    motif_loc_dict[line_arr[0].split('-')[1]][acr].append(int(line_arr[3]))
    return acr_set

def get_motif_loc_dict_full(motif_loc_dict, data, prefix) :
    acr_set = set()
    with open(data, 'r') as acr_motifs :
        curr_acr = ""
        curr_idx = 0
        for line in acr_motifs :
            if 'ACR: ' in line:
                curr_acr = prefix + line.rstrip().split('ACR: ')[1]
                curr_idx = 0
                acr_set.add(curr_acr)

            else :
                motif_loc_dict[line.rstrip()][curr_acr].append(curr_idx)
                curr_idx += 1

    return acr_set

# Takes in two xstreme folders / fimo folders, or files that list motif sequences  
# Gets motif locations
def get_motif_loc_dict(data1, data2) :
    # Holds where each motif is located {MOTIF: {acr: [loc, ..., loc] acr:[loc, ..., loc]}}
    motif_loc_dict = defaultdict(lambda: defaultdict(list))



    print(f"Filling Motif Dict, mode = {motif_mode}")


    if motif_mode == 'only_acr' :
        # Just keeps track of which motifs are in which set
        seta_acrs = get_initial_motif_loc_only_acr(motif_loc_dict, data1, "seta_")
        setb_acrs = get_initial_motif_loc_only_acr(motif_loc_dict, data2, "setb_")
       
        
        # Remove duplicates from repeated sequences (basically remove overlaps)
        for motif, single_motif_dict in motif_loc_dict.items() :
            motif_len = len(motif)
            for acr, loc_list in single_motif_dict.items() :

                loc_list.sort()
                to_remove = set()
                for i in range(len(loc_list) - 1) :
                    if loc_list[i + 1] - loc_list[i] <= motif_len :
                        to_remove.add(loc_list[i])
                single_motif_dict[acr] = [x for x in loc_list if x not in to_remove]
    
    else :
        # Open just the file, format is ACR: ...\n MOTIF \n MOTIF ...
        # Duplicates should have been removed already

        seta_acrs = get_motif_loc_dict_full(motif_loc_dict, data1, "seta_")
        setb_acrs = get_motif_loc_dict_full(motif_loc_dict, data2, "setb_")
        

    return seta_acrs, setb_acrs, motif_loc_dict

def get_scores() :
    score_dict = None
    if CUSTOM_SCORE != None :
        score_dict = {}
        with open(CUSTOM_SCORE, 'r') as cs :
           next(cs)
           for line in cs :
               score_dict[line.split('\t')[0]] = float(line.rstrip().split('\t')[1])

        # Adjust the score
        if SCORE_CENTER != None :
            current_mean = sum(score_dict.values()) / len(score_dict)
            scale = SCORE_CENTER / current_mean
            score_dict = {k: v * scale for k, v in score_dict.items()}

    return score_dict


# ONLY run on only_acr mode
def get_motif_loc_dict_local(inputa, inputb) :

    # This calls the original one, but we must change the indices 
    seta, setb, start_dict = get_motif_loc_dict(inputa, inputb)

    # Stores all locations of each acr as a set for numbering
    # {acr: {loc, loc, loc...}}
    acr_dict = defaultdict(set)

    for motif, single_motif_dict in start_dict.items() :
        for acr, loc_list in single_motif_dict.items() :
            acr_dict[acr].update(loc_list)

    # Matches positions to motif number
    # {acr: {start_pos : motif_no}}
    acr_no_dict = defaultdict(dict)
    for acr, locs in acr_dict.items() :
        acr_no_dict[acr] = {v: i for i, v in enumerate(sorted(locs))}
    
    # The localized motif dict to return
    motif_loc_dict = defaultdict(lambda: defaultdict(list))
    
    for motif, single_motif_dict in start_dict.items() :
        for acr, loc_list in single_motif_dict.items() :
            for loc in loc_list :
                motif_loc_dict[motif][acr].append(acr_no_dict[acr][loc])

    return seta, setb, motif_loc_dict

#####################################################################
# Calculate number of anchors


def init_anchor_array(total_anchors):
    global global_anchor_array, global_anchor_shm
    if CUSTOM_SCORE == None :
        global_anchor_shm = shared_memory.SharedMemory(create=True, size=total_anchors * 2 * 4, name = global_anchor_name)  # 2 ints per anchor, 4 bytes each
        print(f"Allocating {total_anchors * 8} Bytes")
        global_anchor_array = np.ndarray((total_anchors, 2), dtype=np.int32, buffer=global_anchor_shm.buf)
    else :
        global_anchor_shm = shared_memory.SharedMemory(create=True, size=total_anchors * 3 * 4, name = global_anchor_name)  # 3 floats per anchor, 8 bytes each
        print(f"Allocating {total_anchors * 12} Bytes")
        global_anchor_array = np.ndarray((total_anchors, 3), dtype='float32', buffer=global_anchor_shm.buf)


def cleanup_anchor_array():
    global global_anchor_shm
    global_anchor_shm.close()
    global_anchor_shm.unlink()

# Calculate the total number of anchors
# Returns an int, the number of anchors, and a dict, which maps each acr to its start and stop points in the anchor space
def calculate_no_anchors(seta, setb, motif_loc_dict) :
    print("Calculating Anchor Dict Allocation Space")


    # Running total needed room
    total = 0
    # To track individual space needed
    # {(acr1, acr2) : space}
    acr_space_dict = defaultdict(int)
    # To give start and stop indices
    # {(acr1, acr2) : (start, stop)}
    acr_start_stop_dict = defaultdict(lambda: (0, 0))
    for motif_name, single_motif_dict in tqdm(motif_loc_dict.items()) :
        sizes = {acr: len(single_motif_dict[acr]) for acr in list(single_motif_dict.keys())}


        for acr1 in seta & set(single_motif_dict.keys()) :
            
            for acr2 in setb & set(single_motif_dict.keys()) :
                key = (acr1, acr2)
                pair_count = sizes[acr1] * sizes[acr2]
                acr_space_dict[key] += pair_count
                total += pair_count


    curr_start = 0
    for key, space in acr_space_dict.items() :
        acr_start_stop_dict[key] = (curr_start, curr_start + space)
        curr_start += space
    
    return total, acr_start_stop_dict


######################################################################
# Anchors

# Find anchors given a loc dict. For local and global
def find_anchors(seta, setb, anchor_loc_dict, motif_loc_dict, score_dict) :
    print("Finding Anchors")

    curr_filled_per_pair = defaultdict(int)

    for motif_name, single_motif_dict in tqdm(motif_loc_dict.items()) :

        for acr1 in seta & set(single_motif_dict.keys()) :
            for acr2 in setb & set(single_motif_dict.keys()) :

                key = (acr1, acr2)

                anchors1 = single_motif_dict[acr1]
                anchors2 = single_motif_dict[acr2]

                # Write to the correct slice of the global array
                n1, n2 = len(anchors1), len(anchors2)
                n_total = n1 * n2

                if CUSTOM_SCORE == None :
                    start_idx = anchor_loc_dict[key][0] + curr_filled_per_pair[key]

                    i = 0
                    for a1 in anchors1:
                        for a2 in anchors2:
                            global_anchor_array[start_idx + i, 0] = a1
                            global_anchor_array[start_idx + i, 1] = a2
                            i += 1

                    # One write
                    curr_filled_per_pair[key] += n_total

                else :
                    start_idx = anchor_loc_dict[key][0] + curr_filled_per_pair[key]

                    i = 0
                    for a1 in anchors1:
                        for a2 in anchors2:
                            global_anchor_array[start_idx + i, 0] = a1
                            global_anchor_array[start_idx + i, 1] = a2
                            if motif_name in score_dict :
                                global_anchor_array[start_idx + i, 2] = score_dict[motif_name]
                            else :
                                global_anchor_array[start_idx + i, 2] = MISSING_VALUE

                            i += 1

                    curr_filled_per_pair[key] += n_total

####################################################################
# Global chaining


# Parallelized chaining worker
def batch_pair_chaining(args):
    pairs, total_anchors, custom_score, global_anchor_name = args
    shm = shared_memory.SharedMemory(name=global_anchor_name)

    # Create NumPy array wrapper (no copy)
    if custom_score :
        shared_array = np.ndarray((total_anchors, 3), dtype='float32', buffer=shm.buf)
    else :
        shared_array = np.ndarray((total_anchors, 2), dtype=np.int32, buffer=shm.buf)
    results = []
    for acr1, acr2, start, stop in pairs:  
        chain_len = chain_driver_np(shared_array[start:stop], custom_score)
        results.append(f"{acr1[5:]}\t{acr2[5:]}\t{chain_len}\t{stop - start}\n")

    return results

# Finds pairwise chain length (global)
# Uses global anchor dict dict from earlier and the output directory     
# Make sure to clear the output folder first since we are appending to files 
# Printed file in out_file in the format acr1\tacr2\tchain_len\tno_anchors
def chain_all_pairs(anchor_loc_dict, total_anchors, out_file) :

    print("Chaining")

    # Parallelized chaining
    anchor_items = [(a, b, start, stop) for (a, b), (start, stop) in anchor_loc_dict.items()]

    batches = list(chunkify(anchor_items, BATCH_SIZE))
    parallelized_inputs =  [(batch, total_anchors, CUSTOM_SCORE != None, global_anchor_name) for batch in batches]

    with Pool(cpu_count()) as pool, open(out_file, 'w') as out:
        try :
            for lines in tqdm(pool.imap_unordered(batch_pair_chaining, parallelized_inputs), total=len(batches)):
                out.writelines(lines)
        finally:
            pool.close()
            pool.join()           


#######################################################################
# Local chaining


# Parallelized chaining worker
def batch_pair_local_chaining(args):
    pairs, total_anchors, match, mismatch, gap, custom_score, global_anchor_name = args
    shm = shared_memory.SharedMemory(name=global_anchor_name)

    if custom_score :
        shared_array = np.ndarray((total_anchors, 3), dtype='float32', buffer=shm.buf)

    else :
        shared_array = np.ndarray((total_anchors, 2), dtype=np.int32, buffer=shm.buf)

    results = []
    for acr1, acr2, start, stop in pairs:
        chain_results = chain_local_driver_np(shared_array[start:stop], match, mismatch, gap, custom_score)
        results.append(f"{acr1[5:]}\t{acr2[5:]}\t{chain_results[0]}\t{chain_results[1]}\n")

    return results


# Local version
# Printed file in out_file in the format acr1\tacr2\tchain_score\tchain_len\tno_anchors
def chain_all_pairs_local(anchor_loc_dict, total_anchors, match, mismatch, gap, out_file) :

    print("Chaining")
    anchor_items = [(a, b, start, stop) for (a, b), (start, stop) in anchor_loc_dict.items()]
    batches = list(chunkify(anchor_items, BATCH_SIZE))

    # Add match, mismatch, and gap
    parallelized_inputs =  [(batch, total_anchors, match, mismatch, gap, CUSTOM_SCORE != None, global_anchor_name) for batch in batches]
    with Pool(cpu_count()) as pool, open(out_file, 'w') as out:
        try :
            for lines in tqdm(pool.imap_unordered(batch_pair_local_chaining, parallelized_inputs), total=len(parallelized_inputs)):
                out.writelines(lines)
        finally :
            pool.close()
            pool.join()


##############################################################
# Drivers



def chain_global_driver(inputa, inputb, out_file) :
    score_dict = get_scores()
    seta, setb, motif_l_dict = get_motif_loc_dict(inputa, inputb)
    space, anchor_loc_dict = calculate_no_anchors(seta, setb, motif_l_dict)
    init_anchor_array(space)
    find_anchors(seta, setb, anchor_loc_dict, motif_l_dict, score_dict)
    chain_all_pairs(anchor_loc_dict, space, out_file)
    cleanup_anchor_array()

def chain_local_driver(inputa, inputb, out_file, match, mismatch, gap) :
    score_dict = get_scores()
    seta, setb, motif_l_dict = get_motif_loc_dict_local(inputa, inputb)
    space, anchor_loc_dict = calculate_no_anchors(seta, setb, motif_l_dict)
    init_anchor_array(space)
    find_anchors(seta, setb, anchor_loc_dict, motif_l_dict, score_dict)
    chain_all_pairs_local(anchor_loc_dict, space, match, mismatch, gap, out_file)
    cleanup_anchor_array()


if chain_mode == 'local' :
    chain_local_driver(ACR_MOTIF_SET_A, ACR_MOTIF_SET_B, OUTPUT, MATCH, MISMATCH, GAP)

else :
    chain_global_driver(ACR_MOTIF_SET_A, ACR_MOTIF_SET_B, OUTPUT)




