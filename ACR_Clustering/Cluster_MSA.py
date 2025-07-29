# Get a sequence of motifs for each cluster, as a "representative"
# Uses a central star method
import sys
import os
from collections import defaultdict
script_dir = os.path.dirname(os.path.abspath(__file__))
target_dir = os.path.abspath(os.path.join(script_dir, ".."))
sys.path.append(target_dir)
from chaining import chain_local, chain, chain_local_weighted, chain_weighted
from collections import Counter

########################################################################
# Settings

OUTPUT_FILE = '/home/mwarr/Data/Clustering/acr_consensus_DAPv1_clustered_seta_local.txt'

# Result from clustering
CLUSTER_FILE = '/home/mwarr/Data/Clustering/hierarchical_DAPv1_clustered_seta_local.txt'

# List of ACRs
ACR_SEQUENCE_FILE = '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/post_fimo/pairwise_inputs/acr_DAPv1_clustered.txt'

# For determining direction, local or global
MODE = "local"

MATCH = 5
MISMATCH = -1
GAP = -1

# Set to None if not weighted
CUSTOM_SCORE = '/home/mwarr/Data/motif_relevance/scores_DAPv1_clustered.tsv'

# Adjust scores (good for local)
SCORE_CENTER = 5

# Min number of CREs to be put in the output
MIN_CRE_LEN = 5

############################################################################
# Alignment

def central_star_msa(central, others) :
    aligned_central = central
    aligned_others = []

    # Keep all alignments aligned to central
    for seq in others:
        aligned1, aligned2 = needleman_wunsch(aligned_central, seq)

        # Expand previous alignments to match new alignment
        if len(aligned_others) == 0:
            aligned_central = aligned1
            aligned_others.append(aligned2)
        else:
            # Expand all sequences to match current central alignment
            aligned_central, aligned1 = needleman_wunsch(aligned_central, aligned1)
            new_aligned_others = []
            for old in aligned_others:
                _, aligned_old = needleman_wunsch(aligned_central, old)
                new_aligned_others.append(aligned_old)
            _, aligned2 = needleman_wunsch(aligned_central, aligned2)

            aligned_others = new_aligned_others + [aligned2]

    return [aligned_central] + aligned_others

# TODO: weighted scoring
from typing import List, Tuple

GAP_CHAR = "-"


def score(a: str, b: str) -> int:
    if a == b:
        return MATCH
    return MISMATCH

def needleman_wunsch(seq1: List[str], seq2: List[str]) -> Tuple[List[str], List[str]]:
    m, n = len(seq1), len(seq2)
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Init
    for i in range(m + 1):
        dp[i][0] = i * GAP
    for j in range(n + 1):
        dp[0][j] = j * GAP

    # Fill
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = dp[i - 1][j - 1] + score(seq1[i - 1], seq2[j - 1])
            delete = dp[i - 1][j] + GAP
            insert = dp[i][j - 1] + GAP
            dp[i][j] = max(match, delete, insert)

    # Traceback
    aligned1, aligned2 = [], []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] + score(seq1[i-1], seq2[j-1]):
            aligned1.append(seq1[i-1])
            aligned2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and dp[i][j] == dp[i-1][j] + GAP:
            aligned1.append(seq1[i-1])
            aligned2.append(GAP_CHAR)
            i -= 1
        else:
            aligned1.append(GAP_CHAR)
            aligned2.append(seq2[j-1])
            j -= 1

    return aligned1[::-1], aligned2[::-1]

def extract_consensus(msa: List[List[str]], gap_symbol: str = "-") -> List[str]:
    consensus = []
    num_cols = len(msa[0])
    
    for col in range(num_cols):
        column = [seq[col] for seq in msa]
        counts = Counter(column)
        most_common_token, c = counts.most_common(1)[0]
        
        if most_common_token == gap_symbol or c <= (len(msa) / 2):
            continue  # skip this column entirely
        else:
            consensus.append(most_common_token)
    
    return consensus

################################################################################

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

if CUSTOM_SCORE != None:
    score_dict = get_scores()



acr_sequence_dict = defaultdict(list)

with open(ACR_SEQUENCE_FILE, 'r') as f:
    curr_ACR = ""
    for line in f:
        if "ACR: " in line :
            curr_ACR = line.rstrip().split('ACR: ')[1]
        else :
            acr_sequence_dict[curr_ACR].append(line.rstrip())

with open(CLUSTER_FILE, 'r') as f:
    with open(OUTPUT_FILE, 'w') as out:
        for line1, line2 in zip(f, f):
            representative = line1.rstrip().split('Representative: ')[1]
            cluter_name = line1.rstrip().split('Cluster no: ')[1].split(', Representative: ')[0]
            cluster = line2.rstrip().split('\t')
            representative_seq = acr_sequence_dict[representative]

            cluster_seq = [acr_sequence_dict[c] for c in cluster if c != representative]

            # Reverse any that chain the other way
            representative_locations = defaultdict(list)
            for i, cre in enumerate(representative_seq) :

                representative_locations[cre].append(i)


            for seq_idx, seq in enumerate(cluster_seq) :
                anchors = []

                # find anchors
                for i, cre in enumerate(seq) :
                    if cre in representative_locations:
                        for j in representative_locations[cre]:
                            if CUSTOM_SCORE == None:
                                anchors.append((i, j))
                            else :
                                anchors.append((i, j, score_dict[cre]))


                if MODE == 'global' and CUSTOM_SCORE == None:
                    reverse_anchors = [(anchor[0], -1 * anchor[1]) for anchor in anchors]

                    if chain(anchors) < chain(reverse_anchors) :
                        cluster_seq[seq_idx] = seq[::-1]

                if MODE == 'local' and CUSTOM_SCORE == None:
                    reverse_anchors = [(anchor[0], -1 * anchor[1]) for anchor in anchors]

                    if chain_local(anchors) < chain_local(reverse_anchors) :
                        cluster_seq[seq_idx] = seq[::-1]


                if MODE == 'local' and CUSTOM_SCORE != None:
                    reverse_anchors = [(anchor[0], -1 * anchor[1], anchor[2]) for anchor in anchors]

                    if chain_local_weighted(anchors, MATCH, MISMATCH, GAP) < chain_local_weighted(reverse_anchors, MATCH, MISMATCH, GAP) :
                        cluster_seq[seq_idx] = seq[::-1]

                if MODE == 'global' and CUSTOM_SCORE != None:
                    reverse_anchors = [(anchor[0], -1 * anchor[1], anchor[2]) for anchor in anchors]

                    if chain_weighted(anchors, MATCH, MISMATCH, GAP) < chain_weighted(reverse_anchors, MATCH, MISMATCH, GAP) :
                        cluster_seq[seq_idx] = seq[::-1]



            

            new_seq = extract_consensus(central_star_msa(representative_seq, cluster_seq))

            if len(new_seq) >= MIN_CRE_LEN :
                out.write(f"ACR: cluster_{cluter_name}")
                out.write('\n')

                for cre in new_seq :
                    out.write(cre)
                    out.write('\n')

# print(central_star_msa(['a', 'b', 'c'], [['a', 'c'], ['b', 'c']]))