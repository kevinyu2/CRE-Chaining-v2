import pandas as pd
from pathlib import Path
import numpy as np
from collections import defaultdict
import time
from tqdm import tqdm

# To adjust score if necessary (i.e. if scores were averaged to 5, make this 0.2)
SCORE_MULTIPLIER = 0.2

# The motif sequences for each acr
MOTIF_LIST = '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/post_fimo/pairwise_inputs/acr_DAPv1.txt'

CHAINING_OUT = '/home/mwarr/Data/Chaining_one_DAPv1_clustered_global.tsv'

# Output file
DISTANCES_OUT = '/home/mwarr/Data/distances_DAPv1_clustered_global_shortest.tsv'

# Global or local
mode = 'global'

# For Global score alpha is how much we care about the longest, at 1 it is longest - len / longest
ALPHA = 0

# Here the beta is a linear multiplier for local scores: we are doing beta/chain score 
BETA = 5

############################################################

def make_max_length_dict(dict_file) :
    motifs_dict = defaultdict(int)

    with open(dict_file, 'r') as f:
        curr_acr = ""
        for line in f : 
            if 'ACR: ' in line :
                curr_acr = line.rstrip().split('ACR: ')[1]
            else :
                motifs_dict[curr_acr] += 1

    return motifs_dict

# alpha is how much we care about the longest, at 1 it is longest - len / longest

def get_all_distances(alpha, chaining_file, outfile, max_length_dict) :

    with open(outfile, 'w') as out:
        with open(chaining_file, 'r') as input :
            for line in tqdm(input, miniters= 1000000) :
                genome1, genome2, chain, _ = line.rstrip().split('\t')
                max1 = max_length_dict[genome1]
                max2 = max_length_dict[genome2]
                if max1 > max2:
                    longest = max1
                    shortest = max2
                else:
                    longest = max2
                    shortest = max1

                if shortest == 0 :
                    out.write(f"{genome1}\t{genome2}\t1\n")
                else :
                    out.write(f"{genome1}\t{genome2}\t{max(0, alpha * ((longest - (float(chain) * SCORE_MULTIPLIER))/longest) + (1 - alpha) * ((shortest - (float(chain) * SCORE_MULTIPLIER))/shortest))}\n")

# Here the beta is a linear multiplier: we are doing beta/chain score 
def get_all_distances_local(beta, chaining_file, outfile) :
    with open(outfile, 'w') as out:
        with open(chaining_file, 'r') as input :
            for line in tqdm(input, miniters= 1000000) :
                genome1, genome2, chain, _ = line.rstrip().split('\t')
                out.write(f"{genome1}\t{genome2}\t{beta / float(chain)}\n")

#########################################################

if mode == 'global' :
    max_length_dict = make_max_length_dict(MOTIF_LIST)
    # print(max_length_dict)

    get_all_distances(ALPHA, CHAINING_OUT, DISTANCES_OUT, max_length_dict)

else :

    get_all_distances_local(BETA, CHAINING_OUT, DISTANCES_OUT)