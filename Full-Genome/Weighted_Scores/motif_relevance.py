
from collections import defaultdict
from pathlib import Path
import pandas as pd

# Set to none if you don't want to cut off
# ACR_SET = '/home/mwarr/Data/seta_half.txt'
ACR_SET = None

chrs_to_use = ['Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5']

FULL_SEQUENCE_FOLDER = '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/post_fimo/BLAST_inputs/v2ArabidopsisDAPv1_fimo/'

ACR_FILE = '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/post_fimo/BLAST_inputs/v2acr_motif_sequences_ArabidopsisDAPv1.txt'

OUT_FILE = '/home/mwarr/Data/motif_relevance/DAPv1_relevance_seta.tsv'

###########################################################################################
# Parse input files

if ACR_SET != None :
    acrs_to_use = set()
    with open(ACR_SET, 'r') as f :
        for line in f :
            acrs_to_use.add(line.rstrip())


def dict_with_keys_zero():
    return {key: 0 for key in chrs_to_use}

# {motif: {chr1 : num, chr2 : num...}}
full_seq_chr_counts = defaultdict(dict_with_keys_zero)
# {motif: num}
full_seq_total_counts = defaultdict(int)

total_full = 0

for file_path in Path(FULL_SEQUENCE_FOLDER).iterdir():
    chromosome = str(file_path).rstrip().split('/')[-1]
    with open(file_path, 'r') as motif_file :
        for line in motif_file :
            if chromosome in chrs_to_use :
                if 'BREAK' not in line.split('\t')[0] :

                    full_seq_chr_counts[line.split('\t')[0]][chromosome] += 1
                    full_seq_total_counts[line.split('\t')[0]] += 1

                    total_full += 1

# {motif: {chr1 : num, chr2 : num...}}
acr_chr_counts = defaultdict(dict_with_keys_zero)
# {motif: num}
acr_total_counts = defaultdict(int)

total_acr = 0

with open(ACR_FILE, 'r') as acr_file : 

    curr_acr = ""
    for line in acr_file :
        # End of one
        if 'ACR: ' in line:
            # RESET
            curr_acr = line.rstrip().split('ACR: ')[1]

        elif ACR_SET != None :
            if curr_acr in acrs_to_use :
                if curr_acr.split('_')[0] in chrs_to_use :
                    acr_chr_counts[line.rstrip()][curr_acr.split('_')[0]] += 1
                    acr_total_counts[line.rstrip()] += 1
                    total_acr += 1
        else :
            if curr_acr.split('_')[0] in chrs_to_use :

                acr_chr_counts[line.rstrip()][curr_acr.split('_')[0]] += 1
                acr_total_counts[line.rstrip()] += 1

                total_acr += 1

################################################################################
# Output

def z_score(x1, n1, x2, n2):
    p1 = x1 / n1
    p2 = x2 / n2
    p = (x1 + x2) / (n1 + n2)
    numerator = p1 - p2
    denominator = (p * (1 - p) * (1/n1 + 1/n2)) ** 0.5
    return numerator / denominator if denominator != 0 else 0

# print(acr_chr_counts['cluster_14']['Chr1'] + acr_chr_counts['cluster_14']['Chr2'] + acr_chr_counts['cluster_14']['Chr3'] + acr_chr_counts['cluster_14']['Chr4'] + acr_chr_counts['cluster_14']['Chr5'])

# print(full_seq_chr_counts['cluster_14']['Chr1'])
with open(OUT_FILE, 'w') as out :
    out.write('Motif\t')
    # for chr in chrs_to_use :
    #     out.write(f'ACR_{chr}\tFull_{chr}\t')
    out.write('Total_ACR\tTotal_Full\tACR_percent\tFull_percent\tZScore\n')


    sorted_acrs = sorted(
        full_seq_chr_counts.keys(),
        key = lambda acr: z_score(acr_total_counts[acr], total_acr, full_seq_total_counts[acr], total_full),
        # key=lambda acr: (acr_total_counts.get(acr, 0) / total_acr) - (full_seq_total_counts.get(acr, 0) / total_full),
        reverse=True  # Optional: sort descending by difference
    )
    for acr in sorted_acrs :
        out.write(f"{acr}\t")
        # for chr in chrs_to_use :
        #     out.write(f"{acr_chr_counts[acr][chr]}\t{full_seq_chr_counts[acr][chr]}\t")

        out.write(f"{acr_total_counts[acr]}\t{full_seq_total_counts[acr]}\t")
        out.write(f"{acr_total_counts[acr]/total_acr}\t{full_seq_total_counts[acr]/total_full}\t")
        out.write(f"{z_score(acr_total_counts[acr], total_acr, full_seq_total_counts[acr], total_full)}\n")