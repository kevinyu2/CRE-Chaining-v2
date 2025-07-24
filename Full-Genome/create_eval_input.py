# After aligning, cut it down as an input for evalulate blast
# Doing this so I don't have to run a bunch of times

SCORE_CUTOFF = 0
ALIGN_CUTOFF = -0.5

INPUT = '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/post_fimo/BLAST_hits/hits_DAPv1_scorev2l08_aligned.tsv'
OUTPUT = '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/post_fimo/BLAST_hits/hits_DAPv1_scorev2l08_cut.tsv'

##################################################################################

with open(INPUT, 'r') as input, open(OUTPUT, 'w') as output :
    output.write(f"{input.readline().rstrip()}\tFinal_cutoff: {SCORE_CUTOFF}\tFinal_align: {ALIGN_CUTOFF}\n")
    output.write(f"{input.readline()}")

    for line in input :
        line_arr = line.rstrip().split('\t')
        if float(line_arr[1]) >= SCORE_CUTOFF and float(line_arr[5]) <= ALIGN_CUTOFF :
            output.write(line)