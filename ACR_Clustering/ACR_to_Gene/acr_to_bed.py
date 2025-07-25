from pathlib import Path
import sys


def acr_to_bed(output_file) :

    # Get all _xtreme folders
    input_file = Path('/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/dACR.fa')
    base_names = []

    with open(input_file, 'r') as f :
        for line in f :
            if line[0] == '>' :
                base_names.append(line.rstrip().split('>')[1])
    # print(base_names)
    with open(output_file, 'w') as out:
        for name in base_names :
            stop = int(name.split('to')[1])
            start = int(name.split('to')[0].split('_')[1])
            
            # Swap if wrong order
            if stop < start :
                temp = stop
                stop = start
                start = temp
                print("swap")
            chromosome = name.split('_')[0]
            out.write(f"{chromosome}\t{start}\t{stop}\t{name}\n")

acr_to_bed("/home/kyu/CRE-chaining/Single-Genome/acr.bed")