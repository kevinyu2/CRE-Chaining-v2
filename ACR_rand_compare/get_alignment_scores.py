from Bio import Align #type: ignore
import random
import time
import math
import sys
from pathlib import Path
from useful_functions import *

'''
This script can be used to generate a file with all the alignment scores between random regions and reference
ACRs, as well as non-reference ACRs with reference ACRs.

The test (non-ref) ACR set and the random regions should be split into many different fasta files.
Then, use the gen_align_scores.sh bash script to launch the jobs to get the alignment scores between these files
and the fasta file with reference ACR regions. The scores will be output to many different files, which can then be 
concatenated. This allows the 31k arabidopsis ACR alignment scores to be found in a couple of hours.

Time_test can be used to determine how long this will take to run. This will help in determining how many
files to split your fasta files into.
'''

'''
Prints out the time it takes for regions of size <region_len> to be aligned. <align_type> should either be
'global' or 'local'.
'''
def time_test(region_len, align_type):
    bases = ['A', 'C', 'T', 'G']
    target = ""
    query = ""
    for i in range(region_len):
        target += bases[random.randint(0, 3)]
        query += bases[random.randint(0, 3)]
    
    print("Beginning alignment")
    aligner = Align.PairwiseAligner(scoring="blastn")
    aligner.mode = align_type
    
    start_time = time.time()
    score = aligner.score(target, query)
    print(time.time() - start_time)
    print(score)

'''
Takes in a fasta file (input_file) and outputs <num_files> files to <output_dir>.
The output files together will contain the same sequences as the original fasta file.
<num_sequences> is the number of sequences in <input_file>.
'''
def split_fasta(input_file, output_dir, num_files, num_sequences):
    seq_per_file = math.ceil(num_sequences / num_files)
    
    with open(input_file) as input:
        id = input.readline()
        for i in range(num_files):
            with open(f"{output_dir}/temp_{i}", "w") as output:
                for _ in range(seq_per_file):
                    output.write(id)
                    line = input.readline()
                    while line and line != "\n" and line[0] != ">":
                        output.write(line)
                        line = input.readline()
                    if line == "\n" or not line:
                        break
                    id = line

'''
Takes in two fasta files (input_file and ref_file) and outputs the alignment score between
every sequence in <input_file> with every sequence in <ref_file> to <output_file>.
'''
def output_align(input_file, ref_file, output_file):
    aligner = Align.PairwiseAligner(scoring="blastn")
    aligner.mode = "global"
    
    #read in files and create dictionaries
    ref_dict = create_genome_dict(ref_file)
    input_dict = create_genome_dict(input_file)
    
    output = []
    
    #pairwise align
    for id2, seq2 in input_dict.items():
        for id1, seq1 in ref_dict.items():
            score = aligner.score(seq1, seq2)
            output.append((id1, id2, score))

    #output to file
    with open(output_file, "w") as out_file:
        for item in output:
            out_file.write(f"{item[0]}\t{item[1]}\t{item[2]}\n")


if __name__ == "__main__":
    output_align(sys.argv[1], sys.argv[2], sys.argv[3])
