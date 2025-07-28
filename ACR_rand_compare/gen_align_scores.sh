#!/bin/bash

#This script launches 100 jobs to run at once in order to get the alignment scores
#100 different output files will be generated inside <align_dir>
#rand_dir and 

rand_dir="/home/mwarr/Data/One_Genome/experiment2_10-90/alignment/rand_temp"
setb_dir="/home/mwarr/Data/One_Genome/experiment2_10-90/alignment/setb_temp" 
seta_file="/home/mwarr/Data/One_Genome/experiment2_10-90/alignment/seta_90.fa"
align_dir="/home/mwarr/Data/One_Genome/experiment2_10-90/alignment/local/align_temp"
script_dir="./get_alignment_scores.py"

for i in {0..49};
do 
    nohup python $script_dir ${rand_dir}/temp_${i} ${align_dir}/temp_rand_${i} &
    nohup python $script_dir ${setb_dir}/temp_${i} ${align_dir}/temp_${i} &
done;




