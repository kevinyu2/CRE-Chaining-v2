from collections import defaultdict
import os
import time
from useful_functions import *
from pathlib import Path

'''
This file contains functions to create frequency files and histograms of scores for random regions
and ACRs.
'''

'''
Takes in a file with a list of ACRs and creates a set. Only includes ACRs with lengths between
<low> and <high>
'''
def create_ref_set_filter(input_path, low, high):
    ref_set = set()
    with open(input_path, "r") as file:
        for line in file:
            line = line.strip()
            if line != "":
                length = get_length(line)
                if length >= low and length <= high:
                    ref_set.add(line)
    return ref_set

'''
Takes in a chaining file with all random regions and ACRs pairwise chained. Also takes in the
set of "reference" ACRs, i.e the ACRs which should be chained against the random regions and other
ACRs.
Non_ref_num is the number of non-reference and random regions there should be.

To restrict the dictionaries to only ACRs of a certain length, pass in parameters for low (lower bound
on the length) and high (upper bound on the length)

To exclude high alignment pairs, pass in a dictionary in the form {(reg1, reg2) : align_score, ...} and a threshold.
Pairs which have alignment scores strictly greater than <thresh> will be ignored.

Returns two dictionaries (ACR_dict, rand_dict) which contain lists of scores for each non-ref ACR
and random region, respectively.
'''
def create_dicts(chaining_file, reference_set, non_ref_num, low=None, high=None, align_dict=None, thresh=None):
    #For every non-reference ACR, this dictionary contains a list of scores of that ACR
    #chained with every reference ACR
    ACR_dict = defaultdict(list)
    #For every random sequence, this dictionary contains a list of scores of that random sequence
    #chained with every reference ACR
    rand_dict = defaultdict(list)

    with open(chaining_file, "r") as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 3:
                continue
            item1 = min(line_arr[0], line_arr[1])
            item2 = max(line_arr[0], line_arr[1])
            score = float(line_arr[2])

            #Skip if alignment score is too high
            if align_dict != None:
                try:
                    if align_dict[(item1, item2)] > thresh:
                        continue
                except KeyError: 
                #not in align_dict, so either both random, both reference, both non-ref,
                #or one random and one non-ref
                    continue

            if low != None: 
                len1 = get_length(item1)
                len2 = get_length(item2)
                #skip if outside the bounds for the length
                if len1 < low or len1 > high:
                    continue
                if len2 < low or len2 > high:
                    continue

            if item1[-4:] == "rand":
                if item2 in reference_set:
                    rand_dict[item1].append(score)
            elif item1 in reference_set:
                if item2[-4:] == "rand":
                    rand_dict[item2].append(score)
                elif item2 not in reference_set:
                    ACR_dict[item2].append(score)
            else:
                if item2 in reference_set:
                    ACR_dict[item1].append(score)
        
    #Account for any reference ACRs which have chaining scores of 0 with some ACR / random region
    for lst in ACR_dict.values():
        missing_values = len(reference_set) - len(lst)
        lst.extend([0 for _ in range(missing_values)])
        assert(len(lst) == len(reference_set))
    for lst in rand_dict.values():
        missing_values = len(reference_set) - len(lst)
        lst.extend([0 for _ in range(missing_values)])
        assert(len(lst) == len(reference_set))
    
    #Account for any non-reference and random regions with chaining scores of 0 with every ref ACR
    for i in range(len(ACR_dict), non_ref_num):
        ACR_dict[f"unknown_{i}"].extend([0 for i in range(len(reference_set))])
    assert(len(ACR_dict) == non_ref_num)
    
    for i in range(len(rand_dict), non_ref_num):
        rand_dict[f"unknown_rand_{i}"].extend([0 for i in range(len(reference_set))])
    assert(len(rand_dict) == non_ref_num)
    
    return(ACR_dict, rand_dict)

'''
Outputs score frequencies to a file.

<ACR_dict> and <rand_dict> should contain non-ref ACR regions and random regions as keys, respectively.
The values should be lists of chaining scores for the given region. These dictionaries are output by the
create_dicts() function.
'''
def output_all_score_freq(ACR_dict, rand_dict, output_dir):
    ACR = {}
    rand = {}
    for score_list in ACR_dict.values():
        for item in score_list:
            ACR[item] = ACR.get(item, 0) + 1
    for score_list in rand_dict.values():
        for item in score_list:
            rand[item] = rand.get(item, 0) + 1
    
    with open(f"{output_dir}/ACR_vs_ACR_all_freq.tsv", "w") as file:
        for score in ACR:
            file.write(f"{score}\t{ACR[score]}\n")
    
    with open(f"{output_dir}/rand_vs_ACR_all_freq.tsv", "w") as file:
        for score in rand:
            file.write(f"{score}\t{rand[score]}\n")

'''
Does some operation on each non-reference ACR / random region to get some number (e.g. max, average,
min, etc.). Outputs the frequencies of this number to a file.  

<ACR_dict> and <rand_dict> should contain non-ref ACR regions and random regions as keys, respectively.
The values should be lists of chaining scores for the given region. These dictionaries are output by the
create_dicts() function.

<list_op> is the operation which is performed on each list. It should be a function which takes in a list
and returns either a number of a list of numbers.

<op_name> is the name of the operation for the output file.
'''
def output_score_freq(ACR_dict, rand_dict, output_dir, list_op, op_name):
    ACR = {}
    rand = {}
    for score_list in ACR_dict.values():
        value = list_op(score_list)
        if type(value) == list:
            for val in value:
                ACR[val] = ACR.get(val, 0) + 1
        else:
            ACR[value] = ACR.get(value, 0) + 1
    for score_list in rand_dict.values():
        value = list_op(score_list)
        if type(value) == list:
            for val in value:
                rand[val] = rand.get(val, 0) + 1
        else:
            rand[value] = rand.get(value, 0) + 1
    
    with open(f"{output_dir}/ACR_vs_ACR_{op_name}_freq.tsv", "w") as file:
        for score in ACR:
            file.write(f"{score}\t{ACR[score]}\n")
    
    with open(f"{output_dir}/rand_vs_ACR_{op_name}_freq.tsv", "w") as file:
        for score in rand:
            file.write(f"{score}\t{rand[score]}\n")

'''
Takes in the dictionary of ACR scores and random region scores (output by create_dicts() function).
Outputs the max and min values over both dictionaries.
'''
def get_max_min(ACR_dict, rand_dict):
    max = None
    min = None
    for lst in ACR_dict.values():
        for num in lst:
            if max == None or num > max:
                max = num
            if min == None or num < min:
                min = num
    for lst in rand_dict.values():
        for num in lst:
            if num > max:
                max = num
            if num < min:
                min = num
    return {"max": max, "min": min}

'''
Takes in a directory <input_dir> and outputs the lowest alignment score in the frequency files
in that directory. Useful for determining what the adjustment should be to make all scores non-negative.
'''
def get_lowest_align_scores(input_dir):
    min = 0
    for file_path in Path(input_dir).glob("*freq.tsv"):
        with open(file_path) as file:
            for line in file:
                line_arr = line.split("\t")
                if len(line_arr) < 2:
                    continue
                score = float(line_arr[0])
                if score < min:
                    min = score
    return min

'''
Adds <adjustment> to every score in the frequency files in <input_dir>.
Useful for adjusting all scores to be non-negative values.
'''
def adjust_scores(adjustment, input_dir):
    for file_path in Path(input_dir).glob("*freq.tsv"):
        with open(f"{file_path.parent}/{file_path.stem}_2.tsv", "w") as output:
            print(f"{file_path.parent}/{file_path.stem}_2.tsv")
            with open(file_path) as file:
                for line in file:
                    line_arr = line.split("\t")
                    if len(line_arr) < 2:
                        continue
                    output.write(f"{float(line_arr[0]) + adjustment}\t{line_arr[1]}")

'''
Function to be passed into output_score_freq as <list_op>.
Averages the top ten values in <lst>
'''
def list_op_top_avg(lst):
    lst.sort(reverse=True)
    sum = 0
    for i in range(10):
        sum += lst[i]
    return sum / 10

'''
Function to create function to be passed into output_score_freq as <list_op>
'''
def list_op_top_norm(max, min, num):
    '''
    Finds the <num>th largest value in the list and returns the value normalized by the
    range of possible values.
    '''
    def top_norm(lst):
        value = float(sorted(lst)[-(num)])
        return (value - min) / (max - min)
    return top_norm

'''
Reads in score file and reference set file, outputs score frequency files, and
creates histograms.
'''
def driver_all_summary():
    ref_set = create_ref_set("/home/mwarr/Data/seta.txt")
    dicts = create_dicts("/home/mwarr/Data/alignment/global/alignment_scores_glob.tsv", ref_set, 15648)
    print(f"Non reference ACRs: {len(dicts[0])}")
    print(f"Random regions: {len(dicts[1])}")
    output_all_score_freq(dicts[0], dicts[1], "/home/mwarr/Data/alignment/global")
    output_score_freq(dicts[0], dicts[1], "/home/mwarr/Data/alignment/global", max, "max")
    output_score_freq(dicts[0], dicts[1], "/home/mwarr/Data/alignment/global", lambda lst: sum(lst) / float(len(lst)), "avg")
    create_histograms("/home/mwarr/Data/alignment/global/ACR_vs_ACR_all_freq.tsv", "/home/mwarr/Data/alignment/global/rand_vs_ACR_all_freq.tsv", "All Global Alignment Scores", "Score", 11000, 10000)
    create_histograms("/home/mwarr/Data/alignment/global/ACR_vs_ACR_avg_freq.tsv", "/home/mwarr/Data/alignment/global/rand_vs_ACR_avg_freq.tsv", "Average Global Alignment Score", "Average Score", 11000, 10000)
    create_histograms("/home/mwarr/Data/alignment/global/ACR_vs_ACR_max_freq.tsv", "/home/mwarr/Data/alignment/global/rand_vs_ACR_max_freq.tsv", "Max Global Alignment Score", "Max Score", 11000, 10000)

'''
Reads in score file and reference set file, outputs score frequency files, and
creates histograms. Repeats 14 times, each time filtering to only include regions with
lengths within a 100bp interval.
'''
def driver_filter_summary():
    for i in range(650, 2050, 100):
        start_time = time.time()
        dir = f"/home/mwarr/Data/One_Genome/local_freq_{i}-{i+100}"
        os.mkdir(dir)
        ref_set = create_ref_set_filter("/home/mwarr/Data/seta.txt", i, i+100)
        non_ref_set = create_ref_set_filter("/home/mwarr/Data/setb.txt", i, i+100)
        dicts = create_dicts("/home/mwarr/Data/One_Genome/Chaining_one_acr_rand_loc.tsv", ref_set, len(non_ref_set), i, i+100)
        output_all_score_freq(dicts[0], dicts[1], dir)
        output_score_freq(dicts[0], dicts[1], dir, max, "max")
        output_score_freq(dicts[0], dicts[1], dir, lambda lst: sum(lst) / float(len(lst)), "avg")
        create_histograms(f"{dir}/ACR_vs_ACR_all_freq.tsv", f"{dir}/rand_vs_ACR_all_freq.tsv", f"All Local Chaining Scores Lengths {i}-{i+100}", "Score")
        create_histograms(f"{dir}/ACR_vs_ACR_avg_freq.tsv", f"{dir}/rand_vs_ACR_avg_freq.tsv", f"Average Local Chaining Score Lengths {i}-{i+100}", "Average Score")
        create_histograms(f"{dir}/ACR_vs_ACR_max_freq.tsv", f"{dir}/rand_vs_ACR_max_freq.tsv", f"Highest Local Chaining Score Lengths {i}-{i+100}", "Max Score")
        print(f"Finished {i}-{i+100} in {time.time() - start_time} seconds", flush=True)

'''
Outputs frequency files for top 5 score ranks.
<type> should be 'global' or 'local'
<type_short> should be 'glob' or 'loc'
'''
def driver_top(type, type_short):
    ref_set = create_ref_set("/home/mwarr/Data/One_Genome/experiment2_10-90/seta_90.txt")
    dicts = create_dicts(f"/home/mwarr/Data/Chaining_rand_DAPv1_clustered_{type}.tsv", ref_set, 3130)
    out_dir = f"/home/mwarr/Data/One_Genome/experiment2_10-90/{type_short}_freq_weighted"
    for i in range(1, 6):
        output_score_freq(dicts[0], dicts[1], out_dir, lambda lst: sorted(lst)[-(i)], f"score_{i}-highest")

'''
Outputs frequency files for top 5 score ranks. Normalizes the scores.
'''
def driver_top_normal():
    ref_set = create_ref_set("/home/mwarr/Data/One_Genome/experiment2_10-90/seta_90.txt")
    dicts = create_dicts("/home/mwarr/Data/Chaining_rand_DAPv1_clustered_global.tsv", ref_set, 3130)
    out_dir = "/home/mwarr/Data/One_Genome/experiment2_10-90/glob_freq_weighted"
    max_min = get_max_min(dicts[0], dicts[1])
    for i in range(1, 6):
        output_score_freq(dicts[0], dicts[1], out_dir, list_op_top_norm(max_min["max"], max_min["min"], i), f"score_{i}-highest")

if __name__ == "__main__":
    driver_top("local", "loc")
    driver_top("global", "glob")