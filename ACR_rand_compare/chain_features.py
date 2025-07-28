from matplotlib import pyplot as plt #type: ignore
import sys
from useful_functions import *
from frequencies import *
from collections import defaultdict

'''
This file contains functions to look at features of the chains, such as the number of anchors,
the length of the reference ACR, and the number of top chain scores.
'''


'''
Takes in a chaining file with all random regions and ACRs pairwise chained. Also takes in the
set of "reference" ACRs, i.e the ACRs which should be chained against the random regions and other
ACRs.

Returns two dictionaries with the chain scores, number of anchors, and reference ACR lengths
(for exact contents and structure, see in-line comment). The dictionaries are returned as a tuple: 
(ACR_dict, rand_dict), where ACR_dict has all the non-reference ACRs as keys, and rand_dict has
all the random regions as keys.

Note: Unlike the similar function in rand_vs_acr.py, this function does not account for chain
lengths of 0
'''
def create_dicts_mult_feat(chaining_file, reference_set):
    #For every non-reference ACR, this dictionary contains a list of anchor numbers and lengths for that ACR
    #chained with every reference ACR
    #{non-ref-acr: [(score, anchor_num, ref-acr-len), (.., ..)...]}
    ACR_dict = defaultdict(list)
    #For every random sequence, this dictionary contains a list of anchor numbers and lengths of that random sequence
    #chained with every reference ACR
    #{random_reg: [(score, anchor_num, ref-acr-len), (.., ..)...]}
    rand_dict = defaultdict(list)

    with open(chaining_file, "r") as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 4:
                continue
            item1 = line_arr[0]
            item2 = line_arr[1]
            score = float(line_arr[2])
            anchor_num = int(line_arr[3])
            item1_len = get_length(item1)
            item2_len = get_length(item2)

            if item1[-4:] == "rand":
                if item2 in reference_set:
                    rand_dict[item1].append((score, anchor_num, item2_len))
            elif item1 in reference_set:
                if item2[-4:] == "rand":
                    rand_dict[item2].append((score, anchor_num, item1_len))
                elif item2 not in reference_set:
                    ACR_dict[item2].append((score, anchor_num, item1_len))
            else:
                if item2 in reference_set:
                    ACR_dict[item1].append((score, anchor_num, item2_len))
    
    return(ACR_dict, rand_dict)

'''
A helper function. Takes in a multi-feature list (i.e. a list of tuples
where the first item is the chain length, second item is the number of anchors, and third item
is the length of the reference region.) 

Outputs a dictionary indicating the 5th-longest chain length, the number of 5th-longest chain
lengths, the minimum number of anchors between the non-ref region and a 5th-longest chain, and the 
minimum ref-length between a non-ref region and a ref-region.
'''
def get_features_combine_score(lst):
    chain_len_5 = sorted(lst)[-5][0]
    lst_5th = [] #A list of all tuples where score is the 5th highest score
    for item in lst:
        if item[0] == chain_len_5:
            lst_5th.append(item)
    chain_num = len(lst_5th)
    anchor_lst = [item[1] for item in lst_5th]
    ref_len_lst = [item[2] for item in lst_5th]
    anchor_num = min(anchor_lst)
    ref_len_num = min(ref_len_lst)
    return {"chain_len" : chain_len_5, "chain_num": chain_num, "anchor_num": anchor_num,
            "ref_len_num": ref_len_num}

'''
The following four functions are all wrapper functions to create a function for output_score_freq.
'''
def list_op_anchor(chain_score):
    '''
    Function to be passed into output_score_freq. Returns a list of all the anchor numbers
    of the regions with a score of <chain_score> if and only if the 5th highest score is <chain_score>. 
    Otherwise returns 0.
    '''
    def anchor_num(lst):
        score_lst = [item[0] for item in lst]
        if sorted(score_lst)[-5] == chain_score: #5th highest score
            anchor_lst = []
            for item in lst:
                if item[0] == chain_score:
                    anchor_lst.append(item[1])
            return anchor_lst
        return 0
    return anchor_num

def list_op_len_ref(chain_score):
    '''
    Function to be passed into output_score_freq. Returns a list of all the reference region lengths
    of the regions with a score of <chain_score> if and only if the 5th highest score is <chain_score>.
    Otherwise returns 0.
    '''
    def len_ref(lst):
        score_lst = [item[0] for item in lst]
        if sorted(score_lst)[-5] == chain_score: #5th highest score is 10
            len_lst = []
            for item in lst:
                if item[0] == chain_score:
                    len_lst.append(item[2])
            return len_lst
        return 0
    return len_ref

def list_op_count(num_highest):
    '''
    Function to pass into output_score_freq. Takes in a multi-feature list (i.e. a list of tuples
    where the first item is the chain length, second item is the number of anchors, and third item
    is the length of the reference region.) 

    Outputs the number of <num_highest>th-highest scores.
    '''
    def list_count(lst):
        lst_scores = [item[0] for item in lst]
        return lst_scores.count(sorted(lst_scores)[-num_highest])
    
    return list_count

def list_op_combine(chain_max, chain_num_max, anchor_max, ref_num_max, chain_frac):
    '''
    Outputs a score from the multi-feature list. The score is normalized based on the maximums and
    the <chain_frac>. <chain_frac> is the weight that the chain length is given. The remaining
    features are weighted evenly. The score is between 0 and 1.
    '''
    def list_op(lst):
        dict = get_features_combine_score(lst)

        chain_len_5 = float(dict["chain_len"]) / chain_max #higher = better
        chain_num = 1 - (float(dict["chain_num"]) / chain_num_max) #lower = better
        anchor_num = 1 - (float(dict["anchor_num"]) / anchor_max) #lower = better
        ref_len_num = 1 - (float(dict["ref_len_num"]) / ref_num_max) #lower = better

        other_frac = (1-chain_frac) / 3
        return chain_frac*chain_len_5 + other_frac*chain_num + other_frac*anchor_num + other_frac*ref_len_num
    return list_op

'''
Returns a dictionary indicating the max 5th-longest chain length, the max number of 5th-longest
chain lengths, the max min number of anchors, and the max min ref-length. These numbers are used to
create normalized scores between 0 and 1 (a ratio between the number and the max number).
'''
def get_max_dict(ACR_dict, rand_dict):
    max_dict = {"chain_len" : 0, "chain_num": 0, "anchor_num": 0, "ref_len_num": 0}
    for lst in ACR_dict.values():
        dict = get_features_combine_score(lst)
        for feature, val in dict.items():
            if max_dict[feature] < val:
                max_dict[feature] = val
    for lst in rand_dict.values():
        dict = get_features_combine_score(lst)
        for feature, val in dict.items():
            if max_dict[feature] < val:  
                max_dict[feature] = val
    return max_dict

'''
Driver for outputting the frequencies for all the different chain features.
<base> is the base directory which contains subdirectories for the frequency files to
be output to.
<chain_path> is the path to the chaining file.
'''
def driver_frequencies_all(base, chain_path):

    ref_set = create_ref_set("/home/mwarr/Data/One_Genome/experiment2_10-90/seta_90.txt")
    dicts = create_dicts_mult_feat(chain_path, ref_set)

    #output frequencies for number of top scores.
    for i in range(1, 6):
        output_score_freq(dicts[0], dicts[1], f"{base}/counts", list_op_count(i), f"count_{i}")
    
    #output frequencies for chain features with a fixed chain length
    for i in range(2, 20):
        output_score_freq(dicts[0], dicts[1], f"{base}/anchor_num", list_op_anchor(i), f"anchor_num_{i}")
        output_score_freq(dicts[0], dicts[1], f"{base}/ref_len", list_op_len_ref(i), f"ref_len_{i}")
    
    #output frequencies for combined score
    max_dict = get_max_dict(dicts[0], dicts[1])
    for frac in [.25, .4, .79]:
        function = list_op_combine(max_dict["chain_len"], max_dict["chain_num"], max_dict["anchor_num"], max_dict["ref_len_num"], frac)
        output_score_freq(dicts[0], dicts[1], f"{base}/combine/chain_{frac}", function, f"combine_{frac}")

if __name__ == "__main__":
    base = "/home/mwarr/Data/One_Genome/other_features/local"
    chain_path = "/home/mwarr/Data/One_Genome/experiment2_10-90/Chaining_one_acr_rand_10-90_loc.tsv"
    driver_frequencies_all(base, chain_path)