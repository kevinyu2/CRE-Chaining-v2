from useful_functions import *
from frequencies import *

'''
This file contains functions for creating frequency files which exclude pairs with high alignment scores.

It also contains a function to create a score file which takes into account alignment (high alignment=lower score) 
and chain length (longer chain length = higher score). This was an old idea that didn't work very well.
'''

'''
Creates a new score file in the same format as the chaining and alignment files.
Normalizes each score with respect to the max chain/alignment score.
Flips the normalized alignment score so that better alignment = lower score.
Weights the alignment score at <align_frac> and the chain score at 1 - <align_frac>.
Outputs the two regions and score to a line of the file.
'''
def align_chain_scores(align_file, chain_file, output_file, align_frac):
    chain_dict = {}
    min_chain, max_chain = create_score_dict(chain_dict, chain_file)
    
    #Put the alignment score file into a list of tuples (region1, region2, score) and
    #find the max alignment score
    align_lst = []
    max_align = 0
    min_align = 10000
    with open(align_file) as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            align_score = float(line_arr[2])
            reg1 = min(line_arr[0], line_arr[1])
            reg2 = max(line_arr[0], line_arr[1])
            align_lst.append((reg1, reg2, align_score))
            if  align_score > max_align: 
                max_align = align_score
            if align_score < min_align:
                min_align = align_score
    
    chain_frac = 1 - align_frac
    with open(output_file, "w") as file:
        for reg1, reg2, align_score in align_lst:
            chain_score = chain_dict.get((reg1, reg2), 0)
            try:
                align_score_norm = align_frac * (1 - ((align_score - min_align) / (max_align - min_align)))
                chain_score_norm = chain_frac * ((chain_score - min_chain)/ (max_chain - min_chain))
            except ZeroDivisionError:
                print("Scores are all the same value")
                raise
            file.write(f"{reg1}\t{reg2}\t{align_score_norm + chain_score_norm}\n")

'''
Creates a dictionary for ACR test regions and for random regions in the form:
{region: [chain_score, chain_score, ...], ...}

Only includes chaining scores in the bottom <percent> percent if thresh=None.

If a threshold is passed in, <percent> will be ignored. All pairs with alignment scores above
<thresh> will be ignored.
'''
def exclude_high_align(align_file, chain_file, ref_set, percent, thresh=None):
    align_dict = {}
    create_score_dict(align_dict, align_file)
    if thresh == None:
        values_sorted = sorted(align_dict.values())
        thresh_index = int(((percent / float(100)) * len(values_sorted)))
        threshold = values_sorted[thresh_index]
    else:
        threshold = thresh
    return create_dicts(chain_file, ref_set, 3130, align_dict=align_dict, thresh=threshold)

'''
Takes in a normalized score, i.e. (score - min) / (max - min), <score_norm>
Outputs what fraction of the alignment scores have <score_norm> score or lower.
'''
def normalized_score_to_fraction(align_file, score_norm):
    align_dict = {}
    min_align, max_align = create_score_dict(align_dict, align_file)
    values_sorted = sorted(align_dict.values())
    thresh = (score_norm * (max_align - min_align)) + min_align # Turn normalized value back to regular value
    index = 0
    for ind, val in enumerate(values_sorted):
        if val > thresh:
            index = ind
            break
    print("thresh", thresh)
    print("index+1", index + 1)
    print("len_values", len(values_sorted))
    return (index + 1) / float(len(values_sorted))

'''
Driver for creating a score file which scores based on both alignment and chain length.

<type> should either be "global" or "local"
<type_short> should either be "glob" or "loc"
'''
def align_chain_driver(type, type_short, align_frac):
    base_dir = "/home/mwarr/Data/One_Genome/experiment2_10-90"
    align_file = f"{base_dir}/alignment/{type}/alignment_90-10_{type_short}.tsv"
    chain_file = f"{base_dir}/Chaining_one_acr_rand_10-90_{type_short}.tsv"
    output_file = f"{base_dir}/chain_and_align/align-{align_frac}_chain_{type_short}_scores.tsv"

    align_chain_scores(align_file, chain_file, output_file, align_frac)

'''
Driver for creating frequency files for chaining scores, excluding pairs which have alignment scores which are
not in the bottom <percent> percent.

<type> should either be "global" or "local"
<type_short> should either be "glob" or "loc"
'''
def exclude_high_align_driver(type, type_short, percent):
    ref_set = create_ref_set("/home/mwarr/Data/One_Genome/experiment2_10-90/seta_90.txt")
    base_dir = "/home/mwarr/Data/One_Genome/experiment2_10-90"
    dicts = exclude_high_align(f"{base_dir}/alignment/{type}/alignment_90-10_{type_short}.tsv", f"/home/mwarr/Data/Chaining_rand_DAPv1_clustered_{type}.tsv", ref_set, percent)
    for i in range(1, 6):
        lst_op = lambda lst: sorted(lst)[-i]
        output_score_freq(dicts[0], dicts[1], f"{base_dir}/chain_and_align/weighted/{type_short}_freq", lst_op, f"exclude_{round(percent, 2)}%_{i}-highest")


if __name__ == "__main__":
    exclude_high_align_driver("global", "glob", 99.99371685225528)
    exclude_high_align_driver("local", "loc", 99.99371685225528)



            
