from matplotlib import pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import random

'''
Generates box and whisker plots from frequency files. 

Frequency file format:
Each line should have the value, followed by a tab, followed
by the frequency of that value.
'''


'''
Helper function.
Takes in a frequency file and creates a list of all the values. If ignore_zeros is True, zeros will not be included
in the list. If max_value is True, any value greater than <max_value> will not be included
in the list.
'''
def format_data(input_file, ignore_zeros=True, max_value=None):
    data = []
    with open(input_file) as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            value = float(line_arr[0])
            freq = int(line_arr[1])
            if ignore_zeros and value == 0:
                continue
            if max_value != None and value > max_value:
                continue
            extension = [value for _ in range(freq)]
            data.extend(extension)
    return data

'''
Generates a box plot and saves it to <out_path>.

<data> should be a list of lists, where each list contains the values in a distribution. One
box will be created for each list.
<xlabels> should be a list of labels for the x-ticks.
<xtitle> and <ytitle> are the labels for the axes.
<ylim> is the range of the y-axis (optional)
<colors> should either be "two" (alternate blue and yellow) or a list of colors for pair of boxes.
If colors is None, every pair of boxes will be a random color.
'''
def create_box_plot(data, out_path, xlabels, title, xtitle, ytitle, ylim=None, colors=None):
    plt.figure(figsize=(16, 8))
    bplot = plt.boxplot(data, positions=range(1, len(data) + 1), sym="", patch_artist=True, 
                        whiskerprops=dict(linewidth=4), medianprops=dict(linewidth=4, color="#B0413E"),
                        boxprops=dict(linewidth=4), capprops=dict(linewidth=4))
    color_rand = (random.random(), random.random(), random.random())
    for ind, box in enumerate(bplot["boxes"]):
        if colors == None:
            box.set_facecolor(color_rand)
            if ind % 2 == 1: #regenerate every 2 boxes
                color_rand = (random.random(), random.random(), random.random())
        elif colors == "two":
            if ind % 2 == 0:
                box.set_facecolor("#00529b")
            else:
                box.set_facecolor("#FFB302")
        else:
            box.set_facecolor(colors[ind // 2])
    plt.title(title, pad=15, fontsize=25)
    plt.xlabel(xtitle, labelpad=20, fontsize=25)
    plt.ylabel(ytitle, labelpad=15, fontsize=25)
    if len(xlabels) == len(data):
        xtick_pos = range(1, len(data) + 1)
    else:
        xtick_pos = np.linspace(1.5, len(data) - .5, num=int(len(data) / 2))
        plt.legend(fontsize=25 , handles=[Patch(facecolor="#00529b", label="ACR"), Patch(facecolor="#FFB302", label="Random")])

    plt.xticks(xtick_pos, xlabels, fontsize=25)
    plt.yticks(fontsize=25)
    if ylim != None:
        plt.ylim(ylim[0], ylim[1])
    plt.tight_layout(pad=2)
    plt.savefig(out_path, dpi=300)

'''
Reads in the frequency files from <base_dir> which contain frequencies of some feature of
top 5 scores and generates a box plot.

There should be 10 files in <base_dir> in the following format:
<base_dir>/<ACR|rand>_vs_ACR_<op>_<num>-highest_freq.tsv

There should be 5 files with 'rand' and 'num'-highest score, where 'num' is between 1 and 5.
There should be 5 files with 'ACR' and 'num'-highest score, where 'num' is between 1 and 5.

<filename> is the output file name.

<max-value> is the maximum value which should be included in the data (optional)
'''
def top_5_all(base_dir, op, title, xtitle, ytitle, filename, max_value=None, ylim=None):
    data = []
    #labels = ["Max,\nACR", "Max,\nRandom", "2nd Highest,\nACR", "2nd Highest,\nRandom", "3rd Highest,\nACR", "3rd Highest,\nRandom"]
    labels=[]
    for i in range(1, 6):
        acr_data = format_data(f"{base_dir}/ACR_vs_ACR_{op}_{i}-highest_freq.tsv", max_value)
        rand_data = format_data(f"{base_dir}/rand_vs_ACR_{op}_{i}-highest_freq.tsv", max_value)
        data.append(acr_data)
        data.append(rand_data)
        labels.append(f"Rank {i}")
        # if i > 3:
        #     labels.append(f"{i}th Highest,\nACR")
        #     labels.append(f"{i}th Highest,\nrandom")
    colors = ["#00529b", "#FFB302", "#B0413E", "#4B4B4B", "#629DB7"]
    create_box_plot(data, labels, title, xtitle, ytitle, filename, colors="two", ylim=ylim)

'''
Reads in the frequency files from <base_dir> which contain the frequencies of some
feature with a fixed chain length and generates a box and whisker plot.

There should be 2 files for 18 fixed chain lengths in <base_dir> in the following format:
<base_dir>/<ACR|rand>_<op>_<chain-len>_freq.tsv

One file should contain 'rand' and another should contain 'ACR'. Chain-len should range
from 2-19.

<filename> is the output file name.

<max-value> is the maximum value which should be included in the data (optional)
'''
def fixed_chain_all(base_dir, op, title, xtitle, ytitle, filename, max_value=None):

    data = []
    labels = []
    for i in range(2, 20):
        acr_data = format_data(f"{base_dir}/ACR_vs_ACR_{op}_{i}_freq.tsv", max_value)
        rand_data = format_data(f"{base_dir}/rand_vs_ACR_{op}_{i}_freq.tsv", max_value)
        data.append(acr_data)
        labels.append(f"{i}, \nACR")
        data.append(rand_data)
        labels.append(f"{i}, \nrand")
    create_box_plot(data, labels, title, xtitle, ytitle, filename)

'''
Creates a box plot and a histogram for a single feature. The box plot will have two distributions: one
for the random regions and one for the ACRs. The histogram will have two side-by-side subplots: one for the
random regions and one for the ACRs.
'''
def single_compare(freq_file_rand, freq_file_acr, title, out_box, out_hist):
    data = []
    data.append(format_data(freq_file_acr))
    data.append(format_data(freq_file_rand))
    labels = ["ACR", "Random"]
    create_box_plot(data, labels, title, "Type of Region", "Score", out_box)
    
    plt.clf()
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True, sharex=True)
    axes[0].hist(data[0], bins=100)
    axes[0].set_title("ACRs chained with ACRs")
    axes[0].set_ylabel("Frequency")
    
    axes[1].hist(data[1], bins=100)
    axes[1].set_title("ACRs chained with random regions")
    axes[1].set_ylabel("Frequency")

    plt.suptitle(title)
    plt.tight_layout()
    plt.savefig(out_hist)

'''
Generates box and whisker plots for the top chaining scores with weighted scoring and high alignment
excluded from the data.
'''
def driver_weighted_exclude(type, type_short, percent):
    base_dir = f"/home/mwarr/Data/One_Genome/experiment2_10-90/chain_and_align/weighted/{type_short}_freq"
    op = f"exclude_{percent}%"
    title = f"Distribution of Top {type.capitalize()} Chain Scores, Weighted Scoring, High Alignment Excluded"
    xtitle = "Score Rank and Region Type"
    ytitle = "Score"
    filename = f"box_weighted_exclude_{type_short}.png"
    top_5_all(base_dir, op, title, xtitle, ytitle, filename)

'''
Generates box and whisker plots for the top alignment scores
'''
def driver_align_scores(type, type_short, frac, ylim=None):
    base_dir = f"/home/mwarr/Data/One_Genome/experiment2_10-90/chain_and_align/{type_short}_freq_{frac}"
    op = "score"
    title = f"Distribution of Top Scores (with Alignment), Alignment {frac * 100}%, {type}"
    xtitle = "Score Category and Region Type"
    ytitle = "Score"
    filename = f"box_align_chain_scores_{frac}_{type_short}.png"
    top_5_all(base_dir, op, title, xtitle, ytitle, filename, ylim=ylim)

'''
Generates box and whisker plots for the number of top scores
'''
def driver_counts():
    base_dir = "/home/mwarr/Data/One_Genome/other_features/local/counts"
    op = "count"
    title = "Distribution of the Number of Top Scores, Local"
    xtitle = "Score Category and Region Type"
    ytitle = "Number of Top Scores"
    filename = "box_counts_loc.png"
    top_5_all(base_dir, op, title, xtitle, ytitle, filename)

'''
Generates box and whisker plots for the number of anchors for a 5th highest chain, given a fixed chain length
'''
def driver_anchor_num():
    base_dir = "/home/mwarr/Data/One_Genome/other_features/local/anchor_num"
    op = "anchor_num"
    title = "Distribution of Number of Anchors for a 5th-Highest Chain Length"
    xtitle = "5th-Highest Chain Length"
    ytitle = "Number of Anchors"
    filename = "box_anchor_num_loc.png"
    fixed_chain_all(base_dir, op, title, xtitle, ytitle, filename)

'''
Generates box and whisker plots for the reference sequence length for a 5th highest chain, given a fixed chain length
'''
def driver_ref_len():
    base_dir = "/home/mwarr/Data/One_Genome/other_features/local/ref_len"
    op = "ref_len"
    title = "Distribution of the Reference Region Length of the 5th-Highest Chain Score"
    xtitle = "5th-Highest Chain Score"
    ytitle = "Reference Region Length"
    filename = "box_ref_len_loc.png"
    fixed_chain_all(base_dir, op, title, xtitle, ytitle, filename)

