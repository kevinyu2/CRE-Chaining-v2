from matplotlib import pyplot as plt #type: ignore
import numpy as np #type: ignore
from useful_functions import *

'''
Generates bar graphs for averages across various statistics of chains.
'''

'''
Creates a bar graph showing the average ACR v ACR score and ACR v random score,
as well as the average max ACR v ACR score and ACR v random score. Shows both
local and global chaining score averages.

<input_dir_base> should be a directory containing subdirectories "global" and "local" with 
frequency files named "{ACR|rand}_vs_ACR_{all|max}_freq.tsv" for every combination of region
type and score type (all or max).
'''
def create_average_bar(input_dir_base):
    options = ["all", "max"]
    chain_type = ["local", "global"]
    labels = ["All, Local", "Max, Local", "All, Global", "Max, Global"]
    ACR_data = []
    rand_data = []
    for type in chain_type:
        for opt in options:
            ACR_data.append(get_average(f"{input_dir_base}/{type}/ACR_vs_ACR_{opt}_freq.tsv"))
            rand_data.append(get_average(f"{input_dir_base}/{type}/rand_vs_ACR_{opt}_freq.tsv"))
    print(ACR_data)
    print(rand_data)
    width = .4
    x = np.arange(len(labels))
    plt.bar(x - width / 2, ACR_data, label="ACRs", width=width)
    plt.bar(x + width / 2, rand_data, label="Random", width=width)
    plt.xticks(x, labels)
    plt.title("Comparison of ACR Chaining Score Averages")
    plt.xlabel("Data Averaged and Chain Type")
    plt.ylabel("Average Score")
    plt.legend()
    plt.savefig("/home/mwarr/averages_all-sizes.png")

'''
Creates a graph similar to the function above (create_average_graph). However, it will
plot local and global on different subplots so that the axes and scale can be different 
for each.
'''
def create_average_bar_subgraphs(input_dir_base):
    options = ["all", "max"]
    chain_type = ["local", "global"]
    labels = ["All, Local", "Max, Local", "All, Global", "Max, Global"]
    ACR_data = []
    rand_data = []
    for type in chain_type:
        for opt in options:
            ACR_data.append(get_average(f"{input_dir_base}/{type}/ACR_vs_ACR_{opt}_freq.tsv"))
            rand_data.append(get_average(f"{input_dir_base}/{type}/rand_vs_ACR_{opt}_freq.tsv"))
    width = .4
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    x = np.arange(2)

    
    axes[0].bar(x - width / 2, ACR_data[0:2], label="ACR", width=width)
    axes[0].bar(x + width / 2, rand_data[0:2], label="Rand", width=width)
    axes[0].set_title("Local Alignment Scores")
    plt.sca(axes[0])
    plt.xlabel("Data Averaged and Chain Type")
    plt.ylabel("Average Score")
    plt.xticks(x, labels[0:2])
    plt.legend()

    axes[1].bar(x - width / 2, ACR_data[2:4], label="ACR", width=width)
    axes[1].bar(x + width / 2, rand_data[2:4], label="Rand", width=width)
    axes[1].set_title("Global Alignment Scores")
    plt.sca(axes[1])
    plt.xlabel("Data Averaged and Chain Type")
    plt.xticks(x, labels[2:4])
    plt.ylim(11000, 12500)
    plt.legend()
    
    fig.suptitle("Comparison of Alignment Score Averages")
    plt.savefig("/home/mwarr/averages_alignment.png")

'''
Creates a bar graph of the average top scores. There will be bars for the max, 2nd-highest, ..., <top_num>-highest
scores.

<input_dir_base> should be a directory containing files named "{ACR|rand}_vs_ACR_{num}-highest_freq.tsv".

If summary=True, the bar graph will also include a bar averaging the top scores. There should be additional frequency files
named "{ACR|rand}_vs_ACR_top-{top_num}-avg_freq.tsv".
'''
def average_graph_bar_top(input_dir_base, top_num, title, out_file, summary=False):
    acr_data = []
    rand_data = []
    labels = ["Max", "2nd-Highest", "3rd-Highest"]
    for i in range(1, top_num + 1):
        if i > 3:
            labels.append(f"{i}th-Highest")
        acr_data.append(get_average(f"{input_dir_base}/ACR_vs_ACR_{i}-highest_freq.tsv"))
        rand_data.append(get_average(f"{input_dir_base}/rand_vs_ACR_{i}-highest_freq.tsv"))
    if summary:
        labels.append(f"Top-{top_num}-Highest")
        acr_data.append(get_average(f"{input_dir_base}/ACR_vs_ACR_top-{top_num}-avg_freq.tsv"))
        rand_data.append(get_average(f"{input_dir_base}/rand_vs_ACR_top-{top_num}-avg_freq.tsv"))
    print(acr_data)
    print(rand_data)
    width = .4
    x = np.arange(len(labels))
    plt.figure(figsize=(15, 6))
    plt.bar(x - width / 2, acr_data, label="ACRs", width=width)
    plt.bar(x + width / 2, rand_data, label="Random", width=width)
    plt.xticks(x, labels)
    plt.title(title)
    plt.xlabel("Data Averaged")
    plt.ylabel("Average Score")
    plt.legend()
    plt.savefig(out_file)

'''
Reads in the frequency files from <input_dir_base> which contain the frequencies of some
feature with a fixed chain length and generates a bar graph of the averages.

<op_name> will be part of the output filename.
'''
def average_graph_constant_chain(input_dir_base, op_name, ylabel, title):
    acr_data = []
    rand_data = []
    labels = [i for i in range(2, 20)]
    for i in range(2, 20):
        acr_data.append(get_average(f"{input_dir_base}/ACR_vs_ACR_{op_name}_{i}_freq.tsv", True))
        rand_data.append(get_average(f"{input_dir_base}/rand_vs_ACR_{op_name}_{i}_freq.tsv", True))
    width = .4
    x = np.arange(len(labels))
    plt.figure(figsize=(15, 6))
    plt.bar(x - width / 2, acr_data, label="ACRs", width=width)
    plt.bar(x + width / 2, rand_data, label="Random", width=width)
    plt.xticks(x, labels)
    plt.ylim(870, 890)
    plt.title(title)
    plt.xlabel("5th-Highest Chain Length")
    plt.ylabel(ylabel)
    plt.legend()
    plt.savefig(f"/home/mwarr/chainging_exp2_{op_name}_glob_original.png")