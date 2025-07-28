from matplotlib import pyplot as plt #type: ignore


'''
This file contains miscellaneous helper functions for frequent tasks.
'''

'''
Populates <dict> with pairs of regions and a score for that pair. Returns
the min and max chain scores.
<filepath> should be the score file containing pairs of regions and a score for the pair.
'''
def create_score_dict(dict, filepath):
    max_chain = None
    min_chain = None
    #Put the chaining file into a dictionary: {(reg1, reg2): chain_score, ...}
    with open(filepath) as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            chain_score = float(line_arr[2])
            reg1 = min(line_arr[0], line_arr[1])
            reg2 = max(line_arr[0], line_arr[1])
            dict[(reg1, reg2)] = chain_score
            if max_chain == None or chain_score > max_chain:
                max_chain = chain_score
            if min_chain == None or chain_score < min_chain:
                min_chain = chain_score
    return (min_chain, max_chain)

'''
Takes in a file with a list of regions and outputs a fasta file of those regions.
Regions should be formatted chr#_<start>to<end> or chr#_<start>to<end>_rand
'''
def make_region_list_fasta(input_file, genome_dict, output_file):
    with open(output_file, "w") as output:
        with open(input_file) as input:
            for line in input:
                if line != "\n":
                    chr = line[ : line.index("_")]
                    start = int(line[line.index("_") + 1 : line.index("to")])
                    if "_rand" in line:
                        stop = int(line[line.index("to") + 2 : line.index("_rand")])
                    stop = int(line[line.index("to") + 2 : ])
                    seq = genome_dict[chr][start - 1 : stop]
                output.write(f">{line.strip()}\n{seq.strip()}\n")

'''
Gets the average chaining score from a file with the frequencies of chaining scores. 
In this file, each line should have the value, then a tab, then the frequency of the value.
If exclude_zero=True, zeros will not be included in the average.
'''
def get_average(file_path, exclude_zero=False):
    sum = 0
    num_regions = 0
    with open(file_path) as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                break
            if exclude_zero and float(line_arr[0]) == 0:
                continue
            sum += float(line_arr[0]) * int(line_arr[1])
            num_regions += int(line_arr[1])
    return sum / num_regions

'''
Takes in a file with a list of ACRs and creates a set
'''
def create_ref_set(input_path):
    ref_set = set()
    with open(input_path, "r") as file:
        for line in file:
            line = line.strip()
            if line != "":
                ref_set.add(line)
    return ref_set

'''
<genome_file> should be a fasta file of the genome.

Creates a dictionary where the keys are the chromosome identifiers
for the genome and the values are the sequences for the given chromosome
'''
def create_genome_dict(genome_file):
    #create genome dict {chr: seq, chr: seq, ...}
    genome_dict = {}
    with open(genome_file) as genome:
        chromosome = genome.readline().strip()
        while chromosome:
            chromosome = chromosome[1:]
            seq = ""
            line = genome.readline().strip()
            while line != "" and line[0] != ">":
                seq += line
                line = genome.readline().strip()
            genome_dict[chromosome] = seq
            chromosome = line
    return genome_dict

'''
Returns the length of an ACR in the format Chr#_{start_pos}to{end_pos} or Chr#_{start_pos}to{end_pos}_rand
'''
def get_length(item):
    start = int(item[item.index("_") + 1 : item.index("to")])
    if "rand" in item:
        end = int(item[item.index("to") + 2 : item.index("_rand")])
    else:
        end = int(item[item.index("to") + 2 : ])
    length = abs(end - start) + 1
    return length

'''
Returns a dictionary where the keys are the unique numbers in <lst> and the values are the
frequencies of each number in <lst>, i.e. {num: freq, num2: freq, ...}
'''
def get_frequencies(lst):
    count_dict = {}
    for item in lst:
        count_dict[item] = count_dict.get(item, 0) + 1
    return count_dict

'''
Takes in a dictionary where the keys are numbers and the values are the frequencies of the numbers.
Outputs a frequency file to <filepath> where each line has the number, then a tab, then the frequency of the number.
'''
def output_freq_to_file(filepath, freq_dict):
    with open(filepath, "w") as file:
        for val, freq in freq_dict.items():
            file.write(f"{val}\t{freq}\n")

'''
Reads frequencies from files and creates two side-by-side plots.
<file1> should contain ACR chaining with ACR frequency data.
<file2> should contain ACR chaining with random frequency data.
Plots will only include values between <min_index> and <max_index>.
<min_index> defaults to 0.
If equal_freq = True, this function will check that the total number of values in
the ACR file matches the total number of values in the rand file.
'''
def create_histograms(file1, file2, title, x_label, max_index, equal_freq=False, min_index=0):
    data1 = [0 for i in range((max_index - min_index) + 1)]
    acr_outliers = 0
    rand_outliers = 0

    with open(file1, "r") as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            score = round(float(line_arr[0]))
            freq = int(line_arr[1])
            try:
                data1[score-min_index] += freq
            except IndexError:
                acr_outliers += freq

    data2 = [0 for i in range((max_index - min_index) + 1)]
    with open(file2, "r") as file:
        for line in file:
            line_arr = line.split("\t")
            if len(line_arr) < 2:
                continue
            score = round(float(line_arr[0]))
            freq = int(line_arr[1])
            try:
                data2[score-min_index] += freq
            except IndexError:
                rand_outliers += freq

    #loop through data arrays backwards to find max index which is not 0
    for i in range(len(data1) - 1, -1, -1):
        if data1[i] != 0 or data2[i] != 0:
            max_index = i
            data1 = data1[: max_index + 1]
            data2 = data2[: max_index + 1]
            break

    if equal_freq:
        try:
            assert(sum(data1) + acr_outliers == sum(data2) + rand_outliers)
        except:
            print(sum(data1) + acr_outliers)
            print(sum(data2) + rand_outliers)
            raise

    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True, sharex=True)

    axes[0].bar(range(min_index, max_index+1), data1, color='blue')
    axes[0].set_title("ACRs chained with ACRs")
    axes[0].set_ylabel("Frequency")
    axes[0].set_xlabel(x_label)

    axes[1].bar(range(min_index, max_index+1), data2, color='green')
    axes[1].set_title("ACRs chained with random regions")
    axes[1].set_ylabel("Frequency")
    axes[1].set_xlabel(x_label)

    #plt.axis('equal')
    plt.suptitle(title)
    plt.tight_layout()
    file_name = title.replace(" ", "_")
    print("Saving figure")
    plt.savefig(f"/home/mwarr/{file_name}.png")
    print(f"ACR outliers {acr_outliers}")
    print(f"Random outliers {rand_outliers}")
