from collections import defaultdict



CLUSTER_FILE = '/home/mwarr/Data/Clustering/hierarchical_DAPv1_clustered_local.txt'

GENE_BED = '/home/kyu/CRE-chaining/Single-Genome/ACR_to_Gene/acr_to_gene_promoter.bed'

OUT_FILE = '/home/mwarr/Data/Clustering/gene_set_hierarchical_DAPv1_clustered_local.txt'

# Produce matching random file for analysis
RANDOM_OUT_FILE = '/home/mwarr/Data/Clustering/gene_set_random_hierarchical_DAPv1_clustered_local.txt'

ACRS_LIST = '/home/mwarr/Data/all_acr.txt'

acr_gene_dict = defaultdict(set)

acr_gene_biggest_overlap_dict =  {}

# Parse gene match file
with open(GENE_BED, 'r') as gb :
    for line in gb :
        line_arr = line.rstrip().split('\t')
        acr_gene_dict[line_arr[3]].add(line_arr[7])

        end = min(int(line_arr[2]), int(line_arr[6]))
        start = max(int(line_arr[1]), int(line_arr[5]))

        if line_arr[3] in acr_gene_biggest_overlap_dict :
            if end - start >  acr_gene_biggest_overlap_dict[line_arr[3]][1] :
                acr_gene_biggest_overlap_dict[line_arr[3]] = (line_arr[7], end - start)
        else :

            acr_gene_biggest_overlap_dict[line_arr[3]] = (line_arr[7], end - start)

print(acr_gene_biggest_overlap_dict)

print(acr_gene_biggest_overlap_dict['Chr5_26833111to26834211'][0])


cluster_gene_dict = defaultdict(set)

with open(CLUSTER_FILE, 'r') as cf :
    current_cluster = ""
    for line in cf :
        if 'Cluster no:' in line :
            current_cluster = line.rstrip().split('Cluster no: ')[1]


        else :
            line_arr = line.rstrip().split('\t')
            for acr in line_arr :
                # print(acr)
                # print(acr_gene_dict[acr])
                # print(current_cluster)

                cluster_gene_dict[current_cluster] = cluster_gene_dict[current_cluster].union(acr_gene_dict[acr])

with open(OUT_FILE, 'w') as out :
    for cluster, gene_set in cluster_gene_dict.items() :
        out.write(cluster)
        out.write('\n')
        if len(gene_set) > 1:
            out.write('\t'.join(str(item) for item in gene_set))    
        elif len(gene_set) > 0 :
            out.write(next(iter(gene_set)))
        out.write('\n')


cluster_gene_dict_one = defaultdict(set)

with open(CLUSTER_FILE, 'r') as cf :
    current_cluster = ""
    for line in cf :
        if 'Cluster no:' in line :
            current_cluster = line.rstrip().split('Cluster no: ')[1]


        else :
            line_arr = line.rstrip().split('\t')
            for acr in line_arr :
                # print(acr)
                # print(acr_gene_dict[acr])
                # print(current_cluster)
                if len(acr_gene_dict[acr]) >= 1 :
                    cluster_gene_dict_one[current_cluster] = cluster_gene_dict_one[current_cluster].union({acr_gene_biggest_overlap_dict[acr][0]})

with open(f"{OUT_FILE}_one", 'w') as out :
    for cluster, gene_set in cluster_gene_dict_one.items() :
        out.write(cluster)
        out.write('\n')
        if len(gene_set) > 1:
            out.write('\t'.join(str(item) for item in gene_set))    
        elif len(gene_set) > 0 :
            out.write(next(iter(gene_set)))
        out.write('\n')

################################################################################
# RANDOM

import random

acrs = []
with open(ACRS_LIST, 'r') as acr_list_file :
    for line in acr_list_file :
        acrs.append(line.rstrip())

random_cluster_gene_dict = defaultdict(set)

with open(CLUSTER_FILE, 'r') as cf :
    current_cluster = ""
    for line in cf :
        if 'Cluster no:' in line :
            current_cluster = line.rstrip().split('Cluster no: ')[1]


        else :
            line_arr = line.rstrip().split('\t')

            random_acr_list = random.sample(acrs, len(line_arr))
            for acr in random_acr_list :
                # print(acr)
                # print(acr_gene_dict[acr])
                # print(current_cluster)
                if len(acr_gene_dict[acr]) > 0 :
                    random_cluster_gene_dict[current_cluster] = random_cluster_gene_dict[current_cluster].union(acr_gene_dict[acr])


with open(f"{RANDOM_OUT_FILE}", 'w') as out :
    for cluster, gene_set in random_cluster_gene_dict.items() :
        out.write(cluster)
        out.write('\n')
        if len(gene_set) > 1:
            out.write('\t'.join(str(item) for item in gene_set))    
        elif len(gene_set) > 0 :
            out.write(next(iter(gene_set)))
        out.write('\n')







cluster_gene_dict_random_one = defaultdict(set)

with open(CLUSTER_FILE, 'r') as cf :
    current_cluster = ""
    for line in cf :
        if 'Cluster no:' in line :
            current_cluster = line.rstrip().split('Cluster no: ')[1]


        else :
            line_arr = line.rstrip().split('\t')
            random_acr_list = random.sample(acrs, len(line_arr))
            for acr in random_acr_list :
                # print(acr)
                # print(acr_gene_dict[acr])
                # print(current_cluster)
                if len(acr_gene_dict[acr]) >= 1 :
                    cluster_gene_dict_random_one[current_cluster] = cluster_gene_dict_random_one[current_cluster].union({acr_gene_biggest_overlap_dict[acr][0]})

with open(f"{RANDOM_OUT_FILE}_one", 'w') as out :
    for cluster, gene_set in cluster_gene_dict_random_one.items() :
        out.write(cluster)
        out.write('\n')
        if len(gene_set) > 1:
            out.write('\t'.join(str(item) for item in gene_set))    
        elif len(gene_set) > 0 :
            out.write(next(iter(gene_set)))
        out.write('\n')