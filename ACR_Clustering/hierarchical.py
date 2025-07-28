import numpy as np
from sklearn.cluster import AgglomerativeClustering
from collections import defaultdict
from sklearn.metrics import silhouette_score
# import umap
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from tqdm import tqdm
import itertools

#####################################################################

# Hierarchical clustering threshold
distance_threshold = 0.2
distance_file = '/home/mwarr/Data/distances_DAPv1_clustered_local.tsv'
OUTFILE = '/home/mwarr/Data/Clustering/hierarchical_DAPv1_clustered_local.txt'


# Only output clusters of a certain size
CLUSTER_MIN_SIZE = 4

#######################################################################

with open(OUTFILE, 'w') as out:
    out.write(f"Using Distance Threshold {distance_threshold}\n")

# Maps each label to an index
index_dict = {}
curr_idx = 0

print("Reading Distance File (First Pass)")

with open(distance_file, 'r') as distances :
    for line in distances:
        line_arr = line.split('\t')
        if line_arr[0] not in index_dict :
            index_dict[line_arr[0]] = curr_idx
            curr_idx += 1
        if line_arr[1] not in index_dict :
            index_dict[line_arr[1]] = curr_idx
            curr_idx += 1
# Maps indexes to labels
label_dict = {v: k for k, v in index_dict.items()}

labels = index_dict.keys()

distance_matrix = np.full((curr_idx, curr_idx), 1.0)
np.fill_diagonal(distance_matrix, 0)

print("Reading Distance File (Second Pass)")

with open(distance_file, 'r') as distances :
    for line in distances:
        line_arr = line.split('\t')
        if float(line_arr[2]) < 0:
            line_arr[2] = 0
        distance_matrix[index_dict[line_arr[0]]][index_dict[line_arr[1]]] = line_arr[2]
        distance_matrix[index_dict[line_arr[1]]][index_dict[line_arr[0]]] = line_arr[2]

# Cluster
clustering = AgglomerativeClustering(
    metric='precomputed', 
    linkage='average',
    distance_threshold=distance_threshold,
    n_clusters=None
)

print("Clustering...")
labels = clustering.fit_predict(distance_matrix)


cluster_dict = defaultdict(set)
for i, lab in enumerate(labels) :
    cluster_dict[lab].add(i)

print("Outputting...")
with open(OUTFILE, 'w') as out :

    for cluster_no, cluster in tqdm(cluster_dict.items()) :
        if len(cluster) > CLUSTER_MIN_SIZE :

            if len(cluster) > 1:
                # Get the representative
                representative = -1
                min_distances = float('inf')
                for c in cluster :
                    curr_sum = 0
                    # calculate average distances
                    for d in cluster :
                        if c != d :
                            curr_sum += distance_matrix[c][d] 

                    if curr_sum < min_distances :
                        representative = c
                        min_distances = curr_sum
            # if the cluster size is just 1
            elif len(cluster) == 1 :
                representative = cluster[0]
            else :
                representative = -1
                label_dict[-1] = "None"

            out.write(f"Cluster no: {cluster_no}, Representative: {label_dict[representative]}\n")
            for i, c in enumerate(cluster) :
                if i != 0 :
                    out.write("\t")
                out.write(f"{label_dict[c]}")
            out.write("\n")

