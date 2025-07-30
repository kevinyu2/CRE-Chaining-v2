# CRE-Chaining-v2

Suite of tools for chaining with Cis Regulatory Elements (CREs) to match Accessible Chromatin Regions (ACRs)

## Chaining

Global chaining: finds the length of the longest non-overlapping chain

Local chaining: finds the score for the best chaining in a local area (penalizes non-chained elements)

#### ```chaining.py```

Takes in anchors, pairs (a, c) linking a in the first string with c in the second

Contains driver functions for chaining (considering both possible orientations of the sequence): 
- ```chain_driver(anchors, is_weighted)``` -> global chaining, takes in list of [(a, c)] or [(a, c, weight)]
- ```chain_driver_np(anchors_np, is_weighted)``` -> takes in a np array with shape (n, 2) or (n, 3) for weighted
- ```chain_local_driver(anchors, match, mismatch, gap, is_weighted)``` -> if is_weighted, match score ignored
- ```chain_local_driver_np(anchors_np, match, mismatch, gap, is_weighted)``` -> if is_weighted, match score ignored

Contains individual chaining functions (only condiering given orientation)
- ```chain(anchors)```
- ```chain_np(anchors_np)```
- ```chain_local(anchors, match, mismatch, gap)```
- ```chain_local_np(anchors_np, match, mismatch, gap)```
- ```chain_weighted(anchors)```
- ```chain_weighted_np(anchors_np)```
- ```chain_local_weighted(anchors, match, mismatch, gap)```
- ```chain_local_weighted_np(anchors_np, match, mismatch, gap)```

Detailed info in the file header

### ```anchor_chain_driver.py```

Does pairwise chaining for a set of ACRs. Can be given direct fimo results (use the parent directory to all the fimo_out folders) by setting motif_mode to ```only_acr```. If using a list of motifs in the format ```ACR1[\n]MOTIF1[\n]MOTIF2[\n]...ACR2...```, set motif_mode to ```full```.

This can do global or local chaining, with weighted scores. ```SCORE_CENTER``` nudges all weighted scores to average at it. It is useful for local chaining because the gap and mismatch scores are typically fixed, so this allows different types of scoring systems to all average at the same match score.

### ```chain_two_sets.py```

Same as ```anchor_chain_driver```, but does chaining between two sets of ACRs. Faster if that is what is desired. The output always lists the ```ACR_MOTIF_SET_A``` acrs first.

## Full-Genome BLAST-Search (/BLAST_chain_search)

Finds areas with high local chaining scores between a sequence and a set of ACRs

### Step 1 (Optional): Cluster motifs

Use this if the motif set is large, or contains redundancies

- Run TOMTOM (MEME Suite) between a MEME motif file and itself

Ex:

```nohup tomtom -dist kullback -motif-pseudo 0.1 -text -min-overlap 1  /home/jm/meme/motif_databases/ARABD/ArabidopsisDAPv1.meme  /home/jm/meme/motif_databases/ARABD/ArabidopsisDAPv1.meme > /home/kyu/tomtom_JASPAR.txt 2> /home/kyu/tomtom_JASPAR.log &```

#### ```motif_clustering.py```

Outputs a cluster file detailing which motif belongs in which cluster (for post-clustered)
Outputs a cluster meme file which can be used directly in step 2 (for pre-clustered)

- Prefer pre-clustered if there are a large number of motifs (i.e. 500+)

### Step 2: Run FIMO (MEME Suite) on the whole sequence

If you used Step 1, you may put the resulting MEME file here (pre-clustered)
You can also run MEME with the original and provide a CLUSTER_FILE later when asked (post-clustered)

Ex: 

```nohup fimo --oc fimo_out_ArabidopsisDAPv1  /home/jm/meme/motif_databases/ARABD/ArabidopsisDAPv1.meme /home/projects/msu_nsf_pangenomics/pgrp/data/arabidopsis/atacseq/tair10.fa &```

### Step 3: Get motif sequences

These are the inputs for the BLAST search

#### ```get_motif_sequence.py```

Gets the motif sequence for each chromosome of the sequence

Outputs a folder with files for each chromosome

#### ```get_acr_motif_from_full.py```

Outputs the motif sequences for ACRs. 
File format: 
ACR
MOTIF
MOTIF
...
ACR
MOTIF
...

If using a separate set of ACRs, try to replicate this file

### Step 4 (Optional): Weighted scoring (/BLAST_chain_search/Weighted_Scoring)

This is heavily recommended

#### ```motif_relevance.py```

Gets info about motifs in ACR files and in the full background

#### ```create_motif_scoring.py```

Creates a file from the results of motif_relevance.py that can be inputted into Step 5 to get a varied scoring

### Step 5: Blast Search!

Finds locally high chain scores through finding matching k-mers, then extending those hits

#### ```v2_parallel_kmer_CRE_BLAST.py```

Parallelized to improve runtime. 
Some hits are duplicates! If two different seeds get the same extension, they end up with the same score
Output gives info on location of hits and score

### Step 6: Get Alignment Scores

If running ACRs from the same sequence, high aligning hits tend to be duplications. This allows us to remove those

#### ```chain_but_not_align.py```

Gets the alignment scores for all hits above a certain threshold
Also removes duplicate hits

#### ```create_eval_input.py```

Simply filters for high scoring and low alignment

### Step 7: Evaluation

#### ```evaluate_blast.py```

takes in a test set of ACRs, sees how many bases in the hits are from that test set

#### ```difference_blast.py```

Same as evaluate_blast.py, but takes in two files and sees how the differences score

## ACR Clustering (/ACR_Clustering)

Calculate chain scores between ACRs and cluter them
Look at expression correlations within clusters

### Step 1: Get motif sequences & weighted scores

Follow Full-Genome BLAST search steps 1-4. No need to run ```get_motif_sequence.py```.

### Step 2: Get pairwise chain scores

Outputs a file with chain scores (global or local, potentially weighted) for each ACR with each other

#### ```../anchor_chain_driver.py```

Parallelized chaining. Use custom score if you ran Step 4 in Full-Genome BLAST

### Step 3: Get Clusters

Use hierarchical clustering to get clusters of ACRs based on chaining scores

### ```get_distances.py```

There are two different versions, for local and global. For local, BETA is a multiplier. For global, ALPHA determines how much the longest of the two sequence lengths is used, at 1 it is always the longest

### ```hierarchical.py```

Produces a cluster file detailing which thing is in which cluster. If not interested in expression, you can stop here.

### Step 4: Convert Clusters to Gene Sets

Turn from ACRs to the genes themselves

##### ```\ACR_to_Gene\acr_to_bed.py```

Get a bed file from acrs

#### ```\ACR_to_Gene\gff_to_bed.py```

Turn a GFF file of the genes in the genome to a bed file

Now with both bed files you can run ```bedtools intersect -wa -wb -a acr.bed -b gff.bed > acr_to_gene_promoter.bed```

To get an acr_to_gene file

#### ```get_gene_sets_from_clusters.py```

Converts each cluster to sets of genes


### Step 5 (Optional): Evaluate correlation of Expression

Create graphs comparing random expression correlation to the created gene sets

#### ```evaluate_expression.ipynb```


### Step 6 (For Feature Vector): Find Representative Sequence

For each cluster, get a new sequence of CREs that defines it
This uses central star MSA, and takes as the consensus a sequence of CRE's which 
have at least 50% representation in the cluster

#### ```Cluster_MSA.py```

## ACRs vs Random Regions (/ACR_rand_compare)

Compares the chaining scores between random regions and ACRs with the scores between ACRs and ACRs.

### Step 1: Create Region Sets

#### ```generate_rand_sets.py```

Determine which ACRs should be in the TEST group (also referred to as non-ref ACRs in the documentaion)
and which should be in the REFERENCE group. We randomly assigned 10% of the 31k *Arabidopsis thaliana* ACRs to
TEST and the remaining 90% to REFERENCE.

#### ```random_regions.py```

The ```random_region_match_sizes``` function can be used to generate the random regions with the same number and region sizes as TEST. We call this group TEST-RAND.

### Step 2: Run FIMO (MEME Suite) on ACRs and Random Regions

The easiest way to do this is to create a fasta file with all of the ACRs + random regions and then
run a command similar to the one below:

```nohup fimo --oc fimo_out_ArabidopsisDAPv1  /home/jm/meme/motif_databases/ARABD/ArabidopsisDAPv1.meme /home/mwarr/Data/acr_and_rand.fa &```

### Step 3: Pairwise Chain ACRs and Random Regions

#### ```../anchor_chain_driver```

Set ```motif_mode``` to 'only_acr', choose the desired chain mode, and provide the folder containing all the
FIMO folders to ```MOTIFS```. 

### Step 4 (Optional): Find Alignment Scores

#### ```get_alignment_scores.py```

Use ```split_fasta``` to split up the fasta files containing TEST and TEST-RAND. The number of files
depends on the speed you need and the number of cores you have available (we split each into 50).
Then, run ```gen_align_scores.sh``` to launch the program for each file. Once the alignment is finished, concatenate all the files into one file with this command:

```for file in align_temp; do cat $file >> align_scores.tsv; done;```

### Step 5: Generate Desired Frequency Files
Depending on what you would like to analyze, you may run any of the following:

#### ```frequencies.py```

This file contains the following driver functions:
- ```driver_all_summary```: Outputs score frequencies and histograms
- ```driver_filter_summary```: Outputs score frequencies of regions within specific size intervals
- ```driver_top```: Outputs frequency files for top 5 score ranks
- ```driver_top_normal```: Outputs frequency files for top 5 score ranks; normalizes scores.

#### ```freq_with_align.py```

```exclude_high_align_driver``` outputs frequency files for chaining scores, excluding pairs with high alignment.

#### ```chain_features.py```

```driver_frequencies_all``` Outputs frequency files for the number of top scores, the number of anchors, the reference sequence length, and a combined score.

### Step 6: Graph and Analyze

#### ```box_whisker.py```
This is the recommended analysis method. Has functions to create various box and whisker plots, depending on what your frequency data is.

#### ```bar_graph.py```
Various functions for creating bar graphs of averages across various statistics. Not as informative.
