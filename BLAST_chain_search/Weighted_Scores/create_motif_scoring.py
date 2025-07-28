# Given information like motif relevance, create a scoring metric
# Scoring system is log(pacr/pbg) - lambda(log(pbg))
import pandas as pd

# Generated from motif_relevance.py
motif_relevance_file = '/home/mwarr/Data/motif_relevance/PBM_relevance_seta.tsv'

score_out_file = '/home/mwarr/Data/motif_relevance/scores_PBM_lambda08.tsv'

relevance_df = pd.read_csv(motif_relevance_file, sep = '\t')

LAMBDA = 0.8


# Proportions for each base in ACRs and in the background
p_hit = {}
p_bg = {}
for i, row in relevance_df.iterrows() :
    p_hit[row['Motif']] = row['ACR_percent']
    p_bg[row['Motif']] = row['Full_percent']



import math

def compute_match_scores(p_hit, p_bg, lam=0.0):
    """
    Compute match scores for A, C, G, T based on log-likelihood ratio 
    and optional rarity boost.

    Parameters:
        p_hit (dict): Foreground frequencies, e.g., {'A': 0.3, 'C': 0.3, 'G': 0.1, 'T': 0.3}
        p_bg (dict): Background frequencies, e.g., {'A': 0.4, 'C': 0.1, 'G': 0.1, 'T': 0.4}
        lam (float): Weight for rarity adjustment (default: 0.0)

    Returns:
        dict: Match scores for each base
    """
    bases = p_bg.keys()
    scores = {}
    for base in bases:
        if p_hit[base] == 0 or p_bg[base] == 0:
            raise ValueError(f"Probabilities must be > 0 for base {base}")
        enrichment = math.log2(p_hit[base] / p_bg[base])
        rarity = -math.log2(p_bg[base])
        scores[base] = enrichment + lam * rarity
    return scores


scores = compute_match_scores(p_hit, p_bg, LAMBDA)

with open(score_out_file, 'w') as out :
    out.write('Motif\tScore\tACR_percent\tFull_percent\n')
    for motif in p_bg.keys() :
        out.write(f"{motif}\t{scores[motif]}\t{p_hit[motif]}\t{p_bg[motif]}\n")
