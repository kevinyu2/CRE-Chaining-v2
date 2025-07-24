import bisect

import numpy as np

'''
Performs co-linear chainging on anchors.

Driver Functions:

chain_driver(anchors, is_weighted): takes in a list of tuples, and if it is weighted. Weighted tuples are [(a, c, weight)]
    Calculates the global optimal chain length considering BOTH directions (i.e. forwards or backwards)   
chain_local_driver(anchors, match, mismatch, gap, is_weighted): takes in a list of tuples, and if it is weighted. Weighted tuples are [(a, c, weight)]
    If weighted, match scores are based on weights. Mismatch and Gap should be negative. Suggested: 5, -1, -1
    Calculates the global optimal chain length considering BOTH directions (i.e. forwards or backwards) 
chain_driver_np(anchors_np, is_weighted): Same as chain_driver, but takes in a numpy array (n, 2) instead
chain_local_driver_np(anchors_np, is_weighted): Same as chain_local_driver, but takes in a numpy array (n, 2) instead


    
Individual Functions:


These are for lists
chain(anchors)
chain_weighted(anchors)
chain_local(anchors, match, mismatch, gap)
chain_local_weighted(anchors, match, mismatch, gap)
    Match is ignored, kept for consistency. The match scores are the weights instead

These are for np arrays
chain_np(anchors_np)
chain_weighted_np(anchors_np)
chain_local_np(anchors_np, match, mismatch, gap)
chain_local_weighted_np(anchors_np, match, mismatch, gap)
    Match is ignored, kept for consistency. The match scores are the weights instead

'''


############################################################################
# CALL THESE FUNCTIONS

def chain_driver(anchors, is_weighted):
    if is_weighted:
        negative_anchors = [(anchor[0], -1 * anchor[1], anchor[2]) for anchor in anchors]
        return max(chain_weighted(anchors), chain_weighted(negative_anchors))

    else :
        negative_anchors = [(anchor[0], -1 * anchor[1]) for anchor in anchors]
        return (max(chain(anchors), chain(negative_anchors)))


def chain_local_driver(anchors, match, mismatch, gap, is_weighted):
    if is_weighted:
        negative_anchors = [(anchor[0], -1 * anchor[1], anchor[2]) for anchor in anchors]
        return max(chain_local_weighted(anchors, match, mismatch, gap), chain_local_weighted(negative_anchors, match, mismatch, gap), key=lambda x: x[0])

    else :
        negative_anchors = [(anchor[0], -1 * anchor[1]) for anchor in anchors]
        return max(chain_local(anchors, match, mismatch, gap), chain_local(negative_anchors, match, mismatch, gap), key=lambda x: x[0])

def chain_driver_np(anchors_np, is_weighted):
    # Flip second coordinate
    negative_anchors_np = anchors_np.copy()
    negative_anchors_np[:, 1] *= -1
    
    if is_weighted :
        return max(
            chain_np_weighted(anchors_np),
            chain_np_weighted(negative_anchors_np)
        )
    else :
        return max(
            chain_np(anchors_np),
            chain_np(negative_anchors_np)
        )

def chain_local_driver_np(anchors_np, match, mismatch, gap, is_weighted):
    # if weighted, the match score is ignored
    # Flip second coordinate
    negative_anchors_np = anchors_np.copy()
    negative_anchors_np[:, 1] *= -1

    if is_weighted :
        return max(
            chain_local_np_weighted(anchors_np, mismatch, gap),
            chain_local_np_weighted(negative_anchors_np, mismatch, gap)
        )
    else :
        return max(
            chain_local_np(anchors_np, match, mismatch, gap),
            chain_local_np(negative_anchors_np, match, mismatch, gap)
        )




####################################################################################
# Helpers

# Do argsort but break ties by taking the later one in the array first
# This is done by just adding a small amount less than 1 to each in the array,
# Giving a greater value to the earlier ones (so argsort prioritizes later ones)
def argsort_reverse_ties(arr) :
    len_arr = len(arr)

    # Increment each so argsort does it in reverse tiebreaker
    for i in range(len_arr) :
        arr[i] += float((len_arr - i - 1) / len_arr)

    return np.argsort(arr)



def argsort_reverse_ties_np(arr):
    # arr is a 1D numpy array
    len_arr = len(arr)
    # Make a copy so we don't modify original
    arr_adj = arr.astype(float).copy()
    # Add small decreasing increments for reverse tiebreak
    increments = np.linspace((len_arr - 1) / len_arr, 0, len_arr)
    arr_adj += increments
    return np.argsort(arr_adj)


# From a leetcode submission
def lengthOfLIS(nums) :

    sub = []
    for x in nums:
        if len(sub) == 0 or sub[-1] < x:
            sub.append(x)
        else:
            idx = bisect.bisect_left(sub, x)  # Find the index of the first element >= x
            sub[idx] = x  # Replace that number with x

    return len(sub)


class FenwickTree:
    def __init__(self, size):
        self.tree = [0] * (size + 2)

    def update(self, i, value):
        while i < len(self.tree):
            self.tree[i] = max(self.tree[i], value)
            i += i & -i

    def query(self, i):
        res = 0
        while i > 0:
            res = max(res, self.tree[i])
            i -= i & -i
        return res

def weighted_LIS(A, W):
    # Coordinate compression
    sorted_unique = sorted(set(A))
    comp = {v: i+1 for i, v in enumerate(sorted_unique)}  # 1-based for BIT
    A_comp = [comp[x] for x in A]
    
    ft = FenwickTree(len(sorted_unique) + 2)
    max_total_weight = 0
    
    for i in range(len(A)):
        best_prev = ft.query(A_comp[i] - 1)  # max weight for values < A[i]
        curr_weight = best_prev + W[i]
        ft.update(A_comp[i], curr_weight)
        max_total_weight = max(max_total_weight, curr_weight)
    
    return max_total_weight

# BIT

class PrefixMaxBIT:
    def __init__(self, size):
        self.n = size
        self.tree = [float(0)] * (self.n + 1)  # 1-based index

    def update(self, i, val):
        while i <= self.n:
            self.tree[i] = max(self.tree[i], val)
            i += i & -i

    def query(self, i):
        res = float('-inf')
        while i > 0:
            res = max(res, self.tree[i])
            i -= i & -i
        return res




###################################################################################

#anchors should be a list of tuples. Each entry in the list is an MEM. For each MEM, 
#there is a tuple where the first entry is the index of the match in the first sequence
#and the second entry is the index of the match in the second sequence (a, c)
def chain(anchors):
    anchors_len = len(anchors)
    if anchors_len == 0 :
        return 0
    
    #sort anchors by second value (but tiebreak by the first value, descending)
    anchors.sort(key=lambda x: (x[1], -x[0]))
    
    #get list of first values (a)
    anchors_a = [element[0] for element in anchors]

    # Order in terms of a
    order = argsort_reverse_ties(anchors_a)

    return lengthOfLIS(order)


# anchors are [(position1, position2, weight)]
def chain_weighted(anchors):
    anchors_len = len(anchors)
    if anchors_len <= 1 :
        return anchors_len
    
    #sort anchors by second value (but tiebreak by the first value, descending)
    anchors.sort(key=lambda x: (x[1], -x[0]))
    
    #get list of first values (a)
    arr = np.array(anchors) 
    anchors_a = arr[:, 0] 

    # Order in terms of a
    order = argsort_reverse_ties(anchors_a)

    bit = PrefixMaxBIT(anchors_len)
    for i in order:
        maxPrev = bit.query(i + 1)
        bit.update(i + 1, maxPrev + anchors[i][2])
    
    return bit.query(anchors_len)



# anchors are numbered based on how many motifs apart (i.e not counting base pair location, but each motif is 1, 2, etc)
# match, mismatch, and gap scores like alignment
# Returns best local chain in the form of chain_score, chain_len (where the len is the best among all with the best chain score)
def chain_local(anchors, match, mismatch, gap) :
    anchors_len = len(anchors)
    if anchors_len == 0 :
        return 0,0
    
    #sort anchors by second value (but tiebreak by the first value, descending)
    anchors.sort(key=lambda x: (x[1], -x[0]))
    

    #get list of first values (a)
    anchors_a = [element[0] for element in anchors]

    # Order in terms of a
    order = argsort_reverse_ties(anchors_a)

    # basically a tuple, contains for each element (score, length)
    score_list = np.zeros((anchors_len, 2))

    for mem_idx in order :
        max_score = match
        max_length = 1
        for j in range(mem_idx) :
            # Non empty
            if score_list[j][0] != 0 :
                # score is score + gap (start - start - end + end) + mismatch (min(start - start, end - end) - 1)
                score_with_j = match + score_list[j][0] + gap * abs(anchors[mem_idx][0] - anchors[j][0] - anchors[mem_idx][1] + anchors[j][1]) + mismatch * (min(anchors[mem_idx][0] - anchors[j][0], anchors[mem_idx][1] - anchors[j][1]) - 1)
                
                # If we found a new best
                if score_with_j >= max_score :
                    max_score = score_with_j
                    max_length = max(max_length, score_list[j][1] + 1)
                
        score_list[mem_idx][0] = max_score
        score_list[mem_idx][1] = max_length
    
    # Get best score
    max_score = np.max(score_list[:, 0])
    rows_with_max_col0 = score_list[score_list[:, 0] == max_score]
    max_col = rows_with_max_col0[np.argmax(rows_with_max_col0[:, 1])]

    return max_col[0], int(max_col[1])


# anchors are numbered based on how many motifs apart (i.e not counting base pair location, but each motif is 1, 2, etc)
# mismatch and gap scores like alignment. No match score is a multiplier to the weight.
# Returns best local chain in the form of chain_score, chain_len (where the len is the best among all with the best chain score)
def chain_local_weighted(anchors, match, mismatch, gap) :
    anchors_len = len(anchors)
    if anchors_len == 0 :
        return 0
    
    #sort anchors by second value (but tiebreak by the first value, descending)
    anchors.sort(key=lambda x: (x[1], -x[0]))
    
    #get list of first values (a)
    anchors_a = [element[0] for element in anchors]

    # Order in terms of a
    order = argsort_reverse_ties(anchors_a)

    # basically a tuple, contains for each element (score, length)
    score_list = np.zeros((anchors_len, 2))

    for mem_idx in order :
        max_score = match * anchors[mem_idx][2]
        max_length = 1
        for j in range(mem_idx) :
            # Non empty
            if score_list[j][0] != 0 :
                # score is score + gap (start - start - end + end) + mismatch (min(start - start, end - end) - 1)
                score_with_j = match * anchors[mem_idx][2] + score_list[j][0] + gap * abs(anchors[mem_idx][0] - anchors[j][0] - anchors[mem_idx][1] + anchors[j][1]) + mismatch * (min(anchors[mem_idx][0] - anchors[j][0], anchors[mem_idx][1] - anchors[j][1]) - 1)
                
                # If we found a new best
                if score_with_j >= max_score :
                    max_score = score_with_j
                    max_length = max(max_length, score_list[j][1] + 1)
                
        score_list[mem_idx][0] = max_score
        score_list[mem_idx][1] = max_length

    
    # Get best score
    max_score = np.max(score_list[:, 0])
    rows_with_max_col0 = score_list[score_list[:, 0] == max_score]
    max_col = rows_with_max_col0[np.argmax(rows_with_max_col0[:, 1])]

    return max_col[0], int(max_col[1])

################################################################
# Versions to work with numpy arrays 


def chain_np(anchors_np):
    anchors_len = anchors_np.shape[0]
    if anchors_len == 0:
        return 0

    # Sort anchors by second value ascending, first value descending (tiebreak)
    # anchors_np[:,1] is second col, anchors_np[:,0] is first col
    order = np.lexsort((-anchors_np[:,0], anchors_np[:,1]))
    sorted_anchors = anchors_np[order]

    # Extract first values (a)
    anchors_a = sorted_anchors[:, 0]

    # Compute order with reverse tie-breaks
    order_for_lis = argsort_reverse_ties_np(anchors_a)

    # length of longest increasing subsequence
    return lengthOfLIS(order_for_lis)



def chain_local_np(anchors_np, match, mismatch, gap):
    anchors_len = anchors_np.shape[0]
    if anchors_len == 0:
        return 0, 0  # match score 0, length 0
    
    # Sort anchors by second value ascending, first value descending (tiebreak)
    order = np.lexsort((-anchors_np[:, 0], anchors_np[:, 1]))
    sorted_anchors = anchors_np[order]
    
    # Extract first values (a)
    anchors_a = sorted_anchors[:, 0]
    
    # Get order with reverse tie-break
    order = argsort_reverse_ties_np(anchors_a)

    # basically a tuple, contains for each element (score, length)
    score_list = np.zeros((anchors_len, 2))

    for mem_idx in order :
        max_score = match
        max_length = 1
        for j in range(mem_idx) :
            # Non empty
            if score_list[j][0] != 0 :
                # score is score + gap (start - start - end + end) + mismatch (min(start - start, end - end) - 1)
                score_with_j = match + score_list[j][0] + gap * abs(sorted_anchors[mem_idx][0] - sorted_anchors[j][0] - sorted_anchors[mem_idx][1] + sorted_anchors[j][1]) + mismatch * (min(sorted_anchors[mem_idx][0] - sorted_anchors[j][0], sorted_anchors[mem_idx][1] - sorted_anchors[j][1]) - 1)
                
                # If we found a new best
                if score_with_j >= max_score :
                    max_score = score_with_j
                    max_length = max(max_length, score_list[j][1] + 1)
                
        score_list[mem_idx][0] = max_score
        score_list[mem_idx][1] = max_length
    
    # Get best score
    max_score = np.max(score_list[:, 0])
    rows_with_max_col0 = score_list[score_list[:, 0] == max_score]
    max_col = rows_with_max_col0[np.argmax(rows_with_max_col0[:, 1])]

    return max_col[0], int(max_col[1])


def chain_np_weighted(anchors_np):
    anchors_len = anchors_np.shape[0]
    if anchors_len == 0:
        return 0

    # Sort anchors by second value ascending, first value descending (tiebreak)
    # anchors_np[:,1] is second col, anchors_np[:,0] is first col
    order = np.lexsort((-anchors_np[:,0], anchors_np[:,1]))
    sorted_anchors = anchors_np[order]

    # Extract first values (a)
    anchors_a = sorted_anchors[:, 0]

    weights = sorted_anchors[:, 2]

    # Compute order with reverse tie-breaks
    order_for_lis = argsort_reverse_ties_np(anchors_a)

    weights = sorted_anchors[:, 2][order_for_lis]

    # length of longest increasing subsequence
    return weighted_LIS(order_for_lis, weights)





def chain_local_np_weighted(anchors_np, mismatch, gap):
    anchors_len = anchors_np.shape[0]
    if anchors_len == 0:
        return 0, 0  # match score 0, length 0
    
    # Sort anchors by second value ascending, first value descending (tiebreak)
    order = np.lexsort((-anchors_np[:, 0], anchors_np[:, 1]))
    sorted_anchors = anchors_np[order]
    
    # Extract first values (a)
    anchors_a = sorted_anchors[:, 0]
    
    # Get order with reverse tie-break
    order = argsort_reverse_ties_np(anchors_a)

    # basically a tuple, contains for each element (score, length)
    score_list = np.zeros((anchors_len, 2))

    for mem_idx in order :
        max_score = sorted_anchors[mem_idx][2]
        max_length = 1
        for j in range(mem_idx) :
            # Non empty
            if score_list[j][0] != 0 :
                # score is score + gap (start - start - end + end) + mismatch (min(start - start, end - end) - 1)
                score_with_j = sorted_anchors[mem_idx][2] + score_list[j][0] + gap * abs(sorted_anchors[mem_idx][0] - sorted_anchors[j][0] - sorted_anchors[mem_idx][1] + sorted_anchors[j][1]) + mismatch * (min(sorted_anchors[mem_idx][0] - sorted_anchors[j][0], sorted_anchors[mem_idx][1] - sorted_anchors[j][1]) - 1)
                
                # If we found a new best
                if score_with_j >= max_score :
                    max_score = score_with_j
                    max_length = max(max_length, score_list[j][1] + 1)
                
        score_list[mem_idx][0] = max_score
        score_list[mem_idx][1] = max_length
    
    # Get best score
    max_score = np.max(score_list[:, 0])
    rows_with_max_col0 = score_list[score_list[:, 0] == max_score]
    max_col = rows_with_max_col0[np.argmax(rows_with_max_col0[:, 1])]

    return max_col[0], int(max_col[1])


