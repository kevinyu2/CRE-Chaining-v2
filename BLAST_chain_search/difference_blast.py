from collections import defaultdict

########## THE FIRST FILE IS THE ONE WE INCLUDE, THE SECOND IS THE ONE WE EXCLUDE

# Which set of ACRs was used and which set are we looking for
ACRS_TO_FIND = '/home/mwarr/Data/setb.txt'
ACRS_USED = '/home/mwarr/Data/seta.txt'
FULL_GENOME = '/home/projects/msu_nsf_pangenomics/pgrp/data/arabidopsis/atacseq/tair10.fa'
chrs_to_use = {'Chr1', 'Chr2', 'Chr3', 'Chr4', 'Chr5'}

# Set to further give a score cutoff or e-value cutoff
CUTOFF1 = 180

# CRE for the CRE version, NCBI for the NCBI blast version. The outputs are a bit different
MODE1 = 'CRE'

# Set to further give a score cutoff or e-value cutoff
CUTOFF2 = 0

# CRE for the CRE version, NCBI for the NCBI blast version. The outputs are a bit different
MODE2 = 'NCBI'

##################################################################################


class BooleanSegmentTree:
    def __init__(self, n):
        self.n = n
        self.count = [0] * (4 * n)  # how many True in the node's segment
        self.lazy_set = [False] * (4 * n)   # marks segments to be set True
        self.lazy_clear = [False] * (4 * n) # marks segments to be set False

    def _push(self, idx, start, end):
        if self.lazy_set[idx]:
            self.count[idx] = end - start + 1
            if start != end:
                self.lazy_set[2*idx] = True
                self.lazy_set[2*idx+1] = True
                self.lazy_clear[2*idx] = False
                self.lazy_clear[2*idx+1] = False
            self.lazy_set[idx] = False

        if self.lazy_clear[idx]:
            self.count[idx] = 0
            if start != end:
                self.lazy_clear[2*idx] = True
                self.lazy_clear[2*idx+1] = True
                self.lazy_set[2*idx] = False
                self.lazy_set[2*idx+1] = False
            self.lazy_clear[idx] = False


    def _update(self, idx, start, end, l, r):
        self._push(idx, start, end)
        if r < start or l > end:
            return
        if l <= start and end <= r:
            # Mark this node lazy True
            self.lazy_set[idx] = True
            self._push(idx, start, end)
            return
        mid = (start + end) // 2
        self._update(2*idx, start, mid, l, r)
        self._update(2*idx+1, mid+1, end, l, r)
        self.count[idx] = self.count[2*idx] + self.count[2*idx+1]

    def _clear_update(self, idx, start, end, l, r):
        self._push(idx, start, end)
        if r < start or l > end:
            return
        if l <= start and end <= r:
            self.lazy_clear[idx] = True
            self._push(idx, start, end)
            return
        mid = (start + end) // 2
        self._clear_update(2*idx, start, mid, l, r)
        self._clear_update(2*idx+1, mid+1, end, l, r)
        self.count[idx] = self.count[2*idx] + self.count[2*idx+1]

    def update_neg(self, l, r):
        # Set range [l,r] to False
        self._clear_update(1, 0, self.n - 1, l, r)

    def _query(self, idx, start, end, l, r):
        self._push(idx, start, end)
        if r < start or l > end:
            return 0
        if l <= start and end <= r:
            return self.count[idx]
        mid = (start + end) // 2
        return self._query(2*idx, start, mid, l, r) + self._query(2*idx+1, mid+1, end, l, r)

    def update(self, l, r):
        # Set range [l,r] to True
        self._update(1, 0, self.n-1, l, r)

    def query(self, l, r):
        # Count how many True in [l,r]
        return self._query(1, 0, self.n-1, l, r)
    



# Gets a dict matching each chr to it's total length
def get_full_chr_lengths() :
    # {chr: length}
    full_chr_lengths = defaultdict(int)
    curr_chr = ""
    with open(FULL_GENOME, 'r') as full :
        for line in full :
            if '>' in line :
                curr_chr = line.rstrip().split('>')[1]
            else :
                if curr_chr in chrs_to_use :
                    full_chr_lengths[curr_chr] += len(line.rstrip())

    return full_chr_lengths

# Hits file is blast result, full_chr_lengths is the dict from above function, k is a threshold for how much of an ACR is covered for it to be "found"
def coverage(hits_file1, hits_file2, full_chr_lengths, k) :
    bst_dict = {}
    for chr, length in full_chr_lengths.items() :
        bst_dict[chr] = BooleanSegmentTree(length)



    # Fill in where is covered
    if MODE1 == 'CRE' :
        with open(hits_file1, 'r') as hits :
            # Skip first two lines
            next(hits)
            next(hits)
            for line in hits :
                line_arr = line.split('\t')
                if CUTOFF1 == None or int(line_arr[1]) >= CUTOFF1 :
                    if line_arr[2] in chrs_to_use :
                        bst_dict[line_arr[2]].update(int(line_arr[3]), int(line_arr[4]))

    else :
        with open(hits_file1, 'r') as hits :
            for line in hits :
                line_arr = line.split('\t')
                if CUTOFF1 == None or float(line_arr[10]) <= CUTOFF1 :
                    if line_arr[0] in chrs_to_use :
                        bst_dict[line_arr[0]].update(int(line_arr[6]), int(line_arr[7]))

    
    # Get rid of in where is covered
    if MODE2 == 'CRE' :
        with open(hits_file2, 'r') as hits :
            # Skip first two lines
            next(hits)
            next(hits)
            for line in hits :
                line_arr = line.split('\t')
                if CUTOFF2 == None or int(line_arr[1]) >= CUTOFF2 :
                    if line_arr[2] in chrs_to_use :
                        bst_dict[line_arr[2]].update_neg(int(line_arr[3]), int(line_arr[4]))

    else :
        with open(hits_file2, 'r') as hits :
            for line in hits :
                line_arr = line.split('\t')
                if CUTOFF2 == None or float(line_arr[10]) <= CUTOFF2 :
                    if line_arr[0] in chrs_to_use :
                        bst_dict[line_arr[0]].update_neg(int(line_arr[6]), int(line_arr[7]))



    # We get the total amount of the acr space covered, and the average coverage per ACR
    total_acr_coverage = 0
    total_acr_length = 0
    total_acr_coverage_percent = 0
    num_acr_greater_than_k = 0
    num_acrs = 0
    with open(ACRS_TO_FIND, 'r') as to_find :
        for line in to_find :
            num_acrs += 1
            chr = line.split('_')[0]
            start = int(line.rstrip().split('_')[1].split('to')[0])
            stop = int(line.rstrip().split('_')[1].split('to')[1])

            acr_len = stop - start + 1
            total_acr_length += acr_len

            # Amount that is covered
            local_coverage =bst_dict[chr].query(start, stop)

            total_acr_coverage += local_coverage
            total_acr_coverage_percent += local_coverage/acr_len
            if (local_coverage/acr_len) > k :
                num_acr_greater_than_k += 1


    print(f"Average Coverage (Total): {total_acr_coverage / total_acr_length}")
    print(f"Average Coverage (Percent): {total_acr_coverage_percent / num_acrs}")
    print(f"Percent Covered Above {k}: {100 * num_acr_greater_than_k / num_acrs}")

    # Calculate how much to discount
    total_repeat = 0
    with open(ACRS_USED, 'r') as used :
        for line in used :
            start = int(line.rstrip().split('_')[1].split('to')[0])
            stop = int(line.rstrip().split('_')[1].split('to')[1])

            acr_len = stop - start + 1
            total_repeat += bst_dict[chr].query(start, stop)

    total_hit_bases = 0
    for chr, bst in bst_dict.items() :
        total_hit_bases += bst.query(0, full_chr_lengths[chr])

    print(f"True Positive Bases: {total_acr_coverage}")
    print(f"False Positive Bases: {total_hit_bases - total_repeat - total_acr_coverage}")

#################################################################################

fcl = get_full_chr_lengths()
# coverage('/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/hits_511_w2_all.tsv', '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/hits_ncbi_all.tsv', fcl, 0.01)
coverage('/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/hits_ArabidopsisPBM.tsv', '/home/projects/msu_nsf_pangenomics/pgrp/dACRxgenomes/one_genome/fimo_full_genome/hits_ncbi_all.tsv', fcl, 0.01)

