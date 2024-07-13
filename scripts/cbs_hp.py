import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import logging
from scipy.stats import ttest_ind
import argparse  


log = logging.getLogger()
logging.basicConfig(level=logging.WARN)

def calculate_segment_mean(x):
    '''
    Custom function to calculate mean of segment handling low coverage sites
    '''
    neg_cnt, neg_sum, pos_sum = 0, 0, 0
    for ele in x:
        neg_cnt += 1 if (ele<0) else 0
        neg_sum += ele if (ele<0) else 0
        pos_sum += ele if (ele>0) else 0
    if neg_cnt>len(x)/2:
        return float(neg_sum) / neg_cnt
    return float(pos_sum) / (len(x) - neg_cnt)


def transform_pos_neg_array(x):
    '''
    Custom function to remove either negatives, or positives from an array, depending on their frequency. 
    '''
    res=[]
    neg_cnt = sum(1 for ele in x if ele < 0)
    if neg_cnt > len(x)/2:
        res = [ele for ele in x if ele<0]
    else:
        res = [ele for ele in x if ele>0]
    return np.array(res)


def calc_segment_mean_from_metadata(neg_sum, neg_cnt, pos_sum, pos_cnt):
    if neg_cnt > pos_cnt:
        return neg_sum / neg_cnt
    return pos_sum / pos_cnt


def cbs_stat(x, x_pos, min_cpgs=5):
    '''Given x, Compute the subinterval x[i0:i1] with the maximal segmentation statistic t. 
    Returns t, i0, i1'''
    i0, i1 = 0, 0
    n = len(x)
    max_abs_mean_diff_yet = 0
    truth_neg_cnt_ext = sum(1 for ele in x if ele<0)
    truth_pos_cnt_ext = len(x) - truth_neg_cnt_ext
    truth_pos_seg_sum_ext = sum(ele for ele in x if ele>0)
    truth_neg_seg_sum_ext = sum(ele for ele in x if ele<0)

    for i in range(min_cpgs, len(x)-min_cpgs):
        neg_cnt_internal, pos_cnt_internal = 0, 0
        pos_seg_sum_internal, neg_seg_sum_internal = 0, 0
        neg_cnt_ext, pos_cnt_ext = truth_neg_cnt_ext, truth_pos_cnt_ext
        pos_seg_sum_ext, neg_seg_sum_ext = truth_pos_seg_sum_ext, truth_neg_seg_sum_ext
        for j in range(i, len(x)+1-min_cpgs):
            if i==0 and j==len(x):
                continue
            if x[j-1]<0:
                neg_seg_sum_internal += x[j-1]
                neg_seg_sum_ext -= x[j-1]
                neg_cnt_internal += 1
                neg_cnt_ext -= 1
            else:
                pos_seg_sum_internal += x[j-1]
                pos_seg_sum_ext -= x[j-1]
                pos_cnt_ext -= 1
                pos_cnt_internal += 1
            if j-i < min_cpgs:
                continue
            curr_mean_diff = abs(
                calc_segment_mean_from_metadata(neg_seg_sum_internal, neg_cnt_internal, pos_seg_sum_internal, pos_cnt_internal)
                - calc_segment_mean_from_metadata(neg_seg_sum_ext, neg_cnt_ext, pos_seg_sum_ext, pos_cnt_ext)
            )
            if (max_abs_mean_diff_yet < curr_mean_diff or ((max_abs_mean_diff_yet == curr_mean_diff) and ((j-1-i) > (i1-i0)))):
                i0, i1 = i, j-1
                max_abs_mean_diff_yet = curr_mean_diff

    # Perform independent samples t-test
    t_statistic, p_value = ttest_ind(transform_pos_neg_array(x[i0:(i1+1)]), np.concatenate([transform_pos_neg_array(x[i1+1:n]), transform_pos_neg_array(x[0:i0])]))
    return t_statistic, p_value, i0, i1+1


def tstat(x, i):
    '''Return the segmentation statistic t testing if i is a (one-sided) breakpoint in x'''
    n = len(x)
    s0 = calculate_segment_mean(x[:i])
    s1 = calculate_segment_mean(x[i:])
    return (n-i)*i/n*(s0-s1)**2


def cbs(x, x_pos, shuffles=1000, p=.05, min_cpgs=5):
    '''Given x, find the interval x[i0:i1] with maximal segmentation statistic t. Test that statistic against
    given (shuffles) number of random permutations with significance p.  Return True/False, t, i0, i1; True if
    interval is significant, false otherwise.'''

    max_t, p_value, max_start, max_end = cbs_stat(x, x_pos)
    if max_end-max_start == len(x):
        return False, max_t, max_start, max_end
    if max_start < min_cpgs:
        max_start = 0
    if len(x)-max_end < min_cpgs:
        max_end = len(x)
    if p_value > p:
        return False, max_t, max_start, max_end
    return True, max_t, max_start, max_end


def rsegment(x, x_pos, start, end, L_segs=[], shuffles=1000, p=.05, min_cpgs=5):
    '''Recursively segment the interval x[start:end] returning a list L of pairs (i,j) where each (i,j) is a significant segment.
    '''
    threshold, t, s, e = cbs(x[start:end], x_pos[start:end], shuffles=shuffles, p=p, min_cpgs=min_cpgs)
    # log.info('Proposed partition of {} to {} from {} to {} with t value {} is {}'.format(start, end, start+s, start+e, t, threshold))
    # print('Proposed partition of {} to {} from {} to {} with t value {} is {}'.format(x_pos[start], x_pos[end-1] + 1, x_pos[start+s], x_pos[start+e-1] + 1, t, threshold))
    
    if (not threshold) | (e-s < 5) | (e-s == end-start):
        L_segs.append((start, end))
    else:
        if s > 0:
            rsegment(x, x_pos, start, start+s, L_segs, shuffles=shuffles, p=p, min_cpgs=min_cpgs)
        if e-s > 0:
            rsegment(x, x_pos, start+s, start+e, L_segs, shuffles=shuffles, p=p, min_cpgs=min_cpgs)
        if start+e < end:
            rsegment(x, x_pos, start+e, end, L_segs, shuffles=shuffles, p=p, min_cpgs=min_cpgs)
    return L_segs


def segment(x, x_pos, shuffles=1000, p=.05, min_cpgs=5):
    '''Segment the array x, using significance test based on shuffles rearrangements and significance level p
    '''
    start = 0
    end = len(x)
    L_segs = []
    rsegment(x, x_pos, start, end, L_segs, shuffles=shuffles, p=p, min_cpgs=min_cpgs)
    return L_segs


def validate(x, x_pos, L_segs, shuffles=1000, p=.01):
    S_segs = [i[0] for i in L_segs]+[len(x)]
    # S_segs_sites = [sample_positions[i+island_start_cpg_idx] for i in S_segs]
    # print("S_segs: ", S_segs)
    # print("S (before validation): ", S_segs_sites)
    SV = [0]
    left = 0
    for test, s in enumerate(S_segs[:-2]):
        # print("iteration at: ", S_segs_sites[test+1])
        t = tstat(x[S_segs[left]:S_segs[test+2]], S_segs[test+1]-S_segs[left])
        # print('Testing validity of {} in interval from {} to {} yields statistic {}'.format(x_pos[S_segs[test+1]], x_pos[S_segs[left]], S_segs_sites[test+2], t))
        threshold = 0
        thresh_count = 0
        site = S_segs[test+1]-S_segs[left]
        xt = x[S_segs[left]:S_segs[test+2]].copy()
        flag = True
        for k in range(shuffles):
            np.random.shuffle(xt)
            threshold = tstat(xt, site)
            if threshold > t:
                thresh_count += 1
            if thresh_count >= p*shuffles:
                flag = False
                # print('Breakpoint {} rejected'.format(x_pos[S_segs[test+1]]))
                break
        if flag:
            # print('Breakpoint {} accepted'.format(x_pos[S_segs[test+1]]))
            SV.append(S_segs[test+1])
            left = test + 1
    SV.append(S_segs[-1])
    return SV


def read_array_of_methylation_perc(filepath, filter_thresh=5):
    '''
    Read input from sorted bedmethyl file and create an array with sequence of %methylation. 
    Stores methylation values in the array data. 
    Methylation values for low-coverage regions (<3 reads are multiplied by -1. 
    This improves their isolation and segmentation)
    '''
    # print("Reading CpG array for haplotype-1")
    df = pd.read_csv(filepath, sep="\t", header=None)
    chrom_str = df[0][0]
    data_positions = np.array(df.iloc[:, 1], dtype=np.int64)
    filter_condition = df[4] < filter_thresh
    df.loc[filter_condition, 10] = df.loc[filter_condition, 10]*(-1) - 100
    df = df.iloc[:, 10]
    data = np.array(df, dtype=np.float64)
    # print(data[0:5], data_positions[0:5])
    return data, data_positions, chrom_str


def draw_segmented_data(data, S, unmethyl_thresh, methyl_thresh, title=None):
    '''Draw a scatterplot of the data with vertical lines at segment boundaries and horizontal lines at means of 
    the segments. S is a list of segment boundaries.'''
    # Define custom colors based on condition
    colors = ['red' if (methyl_thresh <= y <= 100) else 'blue' if (0 <= y <= unmethyl_thresh) else 'darkgreen' if (y<0) else 'grey' for y in data]
    # print("colors size: ", len(colors))
    j=sns.scatterplot(x=range(len(data)),y=data,color=colors,size=.1,legend=None)
    # print("data size: ", len(data))
    # print("S: ", S)
    for x in S:
        j.axvline(x)
    for i in range(1,len(S)):
        j.hlines(calculate_segment_mean(data[S[i-1]:S[i]]),S[i-1],S[i],color='green')
    j.set_title(title)
    j.set_ylim(-110,110)
    j.get_figure().set_size_inches(16,4)
    return j


def return_methyl_grp(curr_meth_mean, unmeth_thresh, meth_thresh):
    '''
    Returns index of methylation group under which supplied mean categorizes.
    For single hap:-
        low-coverage:- -1, hypo-methylated:- 0, moderate-methylation:- 1, hyper-methylated:- 2
    '''
    if curr_meth_mean < 0:
        return -1
    elif curr_meth_mean <= unmeth_thresh:
        return 0
    elif curr_meth_mean >= meth_thresh:
        return 2
    return 1


def merge_segments(segment_boundaries, data, data_positions, unmeth_thresh, meth_thresh):
    '''
    Merges adjacent segments which fall into the same methylation group.
    '''
    prev_seg_sz, prev_seg_mean, prev_seg_methyl_grp = 0, 0, -2
    merged_segment_boundaries = [0]
    for i in range(1, len(segment_boundaries)):
        curr_seg_mean = calculate_segment_mean(data[segment_boundaries[i-1]:segment_boundaries[i]])
        curr_seg_methyl_grp = return_methyl_grp(curr_seg_mean, unmeth_thresh, meth_thresh)
        # print(" inside loop for end offset: ", data_positions[segment_boundaries[i]-1]+1, " with curr seg mean: ", curr_seg_mean, ", prev seg mean: ", prev_seg_mean)
        if curr_seg_methyl_grp == prev_seg_methyl_grp:
            merged_segment_boundaries.pop()
            merged_segment_boundaries.append(segment_boundaries[i])
            prev_seg_mean = (prev_seg_mean*prev_seg_sz + np.sum(data[segment_boundaries[i-1]:segment_boundaries[i]])) / (prev_seg_sz + (segment_boundaries[i]-segment_boundaries[i-1]))
            prev_seg_sz += (segment_boundaries[i]-segment_boundaries[i-1])
        else:
            # print("reached here for end offset: ", data_positions[segment_boundaries[i]-1]+1, " with curr seg mean: ", curr_seg_mean, ", prev seg mean: ", prev_seg_mean)
            merged_segment_boundaries.append(segment_boundaries[i])
            prev_seg_mean, prev_seg_sz, prev_seg_methyl_grp = curr_seg_mean, (segment_boundaries[i]-segment_boundaries[i-1]), curr_seg_methyl_grp

    return merged_segment_boundaries


if __name__ == '__main__':
    parser = argparse.ArgumentParser()  
  
    # creating two variables using the add_argument method  
    parser.add_argument("-file", metavar='', required=True, help="modkit filtered input file")
    parser.add_argument("-o", metavar='', required=True, help="output file", type=str)
    parser.add_argument("-p", metavar='', required=False, default=0.5, help="p-value, DEFAULT: 0.5", type=float)
    parser.add_argument("-s", metavar='', required=True, help="sample name", type=str)
    #parser.add_argument("-t", metavar='', required=True, help="target name", type=str)
    parser.add_argument("-tfile", metavar='', required=True, help="target file name", type=str)
    parser.add_argument("-ut", metavar='', required=False, default=30, help="unmethylated threshold for merging segments, DEFAULT: 30", type=str)
    parser.add_argument("-mt", metavar='', required=False, default=70, help="methylated threshold for merging segments, DEFAULT: 70", type=str)
    parser.add_argument("-minCG", metavar='', required=False, default=5, help="minimum count of CpG sites per segment, DEFAULT: 5", type=str)
    parser.add_argument("-filt_thresh", metavar='', required=False, default=5, help="filter threshold of read count for rejecting/selecting CPG sites, DEFAULT: 5", type=str)
    parser.add_argument("-plot", metavar='', required=False, default="no", help="make segmentation plots?='yes/no', DEFAULT: no", type=str)

    args = parser.parse_args()

    log.setLevel(logging.INFO)
    # array with delta methylation % from combined
    sample, sample_positions, chrom_str = read_array_of_methylation_perc(filepath=args.file, filter_thresh=args.filt_thresh)
    seg_off_file_obj = open(args.o, "a")
    
    targets_file_path = args.tfile
    targets_df = pd.read_csv(targets_file_path, sep="\t", header=None)
    num_cpg_islands = targets_df.shape[0]
    island_start_cpg_idx, island_end_cpg_idx = 0, 0

    for cpg_idx in range(num_cpg_islands):
        target_name = targets_df[3][cpg_idx]
        print("Running for target {}".format(target_name))
        island_start_pos, island_end_pos = int(targets_df[1][cpg_idx]), int(targets_df[2][cpg_idx])
        while (island_start_cpg_idx < len(sample_positions)) and (sample_positions[island_start_cpg_idx] < island_start_pos):
            island_start_cpg_idx += 1
        island_end_cpg_idx = max(island_end_cpg_idx, island_start_cpg_idx)
        while (island_end_cpg_idx < len(sample_positions)) and (sample_positions[island_end_cpg_idx] <= island_end_pos):
            island_end_cpg_idx += 1
        if island_end_cpg_idx <= island_start_cpg_idx:
            continue

        # print("island number: ", (cpg_idx+1), ", ", sample_positions[island_start_cpg_idx], ", ", sample_positions[island_end_cpg_idx-1] + 1)

        curr_sample = sample[island_start_cpg_idx: island_end_cpg_idx]
        curr_sample_positions = sample_positions[island_start_cpg_idx: island_end_cpg_idx]

        L = segment(curr_sample, curr_sample_positions, shuffles=1000, p=args.p, min_cpgs=int(args.minCG))
        S = validate(curr_sample, curr_sample_positions, L, p=args.p)
        S = merge_segments(S, curr_sample, curr_sample_positions, int(args.ut), int(args.mt))

        for i in range(len(S)-1):
            segment_mean = calculate_segment_mean(curr_sample[S[i]:S[i+1]])
            if segment_mean >= int(args.mt):
                # methylated
                segment_igv_color = [255, 0, 0]
            elif segment_mean <= int(args.ut) and segment_mean >= 0:
                # unmethylated
                segment_igv_color = [0, 0, 255]
            elif segment_mean < 0:
                #low read coverage
                segment_igv_color = [0, 255, 0]
            else:
                segment_igv_color = [192, 192, 192]
            seg_off_file_obj.write(chrom_str + "\t" + str(curr_sample_positions[S[i]]) + "\t" + str(curr_sample_positions[S[i+1]-1]+1) + "\t" + target_name + "_" + str(i+1) + "\t" + "0\t.\t" + str(curr_sample_positions[S[i]]) + "\t" + str(curr_sample_positions[S[i+1]-1]+1) + "\t" + str(segment_igv_color[0]) + "," + str(segment_igv_color[1]) + "," + str(segment_igv_color[2]) + "\t" + str(segment_mean) + "\t" + str(S[i+1]-S[i]) + "\n")
        
        if args.plot == "yes":
            ax = draw_segmented_data(curr_sample,  S, int(args.ut), int(args.mt), title='Circular Binary Segmentation of ' + "cpg" + str(cpg_idx + 1) + ' target region in ' + args.s)
            ax.get_figure().savefig(args.s + "_" + "cpg" + str(cpg_idx + 1) + "_" + chrom_str + "_" + str(curr_sample_positions[S[0]]) + "_" + str(curr_sample_positions[S[-1]-1]+1) + ".png")
    seg_off_file_obj.close()

# Run script:
# python cbs_hp.py -file HG002_R10_2.chr1.bed -o HG002_chr1_test_segments.bed -p 1e-4 -s HG002_hp2 -tfile presegments_HG002_testing.bed -plot yes -minCG 10