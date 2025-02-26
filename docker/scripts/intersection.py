import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import logging
from scipy.stats import ttest_ind
import argparse

def read_segments_and_expand(filenames_ls):
    segments_across_samples_over_target_map = {}
    max_end_offset_per_chrom = {}
    target_to_chrom_mapper = {}
    target_order = []
    
    for fn in filenames_ls:
        df_i = pd.read_csv(fn, sep="\t", header=None)
        print(f"Reading file: {fn}")
        print(df_i.head())
        
        for j in range(len(df_i)):
            target_name = df_i[3][j].split('_')[0]
            if target_name not in segments_across_samples_over_target_map:
                target_order.append(target_name)
                target_to_chrom_mapper[target_name] = df_i[0][j]
                segments_across_samples_over_target_map[target_name] = []
                max_end_offset_per_chrom[target_name] = 0
            
            segments_across_samples_over_target_map[target_name].append(int(df_i[1][j]))
            max_end_offset_per_chrom[target_name] = max(max_end_offset_per_chrom[target_name], int(df_i[2][j]))
    
    for chrom, segments in segments_across_samples_over_target_map.items():
        segments_across_samples_over_target_map[chrom] = sorted(list(set(segments)))
        
    return segments_across_samples_over_target_map, max_end_offset_per_chrom, target_to_chrom_mapper, target_order

def write(outfile, segments_across_samples_over_target_map, max_end_offset_per_chrom, target_to_chrom_mapper, target_order):
    with open(outfile, "w") as intersection_seg_off_file_obj:
        for target_name in target_order:
            segments = segments_across_samples_over_target_map[target_name]
            for seg_idx in range(len(segments)):
                seg_beg_offset = segments[seg_idx]
                seg_end_offset = segments[seg_idx + 1] if (seg_idx < (len(segments) - 1)) else max_end_offset_per_chrom[target_name]
                intersection_seg_off_file_obj.write(f"{target_to_chrom_mapper[target_name]}\t{seg_beg_offset}\t{seg_end_offset}\t{target_name}_{seg_idx+1}\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process BED files and generate intersected segment boundaries.")
    
    # Argument for multiple input BED files
    parser.add_argument(
        "-i", "--infiles",
        nargs="+",  # Accepts multiple files
        required=True,
        help="List of input BED files (space-separated)."
    )
    
    # Argument for output file
    parser.add_argument(
        "-o", "--outfile",
        required=True,
        help="Output BED file name."
    )

    args = parser.parse_args()

    # Read and process input files
    segments_across_samples_over_target_map, max_end_offset_per_chrom, target_to_chrom_mapper, target_order = read_segments_and_expand(args.infiles)

    # Write to output file
    write(args.outfile, segments_across_samples_over_target_map, max_end_offset_per_chrom, target_to_chrom_mapper, target_order)

    print(f"Processed {len(args.infiles)} BED files and saved output to {args.outfile}")
