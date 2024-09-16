# SegMeth
This tool segments genomic regions by CpG methylation values (output of modkit pileup) using circular binary segmentation. CpG methylation values range from 0(unmethylated) to 100 (methylated)

# Overview:
1. first, we create a methylation values array of all CpG sites from the input modkit. Here, we disregard spatial distances between CpGs, and only preserve their ordering. Next, we handle CpG sites having low read coverage (less than threshold). These low-confidence CpG sites can be handled in many ways (eg. estimate missing data using beta distribution; split the regions around the low-confidence CpG sites; filter them out and merge their adjacent segments, etc.). But all the mentioned approaches modify or drop the low-confidence data. Rather, we adopt a simply transformation to segment low-confidence regions as a separate segment class. Specificaly, the methylation value of each low-coverage CpG site is multiplied by a negative one, and then shifted by another -100. Additionally, we also employ a custom segment mean calculation, where we omit either positive methylation values or negative ones (depending on whose relative frequency is higher in the segment) from the segment sum. For example, consider a segment with transformed methylation values [60, 65, -70, 80, 65, -50]. This segment has 4 positive values, and 2 negative ones. So, omitting the negative values, we get segment_sum = 60 + 65 + 80 + 65 = 270, and segment_mean = 270 / 6 = 45. This modification helps ensure that our transformation doesn't affect segment_mean (by +ve and -ve values cancelling each other). This strategy of handling low-coverage CpGs achieves different segments for high and low-coverage regions.

2. Next we run the circular binary segmentation algorithm on each region's methylation array [cite it]. This is a recursive algorithm that scans the array to determine the change points with maximal mean difference (wiz. the difference of the CpG-wise mean methylation levels in the region). The statistical significance of change points in then evaluated by thresholding p-value of the segment (determined using t-test) against a user-defined percentile value (95 by default). The resulting segments are then further split and fed to CBS recursively. There are 2 stopping conditions:- i.) if the segment length drops below 5 CpGs, and ii.) when the subsegment between the detected change-points isn't statistically significant in terms of mean difference signal.

3. Now that we have partitioned each region into segments of methylated regions using the circular method, we validate these segments using the original linear test on each consecutive triple i_{r}, i_{r+1} and i_{r+2} to confirm that i_{r+1} is in fact a stat-sig change point. For this test, we randomly shuffle the data in i_{r} to i_{r+2} some no. of times (say 1000), and compute the t-statistic of the split (i_{r}, i_{r+1}, i_{r+2}) for this random data. If the measured statistic lies in the threshold p-value (default=95th percentile) of the shuffled data, we declare it significant. Else, we drop the break point, and continue.

4. After finding final partition split for each region, in this step, we label the segments into low-coverage, unmethylated, methylated, and uncategorized. Currently, we employ a simple threshold based labelling, with (< 0) as the low-coverage threshold; (0 < methylation mean < unmeth_thresh) as unmethylated segment; and (meth_thresh < methylation mean) as the methylated segment of CpGs. The remaining segments with mean b/w unmeth_thresh and meth_thresh are usually not very valuable, and hence fall within the uncategorized label.
Since we are only concerned with these labels of methylation, this also allows us to merge adjacent segments with the same label into one. This compression also modifies our segment mean values. We create the segmentation output file using these merged segment boundaries and segment mean values.

For delta methylation (hp1 - hp2), we make the following changes to the above steps:-
1. The methylation array is created by using the difference between methylation values between haplotype1 and haplotype2, per CpG site. We filter out CpG sites that have low coverage in atleast one of the haplotypes, and preserve only the valid CpGs in both hps.
As such, methylation array values range between [-100, 100]. 

2. We perform CBS and validation same as before. Next, our threshold-based labelling for delta methylation is the following:- (| methylation mean | > meth_thresh) is labelled as allele specific methylated, elsewhere, we label it as unmethylated. The remaining segment merging remains the same, and the segmentation output file is created.


# Usage 

A command line option is to run the following: -

## `haplotype-specific` mode (for detecting "methylated" and "unmethylated" segments per haplotype)-
```py
python scripts/cbs_delta.py -file test/HG002_1.chr20.bed -file2 test/HG002_2.chr20.bed -o HG002_delta_chr20_segments_out.bed -p 1e-3 -s HG002 -tfile test/targets.bed -ut -70 -mt 70 -delta yes
```
## `delta` mode (for detecting "allele-specific methylated" segments) -
```py
# For haplotype-1
python scripts/cbs_hp.py -file test/HG002_1.chr20.bed -o HG002_1_chr20_segments_out.bed -p 1e-3 -s HG002_1 -tfile test/targets.bed -minCG 5
# For haplotype-2
python scripts/cbs_hp.py -file test/HG002_2.chr20.bed -o HG002_2_chr20_segments_out.bed -p 1e-3 -s HG002_2 -tfile test/targets.bed -minCG 5
```

# SegMeth WDL
SegMeth is also available as a WDL on [dockstore](a link).

## Docker container
The `runsegmeth` task of the workflow uses the [quay.io/repository/shnegi/segmeth](quay.io/repository/shnegi/segmeth). It includes both haplotype-specific and delta segmentation scripts (available at `/opt/scripts/SegMeth-v1.0/cbs_hp.py` and `/opt/scripts/SegMeth-v1.0/cbs_delta.py` respectively in the container). It also contains the segment intersection (`/opt/scripts/SegMeth-v1.0/segment_intersection.py`) script to intersect segments across multiple samples to generate a consistent set for performing DMR analysis.

# Input Parameters

## `haplotype-specific` mode
```
python3 scripts/cbs_hp.py --help
usage: cbs_hp.py [-h] -file  -o  [-p] -s  -tfile  [-ut] [-mt] [-minCG] [-filt_thresh] [-plot]

    -file,  input mokdit bedMethyl pileup
    -o, output file name, to be populated with intersected segments
    -s, sample name
    -tfile, target file name for generating segments in
    -p, OPTIONAL. p-value to be used for significance testing of intermediate segments in CBS algorithm. Increase p-value to get finer segments, and decrease p-value for coarser segments. DEFAULT: 0.5
    -ut, OPTIONAL. threshold for annotating a segment as "unmethylated". DEFAULT: 30
    -mt, OPTIONAL. threshold for annotating a segment as "methylated". DEFAULT: 70
    -filt_thresh, OPTIONAL. valid coverage threshold for rejecting/selecting valid CpG sites. DEFAULT: 5
    -plot, OPTIONAL. make segmentation plots per region? [yes/no], DEFAULT: no
```

## `delta` mode
```
python3 scripts/cbs_delta.py --help
usage: cbs_delta.py [-h] -file  [-file2] -o  [-p] -s  -tfile  [-ut] [-mt] [-minCG] [-filt_thresh] [-plot] [-delta]

    -file,  input mokdit bedMethyl pileup for haplotype-1
    -file2,  input mokdit bedMethyl pileup for haplotype-2
    -o, output file name, to be populated with intersected segments
    -s, sample name
    -tfile, target file name for generating segments in
    -p, OPTIONAL. p-value to be used for significance testing of intermediate segments in CBS algorithm. Increase p-value to get finer segments, and decrease p-value for coarser segments. DEFAULT: 0.5
    -ut, OPTIONAL. threshold for annotating a segment as "allele-specific methylated with hp2-hp1 >= 70". DEFAULT: -70
    -mt, OPTIONAL. threshold for annotating a segment as "allele-specific methylated with hp1-hp2 >= 70". DEFAULT: 70
    -filt_thresh, OPTIONAL. valid coverage threshold for rejecting/selecting valid CpG sites. Only CpG sites with valid coverage in both haplotypes will be included for segmentation. DEFAULT: 5
    -plot, OPTIONAL. make segmentation plots per region? [yes/no], DEFAULT: no
    -delta, OPTIONAL. Run on delta mode? [YES], or haplotype-specific mode? [NO], DEFAULT: no
```

# Test locally
```sh
## Run with miniwdl
miniwdl run --as-me -i test.inputs.json wdl/workflow.wdl

## Test with Toil


```
# Outputs

Outputs a bed file containing info about segments in consecutive lines. Each line represents one segment, with the following parameters:-

| column | name                  | description                                                                    | type  |
|--------|-----------------------|--------------------------------------------------------------------------------|-------|
| 1      | chrom                 | name of reference sequence from BAM header                                     | str   |
| 2      | start position        | 0-based start position                                                         | int   |
| 3      | end position          | 0-based exclusive end position                                                 | int   |
| 4      | target region    | Name of target region, underscore, and index of the current segment (in inc. order of position) in that target region.                                           | str   |
| 5      | score                 | included for compatibility                                                | int   |
| 6      | strand                | '+' for positive strand '-' for negative strand, '.' when strands are combined | str   |
| 7      | start position        | included for compatibility                                                     | int   |
| 8      | end position          | included for compatibility                                                     | int   |
| 9      | color                 | Color (RGB) representing label of the segment. Color code used for methylated and unmethylated segmentation:- unmethylated -> [0,0,255], methylated -> [255,0,0], uncategorized -> [192,192,192], low-coverage -> [0,255,0]; for allele-specific methylation segmentation:- hypomethylated -> [0,76,153], hypermethylated -> [152,0,76], uncategorized -> [192,192,192]                                     | str   |
| 10     | Mean methylation value of segment | calculated using custom mean calculation from valid CpGs; see definitions above                                                         | float   |
| 11     | N<sub>CpG</sub>     | Number of CpGs in segment                                        | int |
| 12     | N<sub>valid-CpG</sub>       | Number of valid CpGs used in mean calculation; see definitions above                                                         | int   |