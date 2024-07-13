
import argparse

# python presegmentation.py -i CpG_coordinates_GCA_000001405.15_GRCh38_no_alt_analysis_set_chr1_chr2.fna.bed -t test_target.bed

def presegment(fasta_infile_obj, target_outfile_obj, dist, min_ncpgs):
    line = fasta_infile_obj.readline().split("\t")
    i, cnt = 1, 1
    chr, l, r = line[:3]
    
    while True:
        line = fasta_infile_obj.readline()
        if not line:
            break
        curr_chr, curr_l, curr_r = line.split("\t")[:3]
        if (int(curr_l) - int(r) > dist - 1) or (curr_chr != chr):
            if cnt >= min_ncpgs:
                target_outfile_obj.write(chr + "\t" + l + "\t" + r + "\t" + "w" + str(i) + "_" + str(cnt) + "\n")
                i+=1
            cnt = 1
            chr, l, r = curr_chr, curr_l, curr_r
        else:
            cnt += 1
            r = curr_r 
    if cnt >= min_ncpgs:
        target_outfile_obj.write(chr + "\t" + l + "\t" + r + "\t" + "w" + str(i) + "_" + str(cnt) + "\n")
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", metavar='', required=True, help="input fasta file", type=str)
    parser.add_argument("-t", metavar='', required=True, help="output targets file", type=str)
    parser.add_argument("-d", metavar='',  required=False, default=300, help="d parameter", type=int)
    parser.add_argument("-mc", metavar='',  required=False, default=5, help="min count of CpGs to create a segment", type=int)
    args = parser.parse_args()
    fasta_infile_obj = open(args.i, "r")
    target_outfile_obj = open(args.t, "w")

    presegment(fasta_infile_obj, target_outfile_obj, args.d, args.mc)
    fasta_infile_obj.close()
    target_outfile_obj.close()




