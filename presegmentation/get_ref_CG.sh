#!/bin/bash
#SBATCH --job-name=get_CpG_coordinates.sh.20240111
#SBATCH --partition=medium
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=2gb
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --output=get_CpG_coordinates.sh.20240111.%j.log

set -e
#set -x

pwd; hostname; date

fasta_file=$1
formatted_file=${fasta_file}.formatted.fa
CpG_coordinates=CpG_coordinates_${fasta_file}.bed

rm -f ${CpG_coordinates}.tmp #delete previous file with CpG coordinates  
seqtk seq -L 0 ${fasta_file} > ${formatted_file}

# Get a list of chromosome names from the FASTA file
chromosomes=$(grep ">" $formatted_file | sed 's/>//')

# Iterate over each chromosome
for chr in $chromosomes; do
    # Extract the sequence of the current chromosome
    chr_sequence=$(grep -A 1 "$chr" $formatted_file | tail -n 1)

    echo "Processing $chr"
    #write C for the positive strand and G for the negative strand, for each CpG
    echo "$chr_sequence" | awk 'BEGIN {OFS="\t"} {gsub(/[^CG]/, "N"); print}' | grep -o -b 'CG' | awk -v chromosome="$chr" -F: '{print chromosome"\t"$1"\t"$1+1"\tC\t.\t+\n"chromosome"\t"$1+1"\t"$1+2"\tG\t.\t-"}' >>${CpG_coordinates}.tmp
done

sort -k1,1 -k2,2n ${CpG_coordinates}.tmp >${CpG_coordinates}

rm -f ${formatted_file} ${CpG_coordinates}.tmp
echo "Done."