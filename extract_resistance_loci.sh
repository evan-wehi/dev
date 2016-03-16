#!/bin/bash
# extract_resistance_loci.sh
# Take stable list of PlasmoDB gene IDs for resistance loci
# At the moment we have included
# ama1, msp1, msp2, eba175, csp, trap, glurp, kelch13, pfmdr, pfcrt, pfdhfr
# pfdhps, pdx1.
# We extract gene positions and create a BED file and interval_list
# for downstream analysis.
# Usage: ./extract_resistance_loci.sh 

input_gff=data/annotations/Pfalciparum.gff3 
input_loci=$(cat data/annotations/drug_resistance_geneIDs.txt)
input_dict=data/fasta/Pfalciparum.genome.dict
output_file=data/annotations/drug_resistance_genes.txt
output_bed=data/annotations/drug_resistance_genes.bed
output_interval_list=data/annotations/drug_resistance_genes.interval_list

# use awk to subset gff3 to just gene features
# first remove all headers and fasta sequences
egrep -v "^[acgt]|^#|^>" $input_gff | awk -F"\t" '{if($3 == "gene") print $0}' > justGenes.txt
# then scan through file and grab the relevant IDs
for id in ${input_loci}; do
	grep -e $id justGenes.txt
done > $output_file

rm justGenes.txt

# construct BED and interval_list files, 
# first for BED file need to subtract position by 1 since BED is 0 based
awk -F'\t' '{$4=$4-1; $5=$5-1; split($9, a, ";"); sub("ID=","", a[1]); print $1,$4,$5,a[1],$6, $7;}' $output_file > $output_bed

# construct interval_list file (except 1-based) and tab-separated
cat $input_dict > $output_interval_list
awk -F'\t' 'BEGIN {OFS = FS} {split($9, a, ";"); sub("ID=","", a[1]); print $1,$4,$5,$7, a[1]}' $output_file >> $output_interval_list


