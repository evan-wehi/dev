#!/bin/bash
# Description: Pairwise alignment between 3D7 and Plasmodium Reichenowi genomes.
set -e 

fquery=data/Pfalciparum.fasta
fsearch=data/Preichenowi.fasta

if [ ! -d data/pairwise_alignments ] then;
    mkdir -p data/pairwise_alignments
fi

# pass 1 - usearch global alignment w 90% identity
usearch -usearch-global $fquery -db $fsearch -id 0.9 -alnout data/pairwise_alignments/usearch_global_Pf_Pr.aln

# pass 2 - usearch 
