#!/bin/sh

SOFTWARE="/usr/local/bioinfsoftware"

############
# Trimming #
############
TRIMMOMATIC="java -jar $SOFTWARE/trimmomatic/current/trimmomatic-0.30.jar"
ADAPTORS="ILLUMINACLIP:$SOFTWARE/trimmomatic/Trimmomatic-0.30/adapters/TruSeq3-SE.fa:2:30:10"

echo --------------------------------------------
echo                   Trimming
echo --------------------------------------------
$TRIMMOMATIC SE ev-extract.fastq ev-trimmed.fastq $ADAPTORS

############
# Aligning #
############
BWA="$SOFTWARE/bwa/current/bin/bwa"
REF_SEQ=~/local/share/bcbio/genomes/Hsapiens/GRCh37/bwa/GRCh37.fa

echo Creating index
$BWA index ev-trimmed.fastq

echo --------------------------------------------
echo                  Aligning
echo --------------------------------------------
$BWA mem -R "@RG\tID:Seq01p\tSM:Seq01\tPL:ILLUMINA\tPI:330" $REF_SEQ ev-trimmed.fastq > ev-aligned.sam

##############
# sam -> bam #
##############
echo --------------------------------------------
echo              Sort and index
echo --------------------------------------------
SAM="$SOFTWARE/samtools/samtools-1.3.1/bin/samtools"
$SAM view -b ev-aligned.sam | $SAM sort - > ev-aligned-sorted.bam
$SAM index ev-aligned-sorted.bam ev-aligned-sorted.bai

####################
# Mutation calling #
####################
REF_SEQ=~/local/share/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
MUTECT2="java -jar $SOFTWARE/gatk/current/GenomeAnalysisTK.jar -T MuTect2"

echo --------------------------------------------
echo                 Calling
echo --------------------------------------------
$MUTECT2 --artifact_detection_mode -R $REF_SEQ -I:tumor ev-aligned-sorted.bam -o ev.vcf 