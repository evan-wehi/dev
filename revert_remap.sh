#!/bin/bash
# revert_remap.bam
# Description: restore BAM file to previous state and output fastq per readgroup for mapping
# create merged BAM output
# Usage: ./revert_remap.bam input.bam
set -e

PICARD_HOME=/usr/local/bioinfsoftware/picard-tools/current/jars
bwaIndex=./data/fasta/Pfalciparum.genome
# extract SM and RG tags
bam=$1
SM=$(basename $1 .bam)
echo $SM
rgline=$(samtools view -H $1 | grep ^@RG )
echo $rgline

if [[ ! -d data/fastq/${SM} ]]; then
  # fastq files
  mkdir -p data/fastq/${SM}
fi

# revert SAM file, mark adapters then clip them
java -Xmx8g -Djava.io.tmpdir=$TMPDIR  -jar $PICARD_HOME/RevertSam.jar \
VALIDATION_STRINGENCY=SILENT \
INPUT=$bam \
OUTPUT=/dev/stdout \
SORT_ORDER=queryname \
COMPRESSION_LEVEL=0 | java -Xmx8g -jar $PICARD_HOME/MarkIlluminaAdapters.jar \
INPUT=/dev/stdin \
OUTPUT=/dev/stdout \
COMPRESSION_LEVEL=0 \
M=${SM}.adapters.txt | java -Xmx8g  -jar $PICARD_HOME/SamToFastq.jar \
INPUT=/dev/stdin \
CLIPPING_ATTRIBUTE=XT \
CLIPPING_ACTION=2 \
OUTPUT_PER_RG=true \
OUTPUT_DIR=data/fastq/${SM} \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=$TMPDIR


# next remap using bwa mem
# find fastq files and sort by name 
# order should be read_1.fq read_2.fa
r1=$(find ./data/fastq/${SM} -name "*_1.fastq" | sort -u)

for p1 in $r1 
do
  echo $p1
  id=$(basename $p1 _1.fastq)
  rg=$(echo -e $rgline | grep -e ID:$id | sed 's/ /\\t/g')
  echo "Aligning SM: ${SM} for ID:${id}"
  p2=$(find ./data/fastq/${SM} -name "${id}_2.fastq")
  echo "RG line: ${rg}"
  # if paired-end
  if [ ! -z $p2 ] 
  then
    echo "Input fastq files: ${p1}, ${p2}"
    echo "Running bwa mem with paired end alignments"
    bwa mem -M -t 8 -R "${rg}" $bwaIndex $p1 $p2 | java -Xmx8g -jar $PICARD_HOME/CleanSam.jar \
    INPUT=/dev/stdin \
    OUTPUT=/dev/stdout \
    VALIDATION_STRINGENCY=SILENT | java -Xmx8g -jar $PICARD_HOME/SortSam.jar \
    INPUT=/dev/stdin \
    OUTPUT=data/fastq/${SM}/${SM}.${id}.sort.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=FALSE \
    TMP_DIR=$TMPDIR
  else
    echo "Running bwa mem with single-end alignments"
    bwa mem -M -t 8 -R $rg $bwaIndex $p1 | java -Xmx8g -jar $PICARD_HOME/CleanSam.jar \
    INPUT=/dev/stdin \
    OUTPUT=/dev/stdout \
    VALIDATION_STRINGENCY=SILENT | java -Xmx8g -jar $PICARD_HOME/SortSam.jar \
    INPUT=/dev/stdin \
    OUTPUT=data/fastq/${SM}/${SM}.${id}.sort.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=FALSE \
    TMP_DIR=$TMPDIR
  fi

done

# merge bams and output to aligned_bams directory
sort_bams=$(find data/fastq/${SM} -name "${SM}*.sort.bam" -print0 | xargs -0 -I {} echo "INPUT={} ")

java -Xmx8g -Djava.io.tmpdir=$TMPDIR  -jar $PICARD_HOME/MergeSamFiles.jar \
  ${sort_bams} \
  USE_THREADING=true \
  VALIDATION_STRINGENCY=LENIENT \
  AS=true \
  SORT_ORDER=coordinate \
  OUTPUT=data/aligned_bams/${SM}.merged.bam

samtools index data/aligned_bams/${SM}.merged.bam