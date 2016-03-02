#!/bin/bash
# revert_remap.bam
# Description: restore BAM file to previous state and output fastq per readgroup for mapping
# create merged BAM output
# Usage: ./revert_remap.bam input.bam
set -e

PICARD_HOME=/usr/local/bioinfsoftware/picard-tools/current/jars

# extract SM and RG tags
bam=$1
SM=$(basename $1 .bam)
echo $SM
rgline=$(samtools view -H $1 | grep ^@RG )

java -Xmx8g -Djava.io.tmpdir=$TMPDIR  -jar $PICARD_HOME/RevertSam.jar \
  VALIDATION_STRINGENCY=SILENT \
  INPUT=$bam \
  OUTPUT=/dev/stdout \
  SORT_ORDER=coordinate | java -jar $PICARD_HOME/SamToFastq.jar \
  INPUT=/dev/stdin \
  OUTPUT_PER_RG=true \
  OUTPUT_DIR=./fastq/${SM} \
  VALIDATION_STRINGENCY=SILENT


# next remap using bwa mem
# find fastq files that the correspond to  SM tag
r1=$(find ./fastq -name "*_1.fastq" | sort -u)
r2=$(find ./fastq -name "*_2.fastq" | sort -u)


bwa mem -M -r $rgline r1.fastq r2.fastq | 
	java -Xmx8g -jar $PICARD_HOME/CleanSam.jar \
      INPUT=/dev/stdin \
      OUTPUT=/dev/stdout \
      VALIDATION_STRINGENCY=SILENT \
    | java -Xmx8g -jar $PICARD_HOME/SortSam.jar \
      INPUT=/dev/stdin \
      OUTPUT=${SM}.${id}.sort.bam \
      SORT_ORDER=coordinate \
      CREATE_INDEX=FALSE \
      TMP_DIR=$TempDir

# merge bams
sort_bams=$(find ./${SM} -name "${SM}*.sort.bam" -print0 | xargs -0 -I {} echo "INPUT={} ")

java -Xmx8g -Djava.io.tmpdir=$TMPDIR  -jar $PICARD_HOME/MergeSamFiles.jar \
        ${sort_bams} \
        USE_THREADING=true \
        VALIDATION_STRINGENCY=LENIENT \ 
        AS=true \
        SORT_ORDER=coordinate \
        OUTPUT=${SM}.merged.bam

samtools index ${SM}.merged.bam