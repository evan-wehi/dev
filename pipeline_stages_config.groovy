////////////////////////////////////////////////////////
//
// Pipeline stage definitions for example WGS variant calling 
// pipeline. See pipeline.groovy for more information.
// 
////////////////////////////////////////////////////////

fastqc = {
    doc "Run FASTQC to generate QC metrics for the reads"
    output.dir = "fastqc"
    transform('.fastq.gz')  to('_fastqc.zip')  {
        exec "fastqc --quiet -o ${output.dir} $inputs.gz"
    }
}

createBWAIndex = {
    doc "Run BWA index on fasta if not available."
    output.dir="$REFBASE/indexes"
    exec "bwa index -a bwtsw -p Pfalciparum.genome $REF"
    
}

// Malaria Gen Consortium uses bwa mem
alignBWAmem = {
    doc "Aligns using BWA mem. Note: assumes input file are zipped"
    output.dir="$REFBASE/aligned_bams"
    exec "bwa mem -t 8 -M $REFBASE/indexes  $input.gz > $output.sam"
}

@transform("sai")
alignBWA = {
    doc "Aligns using BWA. Note: assumes input file are gzipped"
    output.dir="align"
    exec "bwa aln -t 8 $encoding_flag $REF $input.gz > $output.sai"
}

@transform("sam")
alignToSamPE = {
    doc "Create SAM files from BWA alignment. Note that these will be very large."
    output.dir="align"
    branch.lane = (input.sai =~ /.*L([0-9]*)_*R.*/)[0][1].toInteger()
    branch.sample = branch.name
    exec """
        bwa sampe $REF -r "@RG\\tID:1\\tPL:$PLATFORM\\tPU:${branch.lane}\\tSM:${branch.sample}"  $input1.sai $input2.sai $input2.gz $input2.gz > $output.sam
    """
}

@transform("bam")
samToSortedBam = {
    doc "Sort a SAM file so that it is compatible with reference order and convert to BAM file"
    output.dir="align"
    exec """
        java -Xmx2g -Djava.io.tmpdir=$TMPDIR  -jar $PICARD_HOME/SortSam.jar 
                    VALIDATION_STRINGENCY=LENIENT 
                    INPUT=$input.sam 
                    OUTPUT=$output.bam 
                    SORT_ORDER=coordinate
    """
}

@filter("merge")
mergeBams = {
    doc "Merge BAM files from multiple lanes or samples together. BAM files should have unique sample names and / or read groups"
    exec """
            java -Xmx2g -Djava.io.tmpdir=$TMPDIR  -jar $PICARD_HOME/MergeSamFiles.jar
                ${inputs.bam.split().collect { "INPUT="+it }.join(' ')}
                USE_THREADING=true 
                VALIDATION_STRINGENCY=LENIENT 
                AS=true 
                OUTPUT=$output.bam 
    """
}

indexBam = {
    // A bit of a hack to ensure the output file is expected in the
    // same directory as the input bam, no matter where it is
    output.dir=file(input.bam).absoluteFile.parentFile.absolutePath
    transform("bam") to ("bam.bai") {
        exec "samtools index $input.bam"
    }
    forward input
}

flagstat = {
    exec "samtools flagstat $input.bam > $output"
}

indexVCF = {
    exec "./vcftools_prepare.sh $input.vcf"
}

realignIntervals = {
    // Hard-coded to take 2 known indels files right now
    output.dir="align"
    exec """
        java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $REF -I $input.bam --known $GOLD_STANDARD_INDELS --known $INDELS_100G -log $LOG -o $output.intervals
    """
}

realign = {
    output.dir="align"
    exec """
        java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $REF -I $input.bam -targetIntervals $input.intervals -log $LOG -o $output.bam
    """
}

dedup = {
    output.dir="align"
    exec """
        java -Xmx6g -Djava.io.tmpdir=$TMPDIR -jar $PICARD_HOME/lib/MarkDuplicates.jar
             INPUT=$input.bam 
             REMOVE_DUPLICATES=true 
             VALIDATION_STRINGENCY=LENIENT 
             AS=true 
             METRICS_FILE=$LOG 
             OUTPUT=$output.bam
    """
}

baseQualRecalCount = {
    doc "Recalibrate base qualities in a BAM file so that quality metrics match actual observed error rates"
    output.dir="align"
    exec "java -Xmx12g -jar $GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -I $input.bam -R $REF --knownSites $DBSNP -l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -log $LOG -o $output.counts"
}

baseQualRecalTabulate = {
    doc "Recalibrate base qualities in a BAM file so that quality metrics match actual observed error rates"
    output.dir="align"
    exec "java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar -T PrintReads -I $input.bam -BQSR $input.counts -R $REF -l INFO -log $LOG -o $output"
}

callSNPs = {
    doc "Call SNPs/SNVs using GATK Haplotype Caller"
    output.dir="variants"
    exec """
            java -Xmx12g -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper 
               -nt $threads 
               -R $REF 
               -I $input.bam 
               --dbsnp $DBSNP 
               -stand_call_conf 50.0 -stand_emit_conf 10.0 
               -dcov 1600 
               -l INFO 
               -A AlleleBalance -A DepthOfCoverage -A FisherStrand 
               -glm SNP -log $LOG 
               -o $output.vcf
        """
}

callIndels = {
    doc "Call variants using GATK Unified Genotyper"
    output.dir="variants"
    exec """
        java -Xmx12g -jar $GATK/GenomeAnalysisTK.jar -T UnifiedGenotyper 
             -nt $threads
             -R $REF 
             -I $input.bam 
             --dbsnp $DBSNP 
             -stand_call_conf 50.0 -stand_emit_conf 10.0 
             -dcov 1600 
             -l INFO 
             -A AlleleBalance -A DepthOfCoverage -A FisherStrand 
             -glm INDEL 
             -log $LOG -o $output.vcf
    """
}

@filter("filter")
filterSNPs = {
    // Very minimal hard filters based on GATK recommendations. VQSR is preferable if possible.
    output.dir="variants"
    exec """
        java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration 
             -R $REF 
             --filterExpression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' 
             --filterName 'GATK_MINIMAL_FILTER'
             --variant $input.vcf 
             -log $LOG 
             -o $output.vcf
    """
}

@filter("filter")
filterIndels = {
    doc """
            Filter data using very minimal hard filters based on GATK recommendations. VQSR is preferable if possible.
            If you have 10 or more samples GATK also recommends the filter InbreedingCoeff < -0.8
        """
    output.dir="variants"
    exec """
        java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar -T VariantFiltration 
                    -R $REF 
                    --filterExpression 'QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0' 
                    --filterName 'GATK_MINIMAL_FILTER' -log $LOG 
                    --variant $input.vcf 
                    -o $output.vcf
    """
}

@filter("vep")
annotateEnsembl = {
    doc "Annotate variants using VEP to add Ensemble annotations"
    output.dir="variants"
    exec """
        perl $VEP/variant_effect_predictor.pl --cache --dir ./vep_cache -i $input.vcf --vcf -o $output.vcf -species human --canonical --per_gene --protein --sift=b --polyphen=b > $LOG
    """
}

depthOfCoverage = {
    output.dir="qc"
    transform("bam") to ("sample_statistics","sample_interval_summary") {
        exec """
            java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
                    -T DepthOfCoverage 
                    -R $REF -I $input.bam 
                    -omitBaseOutput 
                    -ct 1 -ct 10 -ct 20 -ct 30 
                    -o $output.prefix
        """
    }
}
