////////////////////////////////////////////////////////
//
// Pipeline stage definitions for example WGS variant calling 
// pipeline. See pipeline.groovy for more information.
// 
////////////////////////////////////////////////////////

createBWAIndex = {
    doc "Run BWA index on fasta if not available."
    output.dir="$REFBASE/indexes"
    exec "bwa index -a bwtsw -p Pfalciparum.genome $REF"
    
}

remapBWAmem = {
    doc "Extract fastq from BAM file and remap using bwa mem"
    output.dir="$REFBASE/aligned_bams"
    exec "sh revert_remamp.sh $input"
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

bqsrPass1 = {
    doc "Recalibrate base qualities in a BAM file so that quality metrics match actual observed error rates"
    output.dir="align"
    exec """
            java -Xmx12g -jar $GATK/GenomeAnalysisTK.jar
                -T BaseRecalibrator 
                -R $REF  
                -I $input.bam  
                -knownSites $REFSNP1  
                -knownSites $REFSNP2 
                -knownSites $REFSNP3
                -o $output.pass1.table 
        """
    
}

bqsrPass2 = {
    doc "Recalibrate base qualities in a BAM file so that quality metrics match actual observed error rates"
    output.dir="align"
    exec """
            java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar 
                -T BaseRecalibrator 
                -R $REF 
                -I $input.bam 
                -knownSites $REFSNP1  
                -knownSites $REFSNP2 
                -knownSites $REFSNP3 
                -BQSR $input.table \ 
                -o $output.pass2.table 
        """
}

bqsrCheck = {
    doc "Compare pre and post base-quality scores from recalibration"
    output.dir="qc"
    exec """
        java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
            -T AnalyzeCovariates
            -R $REF
            -before $input.pass1.table \
            -after $input.pass2.table \
            -plots $ouptut.bqsr.pdf
    """
}

bqsrApply = {
    doc "Apply BQSR to input BAM file."
    output.dir="align"
    exec """
        java -jar $GATK/GenomeAnalysisTK.jar  
            -T PrintReads  
            -R $REF 
            -I $input.bam   
            -BQSR $input.pass1.table  
            -o $output.bam
    """ 
}



callVariants = {
    doc "Call SNPs/SNVs using GATK Haplotype Caller, produces .g.vcf"
    output.dir="variants"
    exec """
        java -Xmx2g -jar $GATK/GenomeAnalysisTK.jar
            -T HaplotypeCaller 
            -R $REF
            -I $input.bam
            --emitRefConfidence GVCF
            -gt_mode DISCOVERY
            -o $output.g.vcf    
        """
}

genotype = {
    doc "Jointly genotype gVCF files"
    output.dir="variants"
    exec """
         gvcf_files=$(find ./variants -name "*.g.vcf" -print0 | xargs -0 -I {} echo "--variant {}")
         java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar
            -T GenotypeGVCFs
            $gvcf_files
            -o all_raw.vcf
    """
}

// variant recalibration
vqsrGenerate = {
    doc "Generate variant quality score recalibration"
    output.dir="variants"
    exec """
         java -Xmx12g -jar $GATK/GenomeAnalysisTK.jar 
            -T VariantRecalibrator 
            -R $REF 
            -input all_raw.vcf
            -resource:cross1,training=true $REFSNP1
            -resource:cross2,training=true $REFSNP2
            -resource:cross2,training=true $REFSNP3
            -an QD
            -an MQ
            -an FS 
            -an SOR 
            -an DP
            -mode SNP 
            -mG 8
            -MQCap 70 
            -recalFile $output.recal 
            -tranchesFile $output.tranches 
            -rscriptFile $output.plots.R

    """
}

vqsrApply = {
    doc "Apply variant quality score recalibration"
    output.dir="variants"
    exec """
        java -jar $GATK/GenomeAnalysisTK.jar 
            -T ApplyRecalibration 
            -R $REF  
            -input all_raw.vcf  
            -mode SNP  
            --ts_filter_level 99.0  
            -recalFile $output.recal  
            -tranchesFile $output.tranches  
            -o all_recalibrated.vcf
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


// QC metrics at BAM level
flagstat = {
    exec "samtools flagstat $input.bam > $output"
}

fastqc = {
    doc "Run FASTQC to generate QC metrics for the reads"
    output.dir = "fastqc"
    transform('.fastq.gz')  to('_fastqc.zip')  {
        exec "fastqc --quiet -o ${output.dir} $inputs.gz"
    }
}


depthOfCoverage = {
    output.dir="qc"
    transform("bam") to ("sample_statistics","sample_interval_summary") {
        exec """
            java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
                    -T DepthOfCoverage 
                    -R $REF -I $input.bam 
                    -omitBaseOutput 
                    -ct 1 -ct 5 -ct 10
                    -o $output.prefix
        """
    }
}
