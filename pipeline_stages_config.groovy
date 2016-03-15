////////////////////////////////////////////////////////
//
// Pipeline stage definitions for example WGS variant calling 
// pipeline. See pipeline.groovy for more information.
// 
////////////////////////////////////////////////////////
createBWAIndex = {
    doc "Run BWA index on fasta if not available."

    exec """
        bwa index -a bwtsw -p Pfalciparum.genome $REF
        mv Pfalciparum.genome.* $REFBASE/fasta/

    """

}

fastaIndex = {
    doc "index fasta file if not available"
    exec "samtools faidx $REF"
}

createDictionary = {
    doc "Create Sequence dictionary"
    exec "java -jar $PICARD_HOME/CreateSequenceDictionary.jar R=$REF O=$REFBASE/fasta/Pfalciparum.genome.dict"
}

// this could probably be broken up into multiple steps
// but I'm still learning bpipe intricacies 
remapBWAmem = {
    doc "Extract fastq from BAM file and remap using bwa mem"
    exec "./revert_remap.sh $input.bam"
}

sampleID = {
    bamfile="$input".replaceAll(".bam", ".merged.bam")
    sm=bamfile.replaceAll(/^.*\/{1}/, '')
    branch.sample="$sm"
}

@transform(".intervals")
realignIntervals = {
    output.dir="$REFBASE/aligned_bams"
    exec """
        java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
          -T RealignerTargetCreator
          -R $REF 
          -I $REFBASE/aligned_bams/$sample
          -log $LOG
          -nt $NTHREAD_GATK 
          -o $output.intervals
    """
}

realignIndels = {
    output.dir="$REFBASE/aligned_bams"
    exec """
        java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar 
          -T IndelRealigner 
          -R $REF 
          -I $REFBASE/aligned_bams/$sample
          -targetIntervals $input.intervals 
          -log $LOG
          -nt $NTHREAD_GATK
          -o $output.bam
    """
}

dedup = {
    doc "Mark Duplicates with PicardTools"
    output.dir="$REFBASE/aligned_bams" 
    exec """
        java -Xmx6g -Djava.io.tmpdir=$TMPDIR -jar $PICARD_HOME/MarkDuplicates.jar
             INPUT=$input.bam
             REMOVE_DUPLICATES=true 
             VALIDATION_STRINGENCY=LENIENT 
             AS=true 
             METRICS_FILE=$LOG 
             OUTPUT=$output.bam
    """
}


@transform(".pass1.table")
bqsrPass1 = {
    doc "Recalibrate base qualities in a BAM file so that quality metrics match actual observed error rates"
    output.dir="aligned_bams"
    exec """
            java -Xmx12g -jar $GATK/GenomeAnalysisTK.jar
                -T BaseRecalibrator 
                -R $REF  
                -I $input.bam  
                -knownSites $REFSNP1  
                -knownSites $REFSNP2 
                -knownSites $REFSNP3
                -nt $NTHREAD_GATK
                -o $output.pass1.table 
        """
    
}

@transform(".pass2.table")
bqsrPass2 = {
    doc "Recalibrate base qualities in a BAM file so that quality metrics match actual observed error rates"
    output.dir="aligned_bams"
    exec """
            java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar 
                -T BaseRecalibrator 
                -R $REF 
                -I $input.bam 
                -knownSites $REFSNP1  
                -knownSites $REFSNP2 
                -knownSites $REFSNP3 
                -BQSR $input.pass1.table
                -nt $NTHREAD_GATK  
                -o $output.pass2.table 
        """
}

bqsrApply = {
    doc "Apply BQSR to input BAM file."
    output.dir="aligned_bams"
    exec """
        java -jar $GATK/GenomeAnalysisTK.jar  
            -T PrintReads  
            -R $REF 
            -I $input.bam   
            -BQSR $input.pass1.table  
            -o $output.bam
    """ 
}

bqsrCheck = {
    doc "Compare pre and post base-quality scores from recalibration"
    output.dir="qc"
    from(".pass1.table", ".pass2.table", ".bqsr.pdf") {
        exec """
            java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
              -T AnalyzeCovariates
              -R $REF
              -before $input.pass1.table 
              -after $input.pass2.table 
              -plots $ouptut.bqsr.pdf

        """
    }
}

fastqc = {
    doc "Run FASTQC to generate QC metrics for the reads post-alignment"
    output.dir = "qc"
    transform('.bam')  to('_fastqc.zip')  {
        exec "fastqc --quiet -o ${output.dir} -f bam_mapped $input.bam"
    }
}

// QC metrics at BAM level
flagstat = {
    output.dir = "qc"
    transform('.bam') to('.flagstat') {
        exec "samtools flagstat $input.bam > $output.flagstat"
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
              -o $output
        """ 
    }    
}

indexVCF = {
    exec "./vcftools_prepare.sh $input.vcf"
}


callVariants = {
    doc "Call SNPs/SNVs using GATK Haplotype Caller, produces .g.vcf"
    output.dir="variants"
    transform(".bam") to(".g.vcf") {}
    exec """
        java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar
            -T HaplotypeCaller 
            -R $REF
            -I $input.bam
            --emitRefConfidence GVCF
            -gt_mode DISCOVERY
            -o $output.g.vcf    
        """
}

//genotype = {
//    doc "Jointly genotype gVCF files"
//    output.dir="variants"
//    exec """

//        gvcf_files=$(find ./variants -name '*.g.vcf' -print0 | xargs -0 -I {} echo '--variant {}')
//        java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar
//            -T GenotypeGVCFs
//            $gvcf_files
//            -o all_raw.vcf
//    """
//}

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

// variant annotation using snpEff
annotateVariants = {
    doc "Annotate variants using snpEff"
    exec """ 
        java -Xmx8g -jar  $SNPEFF_HOME/snpEff.jar 
            -no-downstream 
            -no-upstream 
            -v  
            -c $configFile Pf3D7v3 all_recalibrated.vcf > all_recalibrated_anno.vcf
    """
}

addRegionsAnnotation = {
    doc "Include core regions from Pf genetic crosses"
    exec """
        bgzip $CORE_REGIONS
        tabix -s 1 -b 2 -e 3 $CORE_REGIONS.gz
        cat all_recalibrated_anno.vcf | vcf-annotate -a $CORE_REGIONS.gz \
            -d key=INFO,ID=RegionType,Number=1,Type=String,Description='The type of genome region within which the variant is found. 
            SubtelomericRepeat: repetitive regions at the ends of the chromosomes. 
            SubtelomericHypervariable: subtelomeric region of poor conservation between the 3D7 reference genome and other samples. 
            InternalHypervariable: chromosome-internal region of poor conservation between the 3D7 reference genome and other samples. 
            Centromere: start and end coordinates of the centromere genome annotation. Core: everything else.' \
            -c CHROM,FROM,TO,INFO/RegionType > all_recalibrated_anno_regions.vcf

    """
}

// variant filtering
@filter("filter")
filterSNPs = {
    output.dir="variants"
    exec """
        java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
            -T VariantFiltration 
            -R $REF 
            --filterExpression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || HaplotypeScore > 13.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' 
            --filterName 'GATK_MINIMAL_FILTER'
            --filterExpression 'VQSLOD <= 0 or RegionType != "Core"'
            --filterName 'MalariaGen Filter'
            --variant all_recalibrated_anno_regions.vcf
            -log $LOG 
            -o all_recalibrated_anno_regions_filtered.vcf
    """
}


