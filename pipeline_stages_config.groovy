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
    output.dir="$REFBASE/fasta"
    from("fasta") to("dict") {
        exec "java -jar $PICARD_HOME/CreateSequenceDictionary.jar R=$REF O=$output.dict"
    }
}

// only works if PicardTools v2
//updateVCF = {
//    doc "Update VCF with new sequence dictionary"
//    exec """
//        java -jar $PICARD_HOME/SortVcf.jar
//          I=$REFSNP1
//          SEQUENCE_DICTIONARY=$REFBASE/fasta/Pfalciparum.genome.dict 
//          O=$REFSNP1
//    """
//}

// this could probably be broken up into multiple steps i.e. extract, align, merge 
// but I'm still learning bpipe intricacies
// a bit hacky at the moment but for now it'll do.
sampleID = {
    // extract sample name and save as  branch variable
    sm="$input".replaceAll(/^.*\/{1}/, '')
    branch.sample=sm.replaceAll('.bam', '')
    exec """
        echo 'Begin pipeline for isolate: $sample'
    """
}

samToFastq = {
    doc "Extract fastq from BAM file and soft-clip adapters, redirect output."
    output.dir="$REFBASE/fastq/$sample"
    transform(".bam") to(".txt") {
        exec """
        
        mkdir -p $output.dir;

        java -Xmx6g -jar $PICARD_HOME/RevertSam.jar
        VALIDATION_STRINGENCY=SILENT 
        INPUT=$input.bam
        OUTPUT=/dev/stdout 
        SORT_ORDER=queryname
        TMP_DIR=$TMPDIR
        COMPRESSION_LEVEL=0 | java -Xmx6g -jar $PICARD_HOME/MarkIlluminaAdapters.jar 
        INPUT=/dev/stdin 
        OUTPUT=/dev/stdout
        PE=true
        ADAPTERS=PAIRED_END 
        COMPRESSION_LEVEL=0 
        M=$REFBASE/fastq/$sample/adapters.txt | java -Xmx6g -jar $PICARD_HOME/SamToFastq.jar
        INPUT=/dev/stdin 
        CLIPPING_ATTRIBUTE=XT 
        CLIPPING_ACTION=2 
        OUTPUT_PER_RG=true 
        OUTPUT_DIR=$output.dir
        VALIDATION_STRINGENCY=SILENT 
        TMP_DIR=$TMPDIR &> $output.txt
        """
    }
}


@filter("merged")
remapByRG = {
    doc "Map fastq files by RG using bwa mem and merge"
    output.dir="$REFBASE/aligned_bams"
    exec "./revert_remap.sh $input.bam $REFBASE/fastq/$sample $output.bam"  
}


@transform(".intervals")
realignIntervals = {
    output.dir="$REFBASE/aligned_bams"
    exec """
        java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
          -T RealignerTargetCreator
          -R $REF 
          -I $input.bam
          -log $LOG
          -o $output.intervals
    """
}

realignIndels = {
    output.dir="$REFBASE/aligned_bams"
    exec """
        java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar 
          -T IndelRealigner 
          -R $REF 
          -I $input.bam
          -targetIntervals $input.intervals 
          -log $LOG
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
             CREATE_INDEX=true 
             OUTPUT=$output.bam
    """
}

bqsrPass1 = {
    doc "Recalibrate base qualities in a BAM file so that quality metrics match actual observed error rates"
    output.dir="$REFBASE/aligned_bams"
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
    output.dir="$REFBASE/aligned_bams"
    exec """
            java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar 
                -T BaseRecalibrator 
                -R $REF 
                -I $input.bam 
                -knownSites $REFSNP1  
                -knownSites $REFSNP2 
                -knownSites $REFSNP3 
                -BQSR $input.pass1.table
                -o $output.pass2.table 
        """
}

bqsrApply = {
    doc "Apply BQSR to input BAM file."
    output.dir="$REFBASE/aligned_bams"
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
    output.dir="$REFBASE/qc"
    transform('.pass1.table', '.pass2.table') to('.csv') {
        exec """
            java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
              -T AnalyzeCovariates
              -R $REF
              -before $input.pass1.table 
              -after $input.pass2.table 
              -csv $output
        """
    }
}

fastqc = {
    doc "Run FASTQC to generate QC metrics for the reads post-alignment"
    output.dir = "$REFBASE/qc"
    transform('.bam')  to('_fastqc.zip')  {
        exec "fastqc --quiet -o ${output.dir} -f bam_mapped $input.bam"
    }
}

// QC metrics at BAM level
flagstat = {
    output.dir = "$REFBASE/qc"
    transform('.bam') to('.flagstat') {
        exec "samtools flagstat $input.bam > $output.flagstat"
    }
}



depthOfCoverage = {
    output.dir="$REFBASE/qc"
    transform("bam") to("sample_statistics", "sample_summary") {
        exec """
            java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
              -T DepthOfCoverage 
              -R $REF -I $input.bam 
              -omitBaseOutput 
              -ct 1 -ct 5 -ct 10
              -mmq 20 -mbq 20
              -L $RESISTANCE_LOCI
              -o $output.prefix
        """ 
    }    
}



indexVCF = {
    exec "./vcftools_prepare.sh $input.vcf"
}


callVariants = {
    doc "Call SNPs/SNVs using GATK Haplotype Caller, produces .g.vcf"
    output.dir="$REFBASE/variants"
    transform("bam") {
        exec """
            java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar
                -T HaplotypeCaller 
                -R $REF
                -I $input.bam
                --emitRefConfidence GVCF
                -gt_mode DISCOVERY
                -o $sample.g.vcf
        """
    }
}

combineGVCF = {
    doc "Jointly genotype gVCF files"
    output.dir="$REFBASE/variants_combined"
    gvcfs="$inputs.g.vcf"
    exec """
        echo $gvcfs
    """
}

// variant recalibration
vqsrGenerate = {
    doc "Generate variant quality score recalibration. Note requires GATK version 3.5"
    output.dir="$REFBASE/variants_combined"
    exec """
         java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar 
            -T VariantRecalibrator 
            -R $REF 
            -input $input.vcf
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
    output.dir="$REFBASE/variants_combined"
    exec """
        java -jar $GATK/GenomeAnalysisTK.jar 
            -T ApplyRecalibration 
            -R $REF  
            -input $input.vcf  
            -mode SNP  
            --ts_filter_level 99.0  
            -recalFile $input.recal  
            -tranchesFile $input.tranches  
            -o $output.vcf
    """
}

// variant annotation using snpEff
annotate = {
    doc "Annotate variants using snpEff"
    output.dir="$REFBASE/variants_combined"
    exec """ 
        java -Xmx8g -jar  $SNPEFF_HOME/snpEff.jar 
            -no-downstream 
            -no-upstream 
            -v 
            -c $configFile Pf3D7v3 $input.vcf > $output.vcf
    """
}

regions = {
    doc "Include core regions from Pf genetic crosses"
    exec """
        bgzip $CORE_REGIONS
        tabix -s 1 -b 2 -e 3 $CORE_REGIONS.gz
        cat $input.vcf | vcf-annotate -a $CORE_REGIONS.gz \
            -d key=INFO,ID=RegionType,Number=1,Type=String,Description='The type of genome region within which the variant is found. 
            SubtelomericRepeat: repetitive regions at the ends of the chromosomes. 
            SubtelomericHypervariable: subtelomeric region of poor conservation between the 3D7 reference genome and other samples. 
            InternalHypervariable: chromosome-internal region of poor conservation between the 3D7 reference genome and other samples. 
            Centromere: start and end coordinates of the centromere genome annotation. Core: everything else.' \
            -c CHROM,FROM,TO,INFO/RegionType > $output.vcf

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
            --filterExpression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' 
            --filterName 'GATK_MINIMAL_FILTER'
            --filterExpression 'VQSLOD <= 0 or RegionType != "Core"'
            --filterName 'MALARIAGEN_FILTER'
            --variant $input.vcf
            -log $LOG 
            -o $output.vcf
    """
}


