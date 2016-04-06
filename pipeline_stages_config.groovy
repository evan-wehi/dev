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
samToFastq = {
    doc "Extract fastq from BAM file and soft-clip adapters, redirect output."
    def sm="$input".replaceAll(/^.*\/{1}/, '')
    def sample=sm.replaceAll('.bam', '')
    def outdir="$REFBASE/fastq/$sample"
    output.dir="logs"
    transform(".bam") to(".txt") {
        exec """
        echo 'Begin pipeline for isolate: $sample'
        
        mkdir -p $outdir;

        java -Xmx6g -jar $PICARD_HOME/RevertSam.jar
        VALIDATION_STRINGENCY=SILENT 
        INPUT=$input.bam
        OUTPUT=/dev/stdout
        MAX_RECORDS_IN_RAM=250000 
        SORT_ORDER=queryname
        TMP_DIR=$TMPDIR
        COMPRESSION_LEVEL=0 | java -Xmx6g -jar $PICARD_HOME/MarkIlluminaAdapters.jar 
        INPUT=/dev/stdin 
        OUTPUT=/dev/stdout
        PE=true
        ADAPTERS=PAIRED_END
        MAX_RECORDS_IN_RAM=250000 
        QUIET=true 
        COMPRESSION_LEVEL=0 
        M=$outdir/adapters.txt | java -Xmx6g -jar $PICARD_HOME/SamToFastq.jar
        INPUT=/dev/stdin
        MAX_RECORDS_IN_RAM=250000 
        CLIPPING_ATTRIBUTE=XT 
        CLIPPING_ACTION=2 
        OUTPUT_PER_RG=true 
        OUTPUT_DIR=$outdir
        VALIDATION_STRINGENCY=SILENT 
        TMP_DIR=$TMPDIR &> $output.txt
        """
    }
}


@filter("merged")
remapByRG = {
    doc "Map fastq files by RG using bwa mem and merge"
    def sm="$input.bam".replaceAll(/^.*\/{1}/, '')
    def sample=sm.replaceAll('.bam', '')
    def outdir="$REFBASE/fastq/$sample"
    output.dir="$REFBASE/aligned_bams"
    exec "./revert_remap.sh $input.bam $outdir $output.bam"  
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
            java -Xmx6g -jar $GATK/GenomeAnalysisTK.jar
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
            java -Xmx6g -jar $GATK/GenomeAnalysisTK.jar 
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
        java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar  
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
callVariants = {
    doc "Call SNPs/SNVs using GATK Haplotype Caller, produces .g.vcf"
    output.dir="$REFBASE/variants"
    transform(".bam") to(".g.vcf") {
        exec """
            java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar
                -T HaplotypeCaller 
                -R $REF
                -I $input.bam
                --emitRefConfidence GVCF
                -gt_mode DISCOVERY
                -variant_index_type LINEAR -variant_index_parameter 128000
                -o $output

        """
    }
}

combineGVCF = {
    doc "Jointly genotype gVCF files"
    output.dir="$REFBASE/variants_combined"
    def gvcfs = "--variant " + "$inputs.g.vcf".split(" ").join(" --variant ")
    produce("combined_variants.vcf") {
        exec """
            java -Xmx8g -jar $GATK/GenomeAnalysisTK.jar
                -T GenotypeGVCFs
                -R $REF
                $gvcfs
                -o $output
    """
    }
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
            -resource:cross1,training=true,truth=true,known=false $REFSNP1
            -resource:cross2,training=true,truth=true,known=false $REFSNP2
            -resource:cross3,training=true,truth=true,known=false $REFSNP3
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
            -c $SNPEFF_CONFIG
            -no-downstream 
            -no-upstream 
            Pf3D7v3
            $input.vcf > $output.vcf
    """
}

regions = {
    doc "Include core regions from Pf genetic crosses version 1"
    output.dir="$REFBASE/variants_combined"

    exec """
        cat $input.vcf | vcf-annotate -a $CORE_REGIONS \
            -d key=INFO,ID=RegionType,Number=1,Type=String,Description='The type of genome region within which the variant is found. 
            SubtelomericRepeat: repetitive regions at the ends of the chromosomes. 
            SubtelomericHypervariable: subtelomeric region of poor conservation between the 3D7 reference genome and other samples. 
            InternalHypervariable: chromosome-internal region of poor conservation between the 3D7 reference genome and other samples. 
            Centromere: start and end coordinates of the centromere genome annotation. Core: everything else.' \
            -c CHROM,FROM,TO,INFO/RegionType > $output.vcf

    """
}

barcode = {
    doc "Annotate global barcode SNPs from Neafsey et al., 2008"
    output.dir="$REFBASE/variants_combined"

    exec """
        cat $input.vcf | vcf-annotate -a $BARCODE \
            -d key=INFO,ID=GlobalBarcode,Number=1,Type=String,Description='Global Barcode SNP from Neafsey et al., 2008.' \
            -c CHROM,FROM,TO,INFO/GlobalBarcode > $output.vcf 
    """ 
}
// select only biallelic SNPs
keepSNPs= {
    doc "Kepp only SNPs in VCF file"
    output.dir="$REFBASE/variants_combined"
    produce("biallelic_only.vcf") {
        exec """
            java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
                -T SelectVariants \
                -R $REF
                -V $input.vcf \
                -o $output \
                -selectType SNP
                --restrictAllelesTo BIALLELIC
    """
    }
}

// variant filtering annotations
@filter("filter")
filterSNPs = {
    doc "Annotate VCF file with additional filters at the variant level"
    output.dir="$REFBASE/variants_combined"
    exec """
        java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar 
            -T VariantFiltration 
            -R $REF
            --clusterSize 3 
            --clusterWindowSize 30
            --filterName LowQualMG -filter "VQSLOD <= 0.0 && RegionType != 'Core'"
            --variant $input.vcf
            -log $LOG 
            -o $output.vcf
    """
}

// apply GATK SelectVariants to filter Low Qual regions 
cleanVCF = {
    doc "Clean VCF for analysis ready."
    output.dir="cache"
    produce("final_snps.vcf") {
        exec """
            java -Xmx4g -jar $GATK/GenomeAnalysisTK.jar
                -T SelectVariants
                -R $REF
                --variant $input.vcf
                -o $output
                -ef
        """
    }
}



indexVCF = {
    exec "./vcftools_prepare.sh $input.vcf"
}


cleanGDS = {
    doc "Use clean VCF to construct GDS file"
    output.dir="cache"
    def vcf="$input.vcf"+".gz"
    R {"""
        library("SeqArray");
        seqVCF2GDS("$vcf", "$output.gds")
    """
    }
}

extractAnno = {
    doc "Use snpSift to extract annotations as plain text file"
    output.dir="cache"
    produce("final_annotations.txt") {
        exec """
            cat $input.vcf | perl $SNPEFF_HOME/scripts/vcfEffOnePerLine.pl |\
            java -Xmx4g -jar $SNPEFF_HOME/SnpSift.jar extractFields -e "NA" - \
                CHROM POS REF ALT "ANN[*].ALLELE" GlobalBarcode DP "ANN[*].GENEID" "ANN[*].BIOTYPE" "ANN[*].EFFECT"  "ANN[*].HGVS_P" "ANN[*].ERRORS" > $output
        """

    }
}



