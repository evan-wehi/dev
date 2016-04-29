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

// this could probably be broken up into multiple steps i.e. extract, align, merge 
// but I'm still learning bpipe intricacies
// a bit hacky at the moment but for now it'll do.
@filter("reverted")
makeUBAM = {
    doc "Create uBAM from input BAM file using revert sam"
    branch.sample = branch.name
    output.dir="$REFBASE/ubam"
    uses(GB:8) {
        exec """
            $PICARD_HOME/RevertSam \
                    I=$input.bam \
                    O=$output.bam \
                    SANITIZE=true \
                    MAX_DISCARD_FRACTION=0.005 \
                    ATTRIBUTE_TO_CLEAR=XT \
                    ATTRIBUTE_TO_CLEAR=XN \
                    ATTRIBUTE_TO_CLEAR=AS \
                    ATTRIBUTE_TO_CLEAR=OC \
                    ATTRIBUTE_TO_CLEAR=OP \
                    SORT_ORDER=queryname \
                    RESTORE_ORIGINAL_QUALITIES=true \
                    REMOVE_DUPLICATE_INFORMATION=true \
                    REMOVE_ALIGNMENT_INFORMATION=true \
                    TMP_DIR=$TMPDIR 
        """
    }
}

fastqc_unmapped = {
    doc "Run FASTQC to generate QC metrics for the reads post-alignment"
    output.dir = "$REFBASE/qc"
    transform('.bam')  to('_fastqc.zip')  {
        exec "fastqc --quiet -o ${output.dir} -f bam $input.bam"
    }
}

@filter("markilluminaadapters")
markAdapters = {
    doc "Mark Illumina Adapters with PicardTools"
    output.dir="$REFBASE/ubam"
    def adapter_file="$output.dir" + "/" + "$sample" +"adapters.txt"
    uses(GB:8) {
        exec """ 
            $PICARD_HOME/MarkIlluminaAdapters \
                I=$input.bam \
                O=$output.bam \
                PE=true \
                M=$adapter_file \
                TMP_DIR=$TMPDIR

        """
    }
}

@preserve
@filter("mapped")
remapBWA = {
    doc "Create merged BAM alignment from unmapped BAM file"
    output.dir="$REFBASE/aligned_bams"
    indexes="$REFBASE/fasta/Pfalciparum.genome"
    def ubam="$REFBASE/ubam/"+"$sample"+".reverted.bam"
    uses(threads:8,GB:16) {
    exec """
        $PICARD_HOME/SamToFastq I=$input.bam FASTQ=/dev/stdout \
            CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
            TMP_DIR=$TMPDIR | \
        bwa mem -M -t 8 -p $indexes /dev/stdin | \
        $PICARD_HOME/MergeBamAlignment \
            ALIGNED_BAM=/dev/stdin \
            UNMAPPED_BAM=$ubam\
            OUTPUT=$output.bam \
            R=$REF CREATE_INDEX=true ADD_MATE_CIGAR=true \
            CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
            INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
            PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
            TMP_DIR=$TMPDIR


        """
    }

}
samToFastq = {
    doc "Extract fastq from BAM file and soft-clip adapters, redirect output."
    branch.sample = branch.name
    def outdir="$REFBASE/fastq/$sample"
    output.dir="logs"
    transform(".bam") to(".txt") {
        exec """
        echo 'Begin pipeline for isolate: $sample'
        
        mkdir -p $outdir;

        $PICARD_HOME/RevertSam \
            VALIDATION_STRINGENCY=SILENT \
            INPUT=$input.bam \
            OUTPUT=/dev/stdout \
            MAX_RECORDS_IN_RAM=250000 \
            SORT_ORDER=queryname \
            TMP_DIR=$TMPDIR \
            COMPRESSION_LEVEL=0 | $PICARD_HOME/MarkIlluminaAdapters \
            INPUT=/dev/stdin \
            OUTPUT=/dev/stdout \
            PE=true \
            QUIET=true \
            M=$outdir/adapters.txt \
            COMPRESSION_LEVEL=0 | $PICARD_HOME/SamToFastq \
            INPUT=/dev/stdin \
            CLIPPING_ATTRIBUTE=XT \
            CLIPPING_ACTION=2 \
            NON_PF=true \
            OUTPUT_PER_RG=true \
            RG_TAG=ID \
            OUTPUT_DIR=$outdir \
            VALIDATION_STRINGENCY=SILENT 
            TMP_DIR=$TMPDIR &> $output.txt
        """
    }
}


@filter("merged")
remapByRG = {
    doc "Map fastq files by RG using bwa mem and merge"
    def fastqdir="$REFBASE/fastq/$sample"
    output.dir="$REFBASE/aligned_bams"
    exec "./revert_remap.sh $input.bam $fastqdir $output.bam"  
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

@filter("realignindels")
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
@filter("dedup")
dedup = {
    doc "Mark Duplicates with PicardTools"
    output.dir="$REFBASE/aligned_bams"
    def metrics="$REFBASE/qc/" + "$sample" +".dup.metrics.txt" 
    exec """
        $PICARD_HOME/MarkDuplicates
             INPUT=$input.bam
             REMOVE_DUPLICATES=true 
             VALIDATION_STRINGENCY=LENIENT 
             AS=true 
             METRICS_FILE=$metrics
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

@preserve
@filter("final")
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


// QC metrics at BAM level
fastqc_mapped = {
    doc "Run FASTQC to generate QC metrics for the reads post-alignment"
    output.dir = "$REFBASE/qc"
    transform('.bam')  to('_fastqc.zip')  {
        exec "fastqc --quiet -o ${output.dir} -f bam_mapped $input.bam"
    }
}

@transform(".alignment_metrics")
alignment_metrics = {
    doc "Collect alignment summary statistics"
    output.dir="$REFBASE/qc"

    exec """ 
        $PICARD_HOME/CollectAlignmentSummaryMetrics \
            R=$REF \
            I=$input.bam \
            O=$output.alignment_summary
    """
}

@transform(".insert_metrics")
insert_metrics = {
    doc "Collect insert size metrics"
    output.dir="$REFBASE/qc"
    def histogram="$output.dir" + "/" + "$sample" + "is_distribution.pdf"
    exec """
        $PICARD_HOME/CollectInsertSizeMetrics \
            I=$input.bam \
            O=$output.insert_metrics \
            H=$histogram 
    """

}

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
                -o $output

        """
    }
}

@preserve
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
        mv snpEff_summary.html snpEff_genes.txt logs/
    """
}

regions = {
    doc "Include core regions from Pf genetic crosses version 1"
    output.dir="$REFBASE/variants_combined"

    exec """

        bcftools annotate -a $CORE_REGIONS -h $CORE_REGIONS_HDR -Ov -o $output.vcf -c CHROM,FROM,TO,RegionType $input.vcf


    """
}

barcode = {
    doc "Annotate global barcode SNPs from Neafsey et al., 2008"
    output.dir="$REFBASE/variants_combined"

    exec """
        bcftools annotate -a $BARCODE -h $BARCODE_HDR -Ov -o $output.vcf -c CHROM,FROM,TO,GlobalBarcode $input.vcf

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
            --filterName LowQualVQ -filter "VQSLOD <= 0.0"
            --filterName NotCore -filter  "RegionType != 'Core'"
            --variant $input.vcf
            -log $LOG 
            -o $output.vcf
    """
}

// apply GATK SelectVariants to filter Low Qual regions
@preserve 
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
                -select 'vc.isNotFiltered()'
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



