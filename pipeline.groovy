////////////////////////////////////////////////////////////
// GATK-based 'best' practices variant-calling pipeline.
// 
// This pipeline is a port of the VLSCI whole genome variant
// calling pipeline  for from this Git repo to Bpipe:
// 
//     https://github.com/claresloggett/variant_calling_pipeline/
// And modified from
//     https://github.com/ssadedin/variant_calling_pipeline/
//
// There are a number of software requirements, which you should ensure are 
// satisfied before running the pipeline. These need to be configured 
// in the file called 'config.groovy' (which is loaded below). A template
// is provided in config.groovy which you can use to 
// create the file and set the right paths to your tools and reference
// data.
//
// This pipeline attempts to recreate the WGS processing pipeline developed
// by the Malaria Genomics Consortium for Plasmodium falciparum. 
//
// By default this pipeline will attempt to use all the available cores
// on the computer it runs on. If you don't wish to do that, limit the 
// concurrency by running it with the -n flag:
//
//    bpipe run -n 4 pipeline.groovy example_data/input_data_wgs/*.fastq.gz
// 
// Assumes: paired end reads
// Assumes: files in form  *<sample_name>*_..._R1.fastq.gz, *<sample_name>*_..._R2.fastq.gz
//
//
////////////////////////////////////////////////////////////

// Create this file by copying config.groovy.template and editing
load 'config.groovy'

// All the core pipeline stages in the pipeline
load 'pipeline_stages_config.groovy'

run {
    // Create index for reference file
    [
     createBWAIndex + fastaIndex
    ] +
    // Pre-process each unmapped BAM file separately in parallel
    // Map, sort clean then dedup
    // Realign indels
    // Recalibrate base-quality scores
    // Run QC
    // Call variants
    "*.bam" * [
               remapBWAmem +
               dedup +
               realignIntervals +
               realignIndels +
               bqsrPass1 +
               bqsrPass2 +
               bqsrCheck +
               bqsrApply +
               fastqc +
               flagstat +
               depthOfCoverage +
               callVariants

    ] + 
        // For each g.vcf will jointly call genotypes then apply VQSR
      [ genotype + vqsrGenerate + vqsrApply + annotateVariants + filterSNPs

}
