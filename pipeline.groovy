////////////////////////////////////////////////////////////
// GATK-based 'best' practices variant-calling pipeline for 
// Plasmodium falciparum.
// 
// This pipeline is a port of the VLSCI whole genome variant
// calling pipeline  for from this Git repo to Bpipe:
// 
//     https://github.com/claresloggett/variant_calling_pipeline/
// And modified from
//     https://github.com/ssadedin/variant_calling_pipeline/
//
//
// This pipeline attempts to recreate the WGS processing pipeline developed
// by the Malaria Genomics Consortium for Plasmodium falciparum. 
//
////////////////////////////////////////////////////////////

// Create this file by copying config.groovy.template and editing
load 'config.groovy'

// All the core pipeline stages in the pipeline
load 'pipeline_stages_config.groovy'


run {  "*.bam" * [ sampleID + realignIntervals + realignIndels + dedup + bqsrPass1 + bqsrPass2 + bqsrApply + bqsrCheck + fastqc + flagstat + depthOfCoverage + callVariants] }
//run { "*.bam" * [ remapBWAmem + realignIntervals + realignIndels + dedup ] }
//
//run {
    // Create index for reference file
//    [
//     createBWAIndex + fastaIndex
//    ] +
//    // Pre-process each unmapped BAM file separately in parallel
//    // Map, sort clean then dedup
//    // Realign indels
//    // Recalibrate base-quality scores
//    // Run QC
//   // Call variants
//    "*.bam" * [
//               remapBWAmem +
//               dedup +
//               realignIntervals +
//               realignIndels +
//               bqsrPass1 +
//              bqsrPass2 +
            //   bqsrCheck +
          //     bqsrApply +
        //       fastqc +
       //        flagstat +
     //          depthOfCoverage +
     //          callVariants

    //] + 
        // For each g.vcf will jointly call genotypes then apply VQSR
      //[ genotype + vqsrGenerate + vqsrApply + annotateVariants + filterSNPs

//}
