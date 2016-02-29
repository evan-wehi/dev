# Plasmodium falciparum WGS variant calling pipeline
This is my attempt to replicate the pipeline used by the [Malaria Genomics Consortium](https://www.malariagen.net/projects/parasite/pf)
for calling single nucleotide variants (SNVs) from paired-end whole-genome sequencing data.

It is a bpipe pipeline adapted from https://github.com/ssadedin/variant_calling_pipeline/ 

## Overview

This is a port to Bpipe of the basic variant-calling and annotation pipeline developed at the 
Victorian Life Sciences Computation Initiative (VLSCI), University of Melbourne.

See https://github.com/claresloggett/variant_calling_pipeline/

This example is intended as an illustration of how a full pipeline
looks in Bpipe, especially for those who might be familiar with Ruffus.
There are a number of places where the pipeline is non-ideal, but it has
been kept that way to reflect the original as accurately as possible.

One particular point of note is that the alignment output is stored
in SAM format by some of the intermediate stages which would be very
unadvisable for true whole-genome data. It is easy to make it store
output in BAM format instead, so this would be a very advisable 
change if you intend to run this pipeline on a lot of data. 

There are a number of software requirements, which you should ensure are 
satisfied before running the pipeline. These need to be configured 
in the file called 'config.groovy'. A template
is provided in config.groovy.template which you can use to 
create the file and set the right paths to your tools and reference
data.

By default this pipeline will attempt to use all the available cores
on the computer it runs on. If you don't wish to do that, limit the 
concurrency by running it with the -n flag:

    bpipe run -n 4 pipeline.groovy example_data/input_data_wgs/*.fastq.gz
 
Assumes: paired end reads
Assumes: files in form  *<sample_name>*_..._R1.fastq.gz, *<sample_name>*_..._R2.fastq.gz
Assumes: Bpipe Version 0.9.8.5 or later


Note: if you just want to have a look at how the pipeline appears, you can find all the meat 
in the following two files:

  https://github.com/ssadedin/variant_calling_pipeline/blob/master/pipeline_stages_config.groovy
  https://github.com/ssadedin/variant_calling_pipeline/blob/master/pipeline.groovy

