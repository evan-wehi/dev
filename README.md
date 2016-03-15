# Plasmodium falciparum WGS variant calling pipeline
This is my attempt to replicate the pipeline used by the [Malaria Genomics Consortium](https://www.malariagen.net/projects/parasite/pf)
for calling single nucleotide variants (SNVs) from paired-end whole-genome sequencing data.

It is a [bpipe](http://docs.bpipe.org/Overview/Introduction/) pipeline adapted from
https://github.com/ssadedin/variant_calling_pipeline/
https://github.com/claresloggett/variant_calling_pipeline/

## Overview

### What's included?


## Requirements and usage
There are a number of software requirements, which you should ensure are 
satisfied before running the pipeline. These need to be configured 
in the file called 'config.groovy'.

At the moment we use the following software in our pipeline
bwa mem
PicardTools
FastQC
GATK HaplotypeCaller
SnpEff
vcf-annotate

By default this pipeline will attempt to use all the available cores
on the computer it runs on. If you don't wish to do that, limit the 
concurrency by running it with the -n flag:

    bpipe run -n 4 pipeline.groovy isolate_data/*.bam
 
Assumes: files in form  *<sample_name>*.bam *<sample_name>*.bam
Assumes: Bpipe Version 0.9.8.5 or later



