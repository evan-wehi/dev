// Set location of you reference files here (see below for the files required)
REFBASE="./data"

// Set a good location for storing large temp files here (probably not /tmp)
TMPDIR="/usr/local/work/lee.s"

// Set location of Picard tools here
PICARD_HOME="/usr/local/bioinfsoftware/picard-tools/current/jars"

// Set to the reference FASTA file, which must be indexed
REF="$REFBASE/fasta/Pfalciparum.genome.fasta"

// Set to a VCF file containing genetic cross entries
REFSNP1="data/known_sites/3d7_hb3.combined.final.vcf.gz"      
REFSNP2="data/known_sites/7g8_gb4.combined.final.vcf.gz"
REFSNP3="data/known_sites/hb3_dd2.combined.final.vcf.gz"

// Log data from various tools will appear in here
LOG="pipeline.log"

// Set GATK location here
GATK="/usr/local/bioinfsoftware/gatk/current/"
GATK35="~/GenomeAnalysisTK-3.5"
// Set location of snpEff here
SNPEFF_HOME="/usr/local/bioinfsoftware/snpEff/current"
//Set location of snpEff config file
SNPEFF_CONFIG="$REFBASE/annotations/snpEff.config"

// set location for region based annotations
CORE_REGIONS="$REFBASE/annotations/regions-20130225.onebased.txt"

// set resistance genes interval_list
RESISTANCE_LOCI="$REFBASE/annotations/drug_resistance_genes.interval_list"
