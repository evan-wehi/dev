'''
A simple ruffus test 
'''

import os
from socket import gethostname
import subprocess
import sys

import pydevd
from ruffus import suffix, pipeline_run, transform
from ruffus.task import follows

HOST = '10.1.17.131'
SOFTWARE = "/usr/local/bioinfsoftware"
HOME = os.getenv("HOME")

pydevd.settrace(HOST, stdoutToServer=True, stderrToServer=True)
print('Running on host: ' + gethostname())

initial_files = ["data/ev.fastq"]

@transform(initial_files, suffix(".fastq"), "-trimmed.fastq")
def trim(inf, outf):
    TRIMMOMATIC = "java -jar " + SOFTWARE + "/trimmomatic/current/trimmomatic-0.30.jar"
    ADAPTORS = "ILLUMINACLIP:" + SOFTWARE + "/trimmomatic/Trimmomatic-0.30/adapters/TruSeq3-SE.fa:2:30:10"
    
    cmd = TRIMMOMATIC + " SE " + inf + " " + outf + " " + ADAPTORS
    rc = os.system(cmd)
    
    if rc != 0:
        raise ChildProcessError("trimmomatic completed abnormally: rc=" + rc)


@transform(trim, suffix("-trimmed.fastq"), "-aligned.sam")
def align(inf, outf):
    BWA=SOFTWARE + "/bwa/current/bin/bwa"
    REF_SEQ = HOME + "/local/share/bcbio/genomes/Hsapiens/GRCh37/bwa/GRCh37.fa"
    
    """
    Index the FQ file
    """
    rc = os.system(BWA + " index " + inf)
    if rc != 0:
        raise ChildProcessError("bwa index completed abnormally: rc=" + str(rc))
    
    """
    Align to a SAM file
    """
    cmd = [BWA, 'mem', '-R', r'@RG\tID:Seq01p\tSM:Seq01\tPL:ILLUMINA\tPI:330', REF_SEQ, inf]
    outfh = open(outf, 'wb')
    p = subprocess.Popen(cmd, stdout=outfh, stderr=sys.stderr)
    rc = p.wait()
    if rc != 0:
        raise ChildProcessError("bwa mem completed abnormally: rc=" + rc)

    
@transform(align, suffix("-aligned.sam"), "-aligned-sorted.bam")
def convert_to_sorted_bam(inf, outf):
    SAM = SOFTWARE + "/samtools/samtools-1.3.1/bin/samtools"
    
    convertCmd = [SAM, "view",  "-b", inf]
    sortCmd = [SAM, 'sort', '-']
    outfh = open(outf, 'wb')
    pipein, pipeout = os.pipe()
    
    try:
        convertProcess = subprocess.Popen(convertCmd, stderr=sys.stderr, stdout=pipeout)
        sortProcess = subprocess.Popen(sortCmd, stdin=pipein, stdout=outfh, stderr=sys.stderr)
        
        rc = convertProcess.wait()
        if rc != 0:
            raise ChildProcessError("samtool view -b " + inf + ": rc=" + rc)
        
        rc = sortProcess.wait()
        if rc != 0:
            raise ChildProcessError("samtool sort " + inf + ": rc=" + rc)
            
    finally:
        outfh.close()
    
@transform(convert_to_sorted_bam, suffix('-aligned-sorted.bam'), '-aligned-sorted.bai')
def convert_to_indexed_bam(inf, outf):
    SAM = SOFTWARE + '/samtools/samtools-1.3.1/bin/samtools'
    
    cmd = SAM + ' index ' + inf + ' ' + outf
    rc = os.system(cmd)
    if rc != 0:
        raise ChildProcessError(cmd + ' completed abnormally: rc=' + str(rc)) 
    
    
@follows(convert_to_indexed_bam)
@transform(convert_to_sorted_bam, suffix('-aligned-sorted.bam'), '.vcf')
def call(inf, outf):
    REF_SEQ = HOME + '/local/share/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa'
    MUTECT2 = 'java -jar ' + SOFTWARE + '/gatk/current/GenomeAnalysisTK.jar -T MuTect2'
    
    cmd = MUTECT2 + ' --artifact_detection_mode -R ' + REF_SEQ + ' -I:tumor ' +  inf + ' -o ' + outf
    
    rc = os.system(cmd)
    if rc != 0:
        raise ChildProcessError('mutect2 completed abnormally: rc=' + str(rc))
        
pipeline_run()
