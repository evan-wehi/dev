'''
Created on 9Sep.,2016

@author: thomas.e
'''

import sys
import os
import ruffus
import pydevd
import socket
from ruffus.task import Pipeline
import subprocess

HOST = '10.1.17.131'
SOFTWARE = "/usr/local/bioinfsoftware"
HOME = os.getenv("HOME")

def originate(stuff):
    pass

def trim(inf, outf):
    TRIMMOMATIC = "java -jar " + SOFTWARE + "/trimmomatic/current/trimmomatic-0.30.jar"
    ADAPTORS = "ILLUMINACLIP:" + SOFTWARE + "/trimmomatic/Trimmomatic-0.30/adapters/TruSeq3-SE.fa:2:30:10"
    
    cmd = TRIMMOMATIC + " SE " + inf + " " + outf + " " + ADAPTORS
    rc = os.system(cmd)
    
    if rc != 0:
        raise ChildProcessError("trimmomatic completed abnormally: rc=" + rc)


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
    
def convert_to_indexed_bam(inf, outf):
    SAM = SOFTWARE + '/samtools/samtools-1.3.1/bin/samtools'
    
    cmd = SAM + ' index ' + inf + ' ' + outf
    rc = os.system(cmd)
    if rc != 0:
        raise ChildProcessError(cmd + ' completed abnormally: rc=' + str(rc)) 
    
    
def call(inf, outf):
    REF_SEQ = HOME + '/local/share/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa'
    MUTECT2 = 'java -jar ' + SOFTWARE + '/gatk/current/GenomeAnalysisTK.jar -T MuTect2'
    
    cmd = MUTECT2 + ' --artifact_detection_mode -R ' + REF_SEQ + ' -I:tumor ' +  inf + ' -o ' + outf
    
    rc = os.system(cmd)
    if rc != 0:
        raise ChildProcessError('mutect2 completed abnormally: rc=' + str(rc))


def pipeLineFactory (name, initial_files):
    
    pl = Pipeline(name)
    
    pl.originate(originate, initial_files)
    
    pl.transform(task_func = trim,
                 name = 'trim',
                 input = originate,
                 filter = ruffus.suffix(".fastq"),
                 output = '-trimmed.fastq'
                 )
    
    pl.transform(task_func = align,
                 name = 'align',
                 input = ruffus.output_from('trim'),
                 filter = ruffus.suffix('-trimmed.fastq'),
                 output = '-aligned.sam'
                 )
    
    pl.transform(task_func = convert_to_sorted_bam,
                 name = 'convert to sorted sam',
                 input = ruffus.output_from('align'),
                 filter = ruffus.suffix('-aligned.sam'),
                 output = '-aligned-sorted.bam'
                 )
    
    pl.transform(task_func = convert_to_indexed_bam,
                 name = 'convert to indexed bam',
                 input = ruffus.output_from('convert to sorted sam'),
                 filter = ruffus.suffix('-aligned-sorted.bam'),
                 output = '-aligned-sorted.bai'
                 )
    
    pl.transform(task_func = call,
                 name = 'call',
                 input = ruffus.output_from('convert to sorted sam'),
                 filter = ruffus.suffix('-aligned-sorted.bam'),
                 output = '.vcf'
                 ).follows('convert to indexed bam')
    
    return pl

print('Waiting for debug server...\n\n')
pydevd.settrace(HOST, stdoutToServer=True, stderrToServer=True)
print('Running on host: ' + socket.gethostname())

initial_files = ["data/ev.fastq"]

pl = pipeLineFactory('first pipe', initial_files)
pl.run()