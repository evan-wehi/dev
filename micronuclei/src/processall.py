'''
Created on 16Sep.,2016

@author: thomas.e
'''

import glob
import os
import sys
import subprocess
import traceback

HOME = '/wehisan/home/allstaff/t/thomas.e'
QUEUE = 'thomas_e-1'

#SOFTWARE = '/wehisan/general/system/bioinf-software/bioinfsoftware'
SOFTWARE = '/home/thomas.e/software'
REF_SEQ ='/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/reference_genomes/human_new/hg38.fa'

BASE_DIR = '/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/micronuclei'
DATA_DIR = BASE_DIR + '/data'
WORK_DIR = BASE_DIR + '/workdir'
SCRIPTS_DIR = BASE_DIR + '/scripts'

FILE_PATTERN = DATA_DIR + '/' + r'*_R[{}].fastq.gz'

NP_THREADS = 20

def diag(msg):
    sys.__stdout__.write(msg + '\n')
    sys.__stdout__.flush()
    
        
# HOST = '10.1.17.51'
# import pydevd
# pydevd.settrace(HOST, stdoutToServer=True, stderrToServer=True)

# import logger
# LOGGER = logger.Logger(HOME + '/p.log')
# sys.stderr = LOGGER.redirecter('main stderr')
# sys.stdout = LOGGER.redirecter('main stdout')

def processFile(forward, backward):
    makedir(SCRIPTS_DIR)
    makedir(WORK_DIR)

    baseName = getBaseName(forward);
    
    script = initScript(baseName)
    executor = lambda cmd, outfn=None, infn=None: scriptWriter(cmd, script, outfn=outfn, infn=infn)
#     executor = lambda cmd, outfn=None: osExecutor(cmd, LOGGER, outfn)
    
    """
    Trim
    """ 
    trimmedName = WORK_DIR +'/' + baseName + '-trimmed.fastq.gz'
    trim(forward, backward, trimmedName, executor)
    print('echo trim done\n', file=script)
    
    """
    Align
    """
    forward  = WORK_DIR + '/' + baseName + '-trimmed_1P.fastq.gz'
    backward = WORK_DIR + '/' + baseName + '-trimmed_2P.fastq.gz'
    alignedName = WORK_DIR + '/' + baseName + '-aligned.bam'
    align(forward, backward, alignedName, executor)
    print('echo align done\n', file=script)

    """
    Sort and index
    """
    sortedName = WORK_DIR + '/' + baseName + '-sorted.bam'
    samSort(alignedName, sortedName, executor)
    print('echo sort done\n', file=script)
    
    """
    call with GRIDSS
    """
    gridss(sortedName, baseName, executor)
    print('echo call done\n', file=script)
    
def samSort(alignedName, sortedName, executor):
    SAM = SOFTWARE + '/samtools/samtools-1.3.1/bin/samtools'
    sort = [SAM + ' view -b -@ ' + str(NP_THREADS) + ' ' + alignedName, SAM + ' sort -@ ' + str(NP_THREADS) + ' - > ' + sortedName]
    executor(sort, infn=alignedName, outfn=sortedName)


def initScript(baseName):
    fn = SCRIPTS_DIR + '/' + baseName + '.sh'
    script = open(fn, 'w')
    print('#!/bin/sh\n', file=script)
    print('#PBS -l nodes=1:ppn=' + str(NP_THREADS) + '\n', file=script)
    return script
    
def scriptWriter(cmd, script, outfn=None, infn=None):
    outfn = '' if outfn is None else ' > ' + outfn
    infn  = '' if infn  is None else ' < ' + infn
    
    if type(cmd) is not list:
        cmd = [cmd]
        
    fc = cmd[0] + infn
    for c in cmd[1:]:
        fc = fc + ' | ' + c
    fc = fc + outfn
        
    print(fc + '\n', file=script)
    
def osExecutor(cmd, logger, outfn=None):
    if outfn is None:
        sp = subprocess.Popen(cmd.split(), 
                          bufsize=1,
                          stdout=subprocess.PIPE, 
                          stderr=subprocess.PIPE
                          )
    
        logger.redirect('child stdout', sp.stdout)
    else:
        outfn = open(outfn, 'wb')
        sp = subprocess.Popen(cmd.split(), 
                          stdout=outfn, 
                          stderr=subprocess.PIPE
                          )
    

    logger.redirect('child stderr', sp.stderr)
    sp.wait()
    rc = sp.returncode

    if rc:
        raise ChildProcessError(cmd + '\ncompleted abnormally: rc=' + str(rc))
    
def gridss(infile, baseName, executor):
    GRIDSS_JAR = HOME + '/dev//micronuclei/gridss-0.11.7-jar-with-dependencies.jar'
    GRIDSS = 'java -ea -Xmx16g -cp ' + GRIDSS_JAR
    
    outfile = WORK_DIR + '/' + baseName + '.vcf'
    cmd = GRIDSS + ' au.edu.wehi.idsv.Idsv' + \
        ' TMP_DIR=/tmp WORKING_DIR=/tmp REFERENCE=' + REF_SEQ + \
        ' INPUT=' + infile + ' IC=1' + \
        ' OUTPUT=' + outfile + \
        ' THREADS=' + str(NP_THREADS)
    
    executor(cmd)
    
    bedFile = WORK_DIR + '/' + baseName + '.bedpe'
    filteredBedFile = WORK_DIR + '/' + baseName + '.filtered.bedpe'
    cmd = GRIDSS + ' au.edu.wehi.idsv.VcfBreakendToBedpe' + \
    ' INPUT=' + outfile + \
    ' OUTPUT=' + bedFile + \
    ' OUTPUT_FILTERED=' + filteredBedFile + \
    ' REFERENCE=' + REF_SEQ

    executor(cmd)
    
def align(forward, backward, out, executor):
    BWA = SOFTWARE + '/bwa/current/bin/bwa'
    cmd = BWA + ' mem -t ' + str(NP_THREADS) + ' ' + REF_SEQ + ' ' + forward + ' ' + backward 
    executor(cmd, outfn=out)
    
def getBaseName(inf):
    baseName = os.path.basename(inf)
    baseName, ext = os.path.splitext(baseName)
    while ext is not '':
        baseName, ext = os.path.splitext(baseName)
        
    return baseName[:-3]
    
def makedir(wdir):
    try:
        os.mkdir(wdir)
    except FileExistsError:
        pass
    
def trim(forward, backward, out, executor):
    TRIMMOMATIC = 'java -jar ' + SOFTWARE + '/trimmomatic/current/trimmomatic-0.36.jar'
    ADAPTORS = 'ILLUMINACLIP:' + SOFTWARE + '/trimmomatic/current/adapters/TruSeq3-PE.fa:1:30:20:4:true'

    cmd = TRIMMOMATIC + ' PE ' + '-threads ' + str(NP_THREADS) + ' ' + forward +' ' + backward + ' -baseout ' + out + ' ' + ADAPTORS
    
    executor(cmd)

if __name__ == '__main__':
    forwardFiles = glob.glob(FILE_PATTERN.format('1'))
    for forward in forwardFiles:
        backward = forward[:-len('R1.fastq.gz')] + 'R2.fastq.gz'
        try:
            os.lstat(backward)
            exists = True
        except FileNotFoundError:
            exists = False
        try:
            if exists: processFile(forward, backward)
        except Exception as ex:
            traceback.print_exc()

