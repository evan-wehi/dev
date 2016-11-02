'''
Created on 1Nov.,2016

@author: thomas.e
'''

from toil.job import Job
import argparse
import os
import glob

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

NP_THREADS = 10


def build_parser():
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--queue', default=None, action='store_true', help='Queue to run jobs on')
    return parser

def create_job(forward, backward):
    pass
    
def create_jobs():
    forwardFiles = glob.glob(FILE_PATTERN.format('1'))
    jobs = []
    for forward in forwardFiles:
        backward = forward[:-len('R1.fastq.gz')] + 'R2.fastq.gz'
        try:
            os.lstat(backward)
            exists = True
        except FileNotFoundError:
            exists = False
        if exists: jobs.append(create_job(forward, backward))

def main():
    parser = build_parser()
    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()

    jobs = create_jobs()
    
    
if __name__=="__main__":
    main()
