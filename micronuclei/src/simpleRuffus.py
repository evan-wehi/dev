'''
Created on 10Oct.,2016

@author: thomas.e
'''

from ruffus import Pipeline, suffix
from ruffus.drmaa_wrapper import run_job
from drmaa import Session
from glob import glob

input_file = glob('*.in')

def stage1(inf, outf, session):
    run_job(cmd_str = 'cp ' + inf + ' ' + outf,
            job_name = 'stage1',
            drmaa_session = session,
            job_other_options = '-q small')
    
    
def originate(stuff):
    pass

def pipeLineFactory(name, initial_files, session):
    
    pl = Pipeline(name)
    
    pl.originate(originate, initial_files)
    pl.transform(task_func = stage1,
                 input = originate,
                 filter = suffix('.in'),
                 output = '.out',
                 extras = [session])
    
    return pl
    
s = Session()
s.initialize()    

pl = pipeLineFactory('first pipe', input_file, s)
pl.run(multiprocess = 3)