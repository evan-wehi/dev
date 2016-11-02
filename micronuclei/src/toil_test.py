'''
Created on 2Nov.,2016

@author: thomas.e
'''

from toil.job import Job
import os.path
import os

STEPS = 2

def fn(job, vars):
    step, i, fileHandle = vars
    
    if step == 0:
        fn = fileHandle
    else:
        fn = job.fileStore.readGlobalFile(fileHandle)
        
    print('Processing: ' + fn + ' for input ' + i)
    
    last = step == STEPS
    
    if last:
        outfn = str(i) + '.out'
    else:  
        outfn = job.fileStore.getLocalTempFile()
        
    cmd = 'cp ' + fn + ' ' + outfn
    print('Executing cmd=' + cmd)
    os.system('cp ' + fn + ' ' + outfn)
    
    if last:
        outHandle = outfn
    else:
        outHandle = job.fileStore.writeGlobalFile(outfn)
        job.addChildFn(fn, (step+1, i, outHandle))
        
    
def main():
    
    files = ['1.in', '2.in']
    mj = Job()
    for i in range(len(files)):
        file = files[i]
        mj.addChild(Job.wrapJobFn(fn, (0, i, file)))
        
    options = Job.Runner.getDefaultOptions("./toilWorkflowRun")
    options.logLevel = "INFO"
    
    Job.Runner.startToil(mj, options)
        
        
if __name__ == '__main__':
    main()