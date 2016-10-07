'''
Created on 6Oct.,2016

@author: thomas.e
'''

import drmaa
import os

HOME = os.getenv('HOME')
CMD = '/bin/hostname'
WAIT_TIME = 60

s=drmaa.Session()
s.initialize()

jt = s.createJobTemplate()
jt.workingDirectory = drmaa.JobTemplate.HOME_DIRECTORY
jt.remoteCommand = CMD
jt.outputPath = drmaa.JobTemplate.HOME_DIRECTORY
jt.nativeSpecification = '-q small'

ids = s.runBulkJobs(jt, 1, 1, 1)

print('Job started: ' + ids[0])

waiting = True
while waiting:
    try:
        s.synchronize(ids, WAIT_TIME, False)
        waiting = False
    except drmaa.errors.ExitTimeoutException:
        print('waiting...')
        
info=s.wait(ids[0], drmaa.Session.TIMEOUT_WAIT_FOREVER)
        
print( '-' * 76)
print("""
id:                        %(jobId)s
exited:                    %(hasExited)s
signaled:                  %(hasSignal)s
with signal (id signaled): %(terminatedSignal)s
dumped core:               %(hasCoreDump)s
aborted:                   %(wasAborted)s
resource usage:

%(resourceUsage)s
""" % info._asdict()
)