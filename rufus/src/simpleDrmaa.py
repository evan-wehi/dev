'''
Created on 6Oct.,2016

@author: thomas.e
'''

import drmaa
import os

HOME = os.getenv('HOME')
CMD = 'exit 17'
WAIT_TIME = drmaa.Session.TIMEOUT_WAIT_FOREVER
JOB_QUEUE = 'thomas_e-1'

s=drmaa.Session()
s.initialize()

jt = s.createJobTemplate()
jt.workingDirectory = drmaa.JobTemplate.HOME_DIRECTORY
jt.remoteCommand = CMD
jt.outputPath = drmaa.JobTemplate.HOME_DIRECTORY
jt.jobCategory = JOB_QUEUE
jt.nativeSpecification = '-q ' + JOB_QUEUE

jid = s.runJob(jt)

print('Job started: ' + jid)

# waiting = True
# while waiting:
#     try:
#         info=s.wait(jid, 1)
#         waiting = False
#     except drmaa.errors.ExitTimeoutException:
#         print('waiting...')
        
        
info=s.wait(jid, -1)

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
print( '-' * 76)
print(info)