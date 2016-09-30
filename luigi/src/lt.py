'''
Created on 12Sep.,2016

@author: thomas.e
'''

from luigi.task import Task
from luigi.file import LocalTarget
import os
import socket
import pydevd

HOST = '10.1.17.38'
SOFTWARE = "/usr/local/bioinfsoftware"
HOME = os.getenv("HOME")

class InputFiles(Task):
    def output(self):
        return LocalTarget('data/ev.fastq')

class Trim(Task):
    def requires(self):
        return InputFiles()
    
    def output(self):
        return LocalTarget('data/ev-trimmed.fastq')
    
    def run(self):
        TRIMMOMATIC = "java -jar " + SOFTWARE + "/trimmomatic/current/trimmomatic-0.30.jar"
        ADAPTORS = "ILLUMINACLIP:" + SOFTWARE + "/trimmomatic/Trimmomatic-0.30/adapters/TruSeq3-SE.fa:2:30:10"
        
        inf = self.input().path
        outf = self.output().path
        
        cmd = TRIMMOMATIC + " SE " + inf + " " + outf + " " + ADAPTORS
        rc = os.system(cmd)
    
        if rc != 0:
            raise ChildProcessError("trimmomatic completed abnormally: rc=" + rc)



print('Running on host: ' + socket.gethostname())
print('Waiting for debug server...\n\n')
pydevd.settrace(HOST, stdoutToServer=True, stderrToServer=True)
