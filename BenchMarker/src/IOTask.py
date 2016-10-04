"""
    Perform IO in a thread using dd
"""

from multiprocessing import Value
from multiprocessing import Process
from subprocess import Popen
from timeit import default_timer
from tempfile import mkstemp
from os import remove

class IOTask(Process):
    
    def __init__(self, size=10**9, blksize=10*7, dir=None, read=True, write=True):
        Process.__init__(self, target=self.go)
        self.size = size
        self.blksize = blksize
        self.read = read
        self.write = write
        
        self.fn = mkstemp(dir=dir)[1]
        
        self.readTime = Value('d', 0.0)
        self.writeTime = Value('d', 0.0)
        
    def go(self):
        try:
            if self.write:
                cmd = 'dd conv=fdatasync if=/dev/zero of=' + self.fn + ' count=' + str(int(round(self.size/self.blksize))) + ' bs=' + str(self.blksize)
                self.writeTime.value = self.execute(cmd)
            
            if self.read:
                cmd = 'dd of=/dev/null if=' + self.fn + ' count=' + str(int(round(self.size/self.blksize))) + ' bs=' + str(self.blksize)
                self.readTime.value = self.execute(cmd)
        finally:
            self.cleanup()

    def execute(self, cmd):
        t = default_timer()
        sp = Popen(cmd.split(), stdout=open(r'/dev/null'), stderr=open(r'/dev/null'))
        sp.wait()
        return default_timer() - t
            
    def cleanup(self):
        if self.fn is not None:
            remove(self.fn)
            self.fn = None
    
    def readRate(self):
        return self.rate(self.readTime.value)
    
    def writeRate(self):
        return self.rate(self.writeTime.value)
    
    def rate(self, dt):
        if dt == 0:
            return 0.0
        else:
            return float(self.size) / float(dt) / 1000.0