'''
Created on 27Sep.,2016

@author: thomas.e
'''

import threading
import queue
import datetime


class Logger:
    class Writer:
        def __init__(self, logger, level):
            self.level = level
            self.logger = logger
        
        def write(self, data):
            self.logger.ioQueue.put('[' + self.level + '] ' + data)
            
        def flush(self):
            self.logger.flush()
    
    def __init__(self, logFile):
        self.file = open(logFile, 'a')
        now = '======>>>>>>> Log opened at: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
        self.file.write(now + '\n')
        self.file.flush()
        
        self.ioLock = threading.Lock()
        
        self.ioQueue = queue.Queue()
        
        t = threading.Thread(target=self.ioMonitor, args=(self.ioQueue,))
        t.daemon = True
        t.start()
        
    def ioMonitor(self, ioQueue):
        while True:
            line = self.ioQueue.get()
            self.file.write(line)
            self.flush()
            self.ioQueue.task_done()
            
    def redirecter(self, level):
        return self.Writer(self, level)
    
    def redirect(self, level, stream):
        if stream is None: return
        
        t = threading.Thread(target=self.streamCopy, args=(level, stream))
        t.daemon = True
        t.stream = stream
        t.start()
        
    def streamCopy(self, level, stream):
        for line in iter(stream.readline, b''):
            line = line.decode('utf-8')
            self.ioQueue.put('[' + level + '] ' + line)
            
    def log(self, level, msg):
        self.file.write('[' + level+ '] ' + msg)
        self.flush()
        
    def flush(self):
        self.ioLock.acquire()
        try:
            self.file.flush()
        finally:
            self.ioLock.release()
        
    def shutdown(self):
        self.ioQueue.join()
