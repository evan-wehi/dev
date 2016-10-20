'''
Created on 27Sep.,2016

@author: thomas.e
'''

# HOST = '10.1.17.99'
# import pydevd
# pydevd.settrace(HOST, stdoutToServer=True, stderrToServer=True)

from IOTask import IOTask
import argparse
from dataunits import convertToString, convertToInt
from os import getcwd

parser = argparse.ArgumentParser('Benchmark IO performance with dd')
parser.add_argument('--num_threads', type=int, help='The number of current IO threads (defaults to 1).', default=1)
parser.add_argument('--dir',         type=str, help='The directory to create the temporary files in (defaults to current directory).', default=getcwd())
parser.add_argument('--size',        type=str, help='The total size of each file (defaults to 1GB)', default='1GB')
parser.add_argument('--blksize',     type=str, help='The block size to use (defaults to 10MB)', default='10MB')

args = parser.parse_args()

N = args.num_threads
wd = args.dir
size = convertToInt(args.size)
blksize = convertToInt(args.blksize)
    
ts = []
for i in range(N):
    t = IOTask(size=size, blksize=blksize, dir=wd)
    t.start()
    ts.append(t)

for t in ts:
    t.join()
    r = convertToString(t.readRate())
    w = convertToString(t.writeRate())
    print('read={}/s write={}/s'.format(r, w))