#!/usr/bin/env python

from os import getenv
from dispy import JobCluster

HOSTFILE = getenv('PBS_NODEFILE')

def stuff(i):
    import socket
    msg = 'host[' + str(i) + ']=' + socket.gethostname()
    return msg

with file(HOSTFILE) as hf:
    hosts = hf.read().splitlines()

jc = JobCluster(stuff, nodes=hosts)

jobs = []
for i in range(len(hosts)):
    job = jc.submit(i)
    jobs.append(job)
    job.id = i

for job in jobs:
    msg = job()
    print(msg)
