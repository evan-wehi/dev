#!/usr/bin/env python

from os import getenv
from sys import exit, stdout
from os.path import isfile
from subprocess import Popen, PIPE

HOME = getenv('HOME')
PUBLIC_KEY_FILE = HOME + '/.ssh/id_dsa.pub'
PRIVATE_KEY_FILE = HOME + '/.ssh/id_dsa'
AUTHORISED_KEY_FILE = HOME + '/.ssh/authorized_keys'
KNOW_HOST_FILE = HOME + '/.ssh/known_hosts'
HOST_KEY_FILE = '/etc/ssh/ssh_host_ecdsa_key.pub'
SSH_CONFIG = HOME + '/.ssh/config'


def createKeys():
    if isfile(PUBLIC_KEY_FILE) and isfile(PRIVATE_KEY_FILE):
        '''
        Nothing to do.
        '''
        return


    if not isfile(PUBLIC_KEY_FILE) and not isfile(PRIVATE_KEY_FILE):
        CREATE_KEYS_CMD = ['ssh-keygen', '-t', 'dsa',  '-N', '', '-f', PRIVATE_KEY_FILE]
        execute(CREATE_KEYS_CMD)

	# Keys are in an inconsitent state
    if isfile(PRIVATE_KEY_FILE):
        print('The private key exists but the public file does not. Please correct and run again.')

    if isfile(PUBLIC_KEY_FILE):
        print('The public key exists but the private file does not. Please correct and run again.')

	stdout.flush()
	exit(1)

def execute(cmd):
    if type(cmd) is str:
        cmd = cmd.split()

    sp = Popen(cmd,
                stdout=PIPE,
                stderr=PIPE
                )
    sp.wait()
    rc = sp.returncode

    output = sp.stdout.read()
    if rc is not 0:
        print(cmd[0] + ' exited with rc=' + str(rc))
        print('stderr:')
        print(sp.stderr.read())
        print('stdout:')
        print(output)
        exit(rc)

    return output

def authorisePublicKey():
    with file(PUBLIC_KEY_FILE) as pkf:
        publicKey = pkf.read()

    if isfile(AUTHORISED_KEY_FILE):
        with file(AUTHORISED_KEY_FILE) as akf:
            for ak in akf:
                if ak == publicKey:
                    return

    with file(AUTHORISED_KEY_FILE, 'a') as akf:
        akf.write(publicKey)

def configureLenientHostKeyChecking():
    StrictHostKeyChecking = True
    if isfile(SSH_CONFIG):
        with file(SSH_CONFIG) as cf:
            for line in cf:
                words = line.split()
                if words[0].lower() == 'StrictHostKeyChecking'.lower():
                    if words[1].lower() == 'no':
                        StrictHostKeyChecking = False
                    else:
                        StrictHostKeyChecking = True

    if StrictHostKeyChecking:
        with file(SSH_CONFIG, 'a') as cf:
            cf.write('StrictHostKeyChecking no\n')

createKeys()
authorisePublicKey()
configureLenientHostKeyChecking()

