'''
Created on 30Sep.,2016

@author: thomas.e
'''

from math import floor, pow, log
from re import compile, match

__units__ = ("KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
__pattern__ = compile(r'([+-]{0,1})([0-9]+)\ *([a-zA-Z]{0,2})')

def convertToString(size):
    if (size == 0):
        return '0B'
    i = int(floor(log(size,1024)))
    p = pow(1024,i)
    s = round(size/p,1)
    return '%s %s' % (s,__units__[i])

def convertToInt(size):
    result = match(__pattern__, size)

    sign = -1 if result.group(1) == '-' else 1    
    size = int(result.group(2)) * sign
    unit = result.group(3)
    
    if unit == '':
        return size

    if len(unit) == 1:
        raise Exception('Invalid unit: ' + unit)
         
    try:
        i = __units__.index(unit.upper())
    except ValueError:
        raise Exception('Invalid unit: ' + unit)
    
    return int(size * pow(2, (i+1)*10))
    
if __name__ == '__main__':
    import sys
    v = sys.argv[1]
    print(convertToInt(v))