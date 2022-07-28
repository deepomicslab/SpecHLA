"""
Input: the VCF file that saves long InDels
refine the breakpoints, and convert the formate
return breakpoint file
"""

import sys
import re
import os


def convert():
    f = open(output, 'w')
    for line in open(input):
        if line[0] == '#':
            continue
        array = line.strip().split()
        # print (array)
        if array[6] != "PASS":
            continue
        # if re.search("IMPRECISE", line):
        #     continue
        if not re.search('DEL', array[2]) and not re.search('INS', array[2]):
            continue
        if re.search('DUP', array[2]):
            continue
        for_len = re.search('SVLEN=(.*?)$', array[-3])
        sv_len = for_len.group(1)
        if re.search(';', sv_len): 
            for_len = re.search('SVLEN=(.*?);', array[-3])
            sv_len = for_len.group(1)
            
        if abs(int(sv_len)) < 150:
            continue
        # print ('#', sv_len)
        start = int(array[1]) + 1
        if start < 1001:
            start = 1001
        end = start
        seq = array[4][1:]
        if re.search('DEL', array[2]):
            end = start + len(array[3])
            seq = '.'
        if re.search('1/1', array[-1]):
            copy = 2
        elif re.search('0/1', array[-1]):
            copy = 1
        else:
            continue
        if start < 1200 and start == end:
            continue
        result = '%s %s + %s %s + %s %s'%(array[0], start, array[0], end, seq, copy)
        print (result, file = f)
    f.close()

if __name__ == "__main__":  
    input = sys.argv[1]
    output = sys.argv[2]
    if not os.path.isfile(input):
        f = open(output, 'w')
        f.close()
    else:
        convert()

        
         
