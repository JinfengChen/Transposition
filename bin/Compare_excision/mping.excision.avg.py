#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
import subprocess

def usage():
    test="name"
    message='''
python mping.excision.avg.py

Summary average number of excision events according to allele frequency in the population.
    '''
    print message


#Chr7:20204538-20204540 1 0.475
#Chr3:12409837-12409839 1 0.39
def readfile(infile):
    data = defaultdict(lambda: list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r' ',line)
                if float(unit[2]) < 0.05:
                    data['0.025'].append(unit[1])
                elif float(unit[2]) < 0.1:
                    data['0.075'].append(unit[1])
                elif float(unit[2]) < 0.15:
                    data['0.125'].append(unit[1])
                elif float(unit[2]) < 0.2:
                    data['0.175'].append(unit[1])
                elif float(unit[2]) < 0.25:
                    data['0.225'].append(unit[1])                
                elif float(unit[2]) < 0.3:
                    data['0.275'].append(unit[1])
                elif float(unit[2]) < 0.35:
                    data['0.325'].append(unit[1])
                elif float(unit[2]) < 0.4:
                    data['0.375'].append(unit[1])
                elif float(unit[2]) < 0.45:
                    data['0.425'].append(unit[1])
                elif float(unit[2]) < 0.5:
                    data['0.475'].append(unit[1])
                else:
                    data['0.525'].append(unit[1])
    ofile = open('mping.excision.avg.table', 'w')
    for k in sorted(data.keys(), key=float):
        avg = mean(map(float ,data[str(k)]))
        var = std(map(float ,data[str(k)]))
        values = ','.join(data[str(k)])
        print >> ofile, k, avg, var, values
    return data
    ofile.close()

def main():
    s = readfile('mping.excision.table')

if __name__ == '__main__':
    main()




