#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python RIL230.py RILs_275_ping_genotype.table.ping_code.list RILs_230_ping_genotype.table.ping_code.list

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def read_ping(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'GN'): 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = line
    return data

def read_black(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'GN'):
                unit = re.split(r'\t',line) 
                data[unit[0]] = 1
    return data



def main():
    if len(sys.argv) < 3:
        usage()
        sys.exit()   
 
    list1 = read_ping(sys.argv[1])
    list2 = read_black('Bam.Core.list')
    ofile = open(sys.argv[2], 'w')
    for ril in sorted(list1.keys()):
        if not list2.has_key(ril):
            print >> ofile, list1[ril]
    ofile.close()

if __name__ == '__main__':
    main()

