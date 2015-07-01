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
python Sum_unique.py --input RIL275_RelocaTEi.CombinedGFF.characterized.unique_mPing.gff 
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def readtable(infile):
    data = defaultdict(lambda : int())
    r = re.compile(r'RIL(\d+)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                ril  = r.search(unit[1]).groups(0)[0] if r.search(unit[1]) else 'NA'
                data[ril] += 1
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    print 'Sample\tUnique_mPing'
    unique_mping = readtable(args.input)
    for ril in sorted(unique_mping.keys(), key=int):
        print 'RIL%s\t%s' %(ril, unique_mping[ril])

if __name__ == '__main__':
    main()

