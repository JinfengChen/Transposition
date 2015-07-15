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
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Chr1    MSU_osa1r7      mRNA    2903    10817   .       +       .       ID=LOC_Os01g01010.1;Name=LOC_Os01g01010.1;Parent=LOC_Os01g01010
#Chr1	0	200000	Win0	100	+
def readtable(infile):
    data = defaultdict(lambda : int())
    r = re.compile(r'ID=(.*?)\.\d+')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                gene = r.search(unit[8]).groups(0)[0] if r.search(unit[8]) else "NA"
                tss  = unit[3] if unit[6] == '+' else unit[4]
                start = int(tss) - 5000 if int(tss) - 5000 > 0 else 0
                end   = int(tss) + 5000
                print '%s\t%s\t%s\t%s\t%s\t%s' %(unit[0], start, end, gene, 100, unit[6])

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

    readtable(args.input)

if __name__ == '__main__':
    main()

