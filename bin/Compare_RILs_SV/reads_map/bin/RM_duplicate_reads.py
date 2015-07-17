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
python RM_duplicate_reads.py --sam RIL268.sam

Unmapped reads extract from bam may have duplicates, although it is very rare.
We remove these duplicates from sam file before converting to bam
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def rm_sam(infile):
    data  = defaultdict(lambda : int())
    sam   = infile
    sam_t = '%s.temp.sam' %(os.path.splitext(sam)[0])
    os.system('mv %s %s' %(sam, sam_t))
    ofile = open(sam, 'w')
    with open (sam_t, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if line.startswith(r'@'):
                print >> ofile, line
            elif len(line) > 2: 
                unit = re.split(r'\t',line)
                if int(data[unit[0]]) <= 1:
                    data[unit[0]] += 1
                    print >> ofile, line
    ofile.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sam')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.sam) > 0
    except:
        usage()
        sys.exit(2)

    rm_sam(args.sam)

if __name__ == '__main__':
    main()

