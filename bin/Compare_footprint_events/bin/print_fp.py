#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser, createdir

def usage():
    test="name"
    message='''
python print_fp.py --input t.txt --start 40 --end 80

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def readtable(infile, start, end):
    data = defaultdict(lambda : int())
    flag = 0 
    r = re.compile(r'\w+')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r' |\t', line)
                ril  = unit[0]
                part = line[start:end]
                if not ril  == 'Nipponbare' and not ril.startswith(r'RIL') and r.search(line):
                    flag = 1
                if flag == 0:
                    if not data.has_key(part):
                        if not ril == 'Nipponbare':
                            data[part] = 1
                        print '%s\t%s' %('{0:8}'.format(ril), part)
                    flag = 0
                else:
                    print '%s\t%s' %('{0:8}'.format(ril), part)
                    
                
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-s', '--start')
    parser.add_argument('-e', '--end')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    readtable(args.input, int(args.start), int(args.end))

if __name__ == '__main__':
    main()

