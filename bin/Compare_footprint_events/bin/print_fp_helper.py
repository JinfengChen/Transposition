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
from utility import gff_parser, createdir, writefile

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

#Chr1:29494572-29494574  Chr1:29328836-29328838
def readtable(infile, fp):
    data = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r' |\t', line)
                for i in unit:
                    data.append(i)
                    if fp.has_key(i):
                        writefile('\n'.join(fp[i]), '%s.txt' %(i))
                        print 'python ../print_fp.py --input %s.txt --start 50 --end 80 > %s.txt.1' %(i, i)
    return data

#>Chr10:11955070-11955072
#Nipponbare      ACCAAGTAATCAATTTTAAATATACAAATTTCTATGGCGCATGTGCATCTAATTAGTGTCACTACTACCTAAACCCTAAAAGGTGCTAGTTGGAAACCTCAC
def readsequence(infile):
    data = defaultdict(lambda : list())
    mping = ''
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if line.startswith(r'>'):
                unit = re.sub(r'>', r'', line)
                mping = unit
            else:
                data[mping].append(line)
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

    fp = readsequence('Excision_newpipe_version1.footprint.sequence.txt')
    readtable(args.input, fp)

if __name__ == '__main__':
    main()

