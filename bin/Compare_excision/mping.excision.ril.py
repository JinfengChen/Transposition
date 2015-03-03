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
python mping.excision.ril.py --input bamcheck.log

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#> Chr7:20204538-20204540 1
#113 2
#> Chr3:12409837-12409839 1
def readlog(infile, mping):
    data = defaultdict(lambda: int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if not line.startswith(r'>'): 
                unit = re.split(r' ',line)
                data[unit[0]] += 1
    ofile = open('mping.excision.ril.counts', 'w')
    for ril in sorted(data.keys(), key=int):
        print >> ofile, ril, data[ril]
    ofile.close()

    ofile = open('mping.excision.ril.transposition', 'w')
    for ril in sorted(mping.keys(), key=int):
        if data.has_key(ril):
            print >> ofile, '%s\t%s\t%s\t%s\t%s' %(ril, mping[ril][0], mping[ril][1], mping[ril][2], data[ril])
    ofile.close()

#RIL1_0  184     34
def readunique(infile):
    data = defaultdict(lambda: list)
    r = re.compile(r'RIL(\d+)\_\d+')
    with open (infile, 'r') as filehd:       
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                m = r.search(unit[0])
                ril = m.groups(0)[0] if m else 'NA'
                data[ril] = [unit[1], unit[2], unit[3]]
    return data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
    mping = readunique('../input/unique.mping')
    readlog(args.input, mping)

if __name__ == '__main__':
    main()

