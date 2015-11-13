#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
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


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'LOC'): 
                unit = re.split(r'\t',line)
                data[unit.pop(0)] = unit
    return data


def main():
    
    list1 = readtable(sys.argv[1])
    list2 = readtable(sys.argv[2])
    for gene1 in sorted(list1.keys()):
        if list2.has_key(gene1):
            print '%s\t%s\t%s' %(gene1, '\t'.join(list1[gene1]), '\t'.join(list2[gene1]))

if __name__ == '__main__':
    main()

