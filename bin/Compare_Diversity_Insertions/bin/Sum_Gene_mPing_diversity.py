#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser

def usage():
    test="name"
    message='''
python Sum_Gene_mPing_diversity.py --input RIL275_RelocaTEi.mPing.annotation

Read annotation of mPing, summary how many mPing within 3kb range of gene and list them.
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Chr1.58365.spanners     Chr1    58363   58365   Five_primer_UTR LOC_Os01g01115.1        -293
def read_anno(infile, prefix):
    data = defaultdict(lambda : list())
    flank = 3000
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = '%s:%s-%s' %(unit[1], unit[2], unit[3])
                if abs(int(unit[6])) <= flank:
                    data[unit[5]].append(mping)
    ofile = open ('%s.gene_mping.table' %(prefix), 'w')
    for gene in sorted(data.keys()):
        mping_n = len(data[gene])
        mping_l = ';'.join(data[gene])
        print >> ofile, '%s\t%s\t%s' %(gene, mping_n, mping_l)
    ofile.close()


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
    
    prefix = os.path.splitext(args.input)[0]
    read_anno(args.input, prefix)

if __name__ == '__main__':
    main()

