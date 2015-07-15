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

#Chr1    MSU_osa1r7      mRNA    2903    10817   .       +       .       ID=LOC_Os01g01010.1;Name=LOC_Os01g01010.1;Parent=LOC_Os01g01010 Chr1    DHS     DHsites 10537   10645   .       +       .       ID=1;
#
def readtable(infile):
    data = defaultdict(lambda : int())
    rate   = defaultdict(lambda : list())
    region = defaultdict(lambda : int())
    r = re.compile(r'ID=(.*?)\.\d+')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                gene = r.search(unit[8]).groups(0)[0] if r.search(unit[8]) else "NA"
                data[gene] += int(unit[13]) - int(unit[12]) + 1
                region[gene] = int(unit[4]) - int(unit[3]) + 2001

    for g in data.keys():
        rate[g] = [float(data[g])/region[g], data[g]]
    return rate


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

    gene = readtable(args.input)
    for g in sorted(gene.keys()):
        print '%s\t%s\t%s' %(g, gene[g][0], gene[g][1])

if __name__ == '__main__':
    main()

