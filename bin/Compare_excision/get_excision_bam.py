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
python get_excision_bam.py --input mping.excision.draw.example

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


#Chr1:6806761-6806763	10,117,130,134,207,44,264,78,22
#Chr1:36270511-36270513	122,125,127,63,98,89
def readtable(infile):
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = unit[0]
                rils  = re.split(r',', unit[1])
                for ril in rils:
                    bam = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam/GN%s.bam' %(ril)
                    subbam(mping, ril, bam)

def subbam(mping, ril, bam):
    reg = re.compile(r'(Chr\d+):(\d+)\-(\d+)')
    match = reg.search(mping)
    start = int(match.groups(0)[1])
    end   = int(match.groups(0)[2])
    chro    = match.groups(0)[0]
    region  = '%s:%s-%s' %(chro, start-1000, end+1000)
    outdir  = './test_example_draw/%s' %(mping)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    test_bam = '%s/%s_%s.bam' %(outdir, ril, mping)
    os.system('samtools view -hb %s %s > %s' %(bam, region, test_bam))
    os.system('samtools index %s' %(test_bam))

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

