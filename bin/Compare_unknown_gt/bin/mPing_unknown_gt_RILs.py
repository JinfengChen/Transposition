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
python mPing_unknown_gt_RILs.py --input ../input/mping.excision.draw.highexcision.1 --matrix ../input/High_excision_csv_Ping_HEG4_mPing_only/

Find these RILs that have unknown genotype at given mPing loci. We only count these loci on HEG4 block.

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#../input/Bam.Core.blacklist
def read_blacklist(infile):
    data = defaultdict(lambda : str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]): 
                    data[unit[0]] = 1
    return data 

#Chr1:29494572-29494575  133,134,158,185,271,87
def read_mping_table(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = 1
    return data

#Chr5:25474861-25474863,Genotype_Bin,Pseudo_mPing_status_up,Pseudo_mPing_status_down,Ref_mPing_status,Excision_Code,Ping_Number,Ping_Code
#RIL1,HEG4,covered,covered,unknown,Insertion,5,DEFGH
#RIL2,HEG4,unknown,covered,unknown,Insertion,2,CH
def read_mping_gt(infile, blacklist):
    #data = defaultdict(lambda : list())
    data = []
    data1 = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r',',line)
                if unit[1] == 'HEG4' and (unit[5] == 'Unknown' or unit[5] == 'Check'):
                    if blacklist.has_key(unit[0]):
                        continue
                    ril = re.sub(r'RIL', r'', unit[0])
                    data.append('%s\t%s' %(ril, line))
                    data1.append(ril)
    return data, data1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-m', '--matrix')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    ofile0 = open('%s.check.list' %(args.output), 'w')
    ofile1 = open('%s.draw' %(args.output), 'w')
    blacklist = read_blacklist('../input/Bam.Core.blacklist')
    mpings = read_mping_table(args.input)
    for mping in sorted(mpings.keys()):
        mping1 = re.sub(r':', r'_', mping)
        mping1 = re.sub(r'-', r'_', mping1)
        gt     = '%s/%s.matrix.csv' %(args.matrix, mping1)
        mping_gt, rils = read_mping_gt(gt, blacklist)
        print >> ofile0, '>%s %s' %(mping1, len(mping_gt))
        print >> ofile0, '\n'.join(mping_gt)
        print >> ofile1, '%s\t%s' %(mping, ','.join(rils))
    ofile0.close()
    ofile1.close()

if __name__ == '__main__':
    main()

