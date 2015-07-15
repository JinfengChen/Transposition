#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
import subprocess

def usage():
    test="name"
    message='''
python mPing_Boundary_Coverage.py --bam_ref ../input/RIL_ALL_bam --bam_mping ../input/RIL_ALL_mping_bam

Check read coverage at mPing boundary from bam files of read2reference and read2mping_flanking. We will create matrix
of RILs and mPing, which gives information of each end of mPing.

    '''
    print message

#Chr1    PseudoGenome    Transposable_element    1132975 1133405 -       .       .       ID=Chr1_1132975_1133405;Original_ID=Chr1.1132977.spanners;TE=mping;TSD=TAA;
#Chr1    PseudoGenome    Transposable_element    2642232 2647573 +       .       .       ID=Chr1_2642232_2647573;Original_ID=Chr1.2640500;TE=ping;TSD=TAA;
def id_mapping(infile, mping2ID_0, mping2ID_1, mping2flank):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                start = int(unit[3]) 
                end   = int(unit[4])
                chro  = unit[0]
                strand= unit[6]
                temp  = defaultdict(str)
                attrs = re.split(r';', unit[8])
                for attr in attrs:
                    if not attr == '':
                        idx, value = re.split(r'\=', attr)
                        temp[idx] = value
                temp['Original_ID'] = re.sub(r'.spanners', r'', temp['Original_ID'])
                mping2ID_0[temp['ID']]  = temp['Original_ID']
                mping2flank[temp['ID']] = '%s_%s_%s' %(unit[0], str(int(unit[3])-1), unit[4])
                mping2ID_1[temp['ID']]  = '%s:%s-%s' %(mping2flank[temp['ID']], 10000, 10000+(int(unit[4])-int(unit[3]))+int(len(temp['TSD'])))
                #print '%s\t%s\t%s' %(temp['ID'], mping2ID_0[temp['ID']], mping2ID_1[temp['ID']])
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_ref')
    parser.add_argument('--bam_mping')
    parser.add_argument('--gff_pseudo')
    parser.add_argument('--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.bam_mping) > 0 and len(args.gff_pseudo) > 0
    except:
        usage()
        sys.exit(2)

    if not args.project:
        args.project = 'mPing_boundary'

    #we use mping gff from pseudogenome as key to project to everything 
    mping2ID_0  = defaultdict(lambda : str()) #ID_0 is the mping id from original call in HEG4, Chr1.1132977
    mping2ID_1  = defaultdict(lambda : str()) #ID_1 is the mping id in 10kb_flank, Chr1_1122974_1143405:10000-10430 
    mping2flank = defaultdict(lambda : str()) #flank is the seq id of 10kb_flank, Chr1_1122974_1143405
    id_mapping(args.gff_pseudo, mping2ID_0, mping2ID_1, mping2flank)

    

if __name__ == '__main__':
    main()

