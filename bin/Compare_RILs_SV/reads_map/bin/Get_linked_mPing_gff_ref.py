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
python Get_linked_mPing_gff.py --gff ../input/HEG4.ALL.mping.non-ref.gff --distance ../input/mPing_dist.100kb.list.sorted

Version for default genome, not psedugenome
Get sub gff of these linked mPing with 100kb. We deal with these linked mPing first or give particular notice.
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Chr1    PseudoGenome    Transposable_element    1132975 1133405 -       .       .       ID=Chr1_1132975_1133405;Original_ID=Chr1.1132977.spanners;TE=mping;TSD=TAA;
def sub_gff(gff, mpings):
    mpings_cp = mpings
    link_100kb_gff = '%s.linked_100kb.gff' %(os.path.splitext(gff)[0])
    other_gff      = '%s.other.gff' %(os.path.splitext(gff)[0])
    ofile1 = open(link_100kb_gff, 'w')
    ofile2 = open(other_gff, 'w')
    with open (gff, 'r') as filehd:
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
                mping_idx = '%s.%s' %(unit[0], unit[4])
                if mpings.has_key(mping_idx):
                    del mpings_cp[mping_idx]
                    print >> ofile1, line
                else:
                    print >> ofile2, line
    ofile1.close()
    ofile2.close()
    for m in sorted(mpings_cp.keys()):
        print 'Not in gff: %s' %(m)

#Chr3.29404858   Chr3.29404901   43      -       +
def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping1 = re.split(r'\.', unit[0])
                mping2 = re.split(r'\.', unit[1])
                mping1_idx = '%s.%s' %(mping1[0], str(int(mping1[1]) + 2))
                mping2_idx = '%s.%s' %(mping2[0], str(int(mping2[1]) + 2))
                if data.has_key(mping1_idx):
                    print 'present in more than 1 pairs: %s' %(mping1_idx)
                if data.has_key(mping1_idx):
                    print 'present in more than 1 pairs: %s' %(mping2_idx)
                data[mping1_idx] = 1
                data[mping2_idx] = 1
    print 'linked mPing: %s' %(str(len(data.keys())))
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gff')
    parser.add_argument('-d', '--distance')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.gff) > 0
    except:
        usage()
        sys.exit(2)

    mpings = readtable(args.distance)
    sub_gff(args.gff, mpings)

if __name__ == '__main__':
    main()

