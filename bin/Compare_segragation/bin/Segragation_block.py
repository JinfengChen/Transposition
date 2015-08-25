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
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data

#1 222046 68 50 1.36 59 0.0975125381781
def segragation_block(infile):
    last_pos = defaultdict(lambda : int())
    prefix   = os.path.splitext(infile)[0]
    ofile  = open('%s.distort.txt' %(prefix), 'w')
    ofile1 = open('%s.distort.bed' %(prefix), 'w')
    ofile2 = open('%s.bed' %(prefix), 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r' |\t',line)
                unit[0] = re.sub(r'Chr', '', unit[0])
                if last_pos.has_key(int(unit[0])):
                    unit.insert(1, last_pos[int(unit[0])]+1)
                else:
                    unit.insert(1, 0)
                last_pos[int(unit[0])] = int(unit[2])

                #print >> ofile, '\t'.join(map(str, unit))
                idx = 'Chr%s:%s-%s' %(unit[0], unit[1], unit[2])
                heg4_ratio = float(unit[4])/(float(unit[3]) + float(unit[4]))
                print >> ofile2, 'Chr%s\t%s\t%s\t%s\t%s:%s\t+' %(unit[0], unit[1], unit[2], idx, unit[7], heg4_ratio)
                if float(unit[7]) <= 0.05:
                    print >> ofile, '\t'.join(map(str, unit))
                    #idx = 'Chr%s:%s-%s' %(unit[0], unit[1], unit[2])
                    #heg4_ratio = float(unit[4])/(float(unit[3]) + float(unit[4]))
                    print >> ofile1, 'Chr%s\t%s\t%s\t%s\t%s:%s\t+' %(unit[0], unit[1], unit[2], idx, unit[7], heg4_ratio)
                
    ofile.close()
    ofile1.close()

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

    segragation_block(args.input)
  
if __name__ == '__main__':
    main()

