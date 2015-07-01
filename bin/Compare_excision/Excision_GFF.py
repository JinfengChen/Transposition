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
python Excision_GFF.py --input excision.table --gff HEG4_mPing.gff
Read excision.table and get sub gff from HEG4_mPing.ff

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def gff_parse(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd: 
            line = line.rstrip()
            if len(line) > 2 and not line.startswith('#'): 
                unit  = re.split(r'\t',line)
                start = int(unit[3]) 
                end   = int(unit[4])
                chro  = unit[0]
                strand= unit[6]
                temp  = defaultdict(str)
                attrs = re.split(r';', unit[8])
                for attr in attrs:
                    #print attr
                    if not attr == '':
                        #print 'yes'
                        idx, value = re.split(r'\=', attr)
                        temp[idx] = value
                index = '%s:%s-%s' %(unit[0], start, end)
                data[index] = unit
                #print '%s\t%s\t%s\t%s\t%s\t%s' %(repid, chro, start, end, repname, repfam)
                #print '%s\t%s\t%s' %(chro, start, end)
    return data



#Chr2:21319934-21319936 3 0.48
def readtable(infile, gff):
    data = defaultdict(str)
    outfile = '%s.gff' %(os.path.splitext(infile)[0])
    ofile = open(outfile, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r' ',line)
                if gff.has_key(unit[0]):
                    idx = '%s_%s' %(gff[unit[0]][0], gff[unit[0]][4])
                    attrs = re.split(r';', gff[unit[0]][8])
                    flag = 0
                    for attr in attrs:
                        if attr == 'ID':
                            flag = 1
                    if flag == 0:
                        gff[unit[0]][8] = 'ID=%s;%s' %(idx, gff[unit[0]][8])
                    print >> ofile, '\t'.join(gff[unit[0]])
                else:
                    print 'not found: %s' %(unit[0])
    ofile.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-g', '--gff')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    refgff = gff_parse(args.gff)
    readtable(args.input, refgff)
    

if __name__ == '__main__':
    main()

