#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/software/ProgramPython/lib')
from utility import gff_parser

def usage():
    test="name"
    message='''
python avg_interval.py --input mPing_boundary.linked_50Mb_debug2.table_clean.txt
 
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Chr1.10903901:Chr1.11226135     322234  -       +       136     4       128     1       3       0
def readtable(infile):
    data1 = defaultdict(lambda : list())
    data2 = defaultdict(lambda : list())
    ofile = open('%s.avg.list' %(os.path.splitext(infile)[0]), 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                index1 = int(float(unit[1])/100000)
                index2 = int(float(unit[1])/10000)
                num   = []
                for i in [7, 8, 9]:
                    if int(unit[i]) < 20:
                        num.append(int(unit[i])) 
                    #print '%s\t%s\t%s' %(unit[0], unit[1], unit[i])
                #print '%s\t%s' %(index, np.sum(num))
                unit.append(str(np.sum(num)))
                print >> ofile, '\t'.join(unit)
                data1[index1].append(np.sum(num))
                data2[index2].append(np.sum(num))
    ofile.close()
    ofile1 = open('%s.sum100kb.txt' %(os.path.splitext(infile)[0]), 'w')
    ofile2 = open('%s.sum10kb.txt' %(os.path.splitext(infile)[0]), 'w')
    for i in sorted(data1.keys(), key=int):
        print >> ofile1, '%s\t%s\t%s\t%s\t%s' %((i+1)*100000, len(data1[i]), np.sum(data1[i]), np.mean(data1[i]), np.std(data1[i]))
    for i in sorted(data2.keys(), key=int):
        print >> ofile2, '%s\t%s\t%s\t%s\t%s' %((i+1)*10000, len(data2[i]), np.sum(data2[i]), np.mean(data2[i]), np.std(data2[i]))
    ofile1.close()
    ofile2.close()
    #return data


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

