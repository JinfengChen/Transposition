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
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message


def gff_parse(infile):
    data = []
    with open (infile, 'r') as filehd:
        for line in filehd: 
            line = line.rstrip()
            if len(line) > 2 and not line.startswith('#'):
                #print line 
                unit  = re.split(r'\t',line)
                #start = int(unit[3+8]) 
                #end   = int(unit[4+8])
                #chro  = unit[0+8]
                #strand= unit[6+8]
                temp  = defaultdict(str)
                attrs = re.split(r';', unit[17])
                for attr in attrs:
                    #print attr
                    if not attr == '':
                        #print 'yes'
                        idx, value = re.split(r'\=', attr)
                        temp[idx] = value[:-2]
                repid   = temp['ID']
                data.append(repid)
                #print '%s\t%s\t%s\t%s\t%s\t%s' %(repid, chro, start, end, repname, repfam)
                #print '%s\t%s\t%s' %(chro, start, end)
    return data

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Gene    NB      EG4     FoldChange      P-value FDR     NB      HEG4    FoldChange      P-value FDR
#LOC_Os01g01010  1810.32 1390.17 0.77    0.02    0.39    1318.05 1414.15 1.07    0.65    1.00    Up
def expr_parse(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[0]] = line
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gff')
    parser.add_argument('-e', '--expression')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.gff) > 0
    except:
        usage()
        sys.exit(2)

    gene_list = gff_parse(args.gff)
    expr_list = expr_parse(args.expression)
    ofile = open('%s.expression.txt' %(os.path.splitext(args.gff)[0]), 'w')
    for g in gene_list:
        if expr_list.has_key(g):
            print >> ofile, '%s\t%s' %(g, expr_list[g])
        else:
            print >> ofile, '%s\t%s' %(g, 'NA')
    ofile.close()

if __name__ == '__main__':
    main()

