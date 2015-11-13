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
import gzip

def usage():
    test="name"
    message='''
python vcf2gff.py --input GN22.sv.vcf.gz

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#1       792446  1       N       <DEL>   .       .       SVTYPE=DEL;SVLEN=-1021;END=793467;STRANDS=+-:5;IMPRECISE;CIPOS=-10,10;CIEND=-10,6;CIPOS95=-1,1;CIEND95=-1,0;SU=5;PE=2;SR=3      GT:SU:PE:SR:CN  ./.:5:2:3:2.07
#Chr1    Pindel1 Deletion        1033184 1033405 .       .       .       Size=222;
def vcf2gff(infile, outfile):
    ofile = open(outfile, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'): 
                unit = re.split(r'\t',line)
                data = defaultdict(lambda : str())
                annos = re.split(r';', unit[7])
                for anno in annos:
                    print anno
                    #feature, value = re.split(r'\=', anno)
                    #data[feature] = value
                #if data['SVTYPE'] == 'DEL':
                #    start = unit[1]
                #    end   = data['END']
                #    size  = abs(data['SVLEN'])
                #    print >> ofile, 'Chr%s\tSpeedseq\tDeletion\t%s\t%s\t.\t.\t.\tSize=%s;' %(unit[0], start, end, size)
    ofile.close()


#1       792446  1       N       <DEL>   .       .       SVTYPE=DEL;SVLEN=-1021;END=793467;STRANDS=+-:5;IMPRECISE;CIPOS=-10,10;CIEND=-10,6;CIPOS95=-1,1;CIEND95=-1,0;SU=5;PE=2;SR=3      GT:SU:PE:SR:CN  ./.:5:2:3:2.07
#Chr1    Pindel1 Deletion        1033184 1033405 .       .       .       Size=222;
def vcf2gff_gz(infile):
    cutoff = 500
    filehd = ''
    outfile= ''
    if os.path.splitext(infile)[1] == '.gz':
        outfile = re.sub(r'.vcf.gz', r'.gff', infile)
        filehd  = gzip.open (infile, 'rb')
    else:
        outfile = re.sub(r'.vcf', r'.gff', infile)
        filehd  = open (infile, 'rb')

    ofile = open(outfile, 'w')
    for line in filehd:
        line = line.rstrip()
        if len(line) > 2 and not line.startswith(r'#'): 
            unit = re.split(r'\t',line)
            data = defaultdict(lambda : str())
            annos = re.split(r';', unit[7])
            for anno in annos:
                #print anno
                try:
                    feature, value = re.split(r'=', anno)
                    data[feature] = value
                except:
                    continue
            if data['SVTYPE'] == 'DEL':
                start = unit[1]
                end   = data['END']
                size  = abs(int(data['SVLEN']))
                if size >= cutoff:
                    print >> ofile, 'Chr%s\tSpeedseq\tDeletion\t%s\t%s\t.\t.\t.\tSize=%s;' %(unit[0], start, end, size)
    ofile.close()
    filehd.close()

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

    vcf2gff_gz(args.input)

    #if os.path.splitext(args.input)[1] == '.gz':
    #    outfile = re.sub(r'.vcf.gz', r'.gff', args.input) 
    #    vcf2gff_gz(args.input, outfile)
    #else:
    #    outfile = re.sub(r'.vcf', r'.gff', args.input)
    #    vcf2gff(args.input, outfile)


if __name__ == '__main__':
    main()

