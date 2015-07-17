#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
import glob
from Bio import SeqIO

def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

def fasta_id(fastafile):
    fastaid = []
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid.append(record.id)
    return fastaid

def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=30G --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-r', '--ref')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
   
    if not args.ref:
        args.ref = '/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/input/MSU_r7.Pseudo_mPing.10kb_flank.fa'
    #split fasta by sequence
    fas = []
    os.system('perl /rhome/cjinfeng/software/bin/fastaDeal.pl --cuts 1 %s' %(args.ref))
    fas = glob.glob('%s.cut/*.fa.*' %(os.path.split(args.ref)[1]))
    for fa in sorted(fas):
        os.system('bwa index %s' %(fa))

    #write shell for bwa
    fqs = glob.glob('%s/RIL*/*_1.fq' %(args.input))
    ofile = open('RIL_bwa.sh', 'w')
    for fq1 in sorted(fqs):
        fq1 = os.path.abspath(fq1)
        fq2 = re.sub(r'_1', r'_2', fq1)
        fq1_dirs = re.split(r'/', fq1)
        ril = fq1_dirs[-2]
        for fa in sorted(fas):
            fa   = os.path.abspath(fa)
            faid = fasta_id(fa)
            print fa, faid
            prefix = '%s.%s' %(ril, faid[0])
            cmd = 'perl /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/step1_Mapping.pl -ref %s -1 %s -2 %s -min 0 --max 500 -cpu 12 --tool bwa --project %s' %(fa, fq1, fq2, os.path.abspath(prefix))
            print >> ofile, cmd
    ofile.close()

    #run bwa
    #runjob('RIL_bwa.sh', 174)
if __name__ == '__main__':
    main()

