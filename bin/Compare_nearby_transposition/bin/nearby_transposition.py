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
import glob

def usage():
    test="name"
    message='''
python nearby_transposition.py --frequency ../input/mping.ril.frquency

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Chr9    7356518 7356520 Chr9:7356518_7356520    +       1       0.00363636363636
def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                index= '%s:%s-%s' %(unit[0], unit[1], unit[2])
                data[index] = unit[6]
    return data

#Chr1    2129220 6350512 Non_Ref
def readdist(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                index= '%s:%s-%s' %(unit[0], unit[1], int(unit[1])+2)
                data[index] = int(unit[2])
    return data

def sum_dist(dist, frq):
    mping_new    = []
    mping_shared = []
    for mping in sorted(dist.keys()):
        if frq.has_key(mping):
            if float(frq[mping]) < 0.01:
                mping_new.append(dist[mping])
            else:
                mping_shared.append(dist[mping])
        else:
            print '%s: frquency not avialiable' %(mping)
    return [len(mping_new), np.mean(mping_new), len(mping_shared), np.mean(mping_shared)]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--frequency')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.frequency)
    except:
        usage()
        sys.exit(2)

    frq    = readtable(args.frequency)
    rils   = glob.glob('../../RIL275_RelocaTEi/RelocaTEi_GN*')
    for ril in rils:
        ril_n = os.path.split(ril)[1]
        ril_n = re.sub(r'RelocaTEi_GN', r'', ril_n)
        #gff = '%s/repeat/results/ALL.all_nonref_insert.characTErized.gff' %(ril)
        gff = 'RIL%s.mping.gff' %(ril_n)
        os.system('grep "RIL%s_" ../input/RIL275_RelocaTE.CombinedGFF.Ref_only.gff >> %s' %(ril_n, gff))
        os.system('grep "RIL%s_" ../input/RIL275_RelocaTE.CombinedGFF.Shared.gff >> %s' %(ril_n, gff))
        os.system('grep "RIL%s_" ../input/RIL275_RelocaTE.CombinedGFF.Non_ref.gff >> %s' %(ril_n, gff))
        print ril_n
        print gff
        os.system('perl mPing_dist_ref.pl --input %s > RIL%s.dist.txt' %(gff, ril_n))
        ril_dist = readdist('RIL%s.dist.txt' %(ril_n))
        summary = sum_dist(ril_dist, frq)
        print '%s' %('\t'.join(map(str, summary)))

if __name__ == '__main__':
    main()

