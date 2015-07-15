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
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Chr1:36267659-36267661  D5      127;143;151;155;156;157;158;160;165;168;169;173;175;198;200;204;207;208;209;210;214;215;229;234;241;243;247;249;251;254;255;264;267;271;274;277
def read_early_events(infile):
    data = defaultdict(lambda : defaultdict(lambda : str()))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t', line)
                #print unit[1]
                rils = re.split(r';', unit[2])
                for r in rils:
                    data[unit[0]][r] = unit[1]
    return data

#mPing   Presence        Absence Unknown Footprint
#Chr10:11955070-11955072 4       0       10      0               102,124,131,147,185,196,232,261,274,73  
#Chr10:13102744-13102746 6       0       10      0               124,131,145,146,147,172,196,232,275,35  
#Chr10:13370541-13370543 4       1       3       0       84      261,264,85
def read_all_events(infile, early):
    data = defaultdict(lambda : defaultdict(lambda : str()))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'Chr'): 
                unit = re.split(r'\t',line)
                if int(unit[2]) > 0:
                    rils = re.split(r',', unit[5])
                    for r in rils:
                        if early.has_key(unit[0]):
                            if not early[unit[0]].has_key(r):
                                data[r][unit[0]] = 1
                        else:
                            data[r][unit[0]] = 1
    ofile = open('mping.excision.non_ref.shared_rils', 'w')
    ofile1 = open('mping.excision.non_ref.shared_rils_dist', 'w')
    for ril in sorted(data.keys()):
        if len(data[ril].keys()) > 0:
            print >> ofile, '%s\t%s\t%s' %(ril, len(sorted(data[ril].keys())), ','.join(sorted(data[ril].keys())))
            mPing_dist(sorted(data[ril].keys()), ril, ofile1)
    ofile.close()
    ofile1.close()
    #return data

#Chr3:28019800-28019802
def mPing_dist(mpings, ril, ofile1):
    data = defaultdict(lambda : list())
    r = re.compile(r'(\w+):(\d+)\-(\d+)')
    for mping in mpings:
        if r.search(mping):
            #print mping
            chro = r.search(mping).groups(0)[0]
            start= r.search(mping).groups(0)[1]
            data[chro].append(int(start))

    #ofile = open('mping.excision.non_ref.shared_rils_dist', 'w')
    for c in data.keys():
        if len(data[c]) > 1:
            starts = sorted(data[c], key=int)
            for i in range(1,len(starts)):
                #print >> ofile1, '%s,%s' %(starts[i],starts[i-1])
                d = starts[i] - starts[i-1]
                #dist\tril\tchr\tmping1\tmping2
                print >> ofile1, '%s\t%s\t%s\t%s\t%s' %(d, ril, c, starts[i-1], starts[i])
    #ofile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
   
    early = read_early_events('mping.excision.non_ref.early_events')
    read_all_events('mping.excision.non_ref.number.no_low_snp_region', early)
     

if __name__ == '__main__':
    main()

