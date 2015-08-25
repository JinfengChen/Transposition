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
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Sample  InsertSize      Size_STD        Depth   Mapped_Depth    Mapped_Rate     Non_Ref_mPing   Confident       Evidence_From_One_End   Characterized   Homozygous      Heterozygous    Somatic Ping    Excision
#RIL1    141.39  25.74   6.35    6.19    0.98    336     228     108     228     211     17      0       5       5
def readtable(infile):
    data = defaultdict(lambda : list())
    ofile = open('Excision_transposases.ping_number.txt', 'w')
    #print >> ofile, '#Ping\t#Count\t#Excision_avg\t#Excision_std'
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Sample'): 
                unit = re.split(r'\t',line)
                if not unit[13] == 'NA':
                    data[unit[13]].append(float(unit[14]))
    for ping in sorted(data.keys(), key=int):
        excision_mean = np.mean(data[ping])
        excision_std  = np.std(data[ping])
        print >> ofile, '%s\t%s\t%s\t%s' %(ping, len(data[ping]), excision_mean, excision_std)
    ofile.close()
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    
    readtable('Excision_transposases.table.txt')

if __name__ == '__main__':
    main()

