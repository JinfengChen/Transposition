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
python mping.excision.ril.py --input bamcheck.log

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#ril excision
#1 5
def read_excision(infile):
    data = defaultdict(lambda: int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r' ',line)
                data[unit[0]] = unit[1]
    return data

#Sample  Parental_mPing
#RIL1    152
def read_parental_mping(infile):
    data = defaultdict(lambda: int)
    r = re.compile(r'RIL(\d+)')
    with open (infile, 'r') as filehd:       
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Sample'): 
                unit = re.split(r'\t',line)
                m = r.search(unit[0])
                ril = m.groups(0)[0] if m else 'NA'
                #print ril, unit[0] 
                data[ril] = unit[1]
    return data


#Pings   Ping_Code       RIL
#NA      NA      RIL158
#NA      NA      RIL242
#0       NA      RIL39
def read_ping(infile):
    data = defaultdict(lambda: list)
    r = re.compile(r'RIL(\d+)')
    with open (infile, 'r') as filehd:       
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Ping'): 
                unit = re.split(r'\t',line)
                m = r.search(unit[2])
                ril = m.groups(0)[0] if m else 'NA'
                #print ril, unit[0] 
                data[ril] = unit[0]
    return data

#Sample  InsertSize      Size_STD        Depth   Mapped_Depth    Mapped_Rate     Non_Ref_mPing   Confident       Evidence_From_One_End   Characterized   Homozygous      Heterozygous    Somatic
#RIL1    141.39  25.74   6.35    6.19    0.98    336     228     108     228     211     17      0
def read_mping(infile, ping, p_mping, p_mping_s, excision):
    data = defaultdict(lambda: list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Sample'):
                unit = re.split(r'\t',line)
                ril  = re.sub(r'RIL', r'', unit[0])
                if ping.has_key(ril):
                    unit.append(ping[ril])
                else:
                    unit.append('NA')
                if excision.has_key(ril):
                    unit.append(excision[ril])
                else:
                    unit.append('0')
                if p_mping.has_key(ril):
                    unit.append(p_mping[ril])
                else:
                    unit.append('NA')
                if p_mping_s.has_key(ril):
                    unit.append(p_mping_s[ril])
                else:
                    unit.append('NA')
                print '\t'.join(unit)
            else:
                print '%s\tPing\tExcision\tParental_mPing\tParental_mPing_close' %(line)
    return data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    excision = read_excision('../input/mping.excision.ril.counts')
    ping     = read_ping('../input/RIL275_RelocaTE.sofia.ping_code.table')
    p_mping   = read_parental_mping('RIL275.Parental_mPing.table.txt')
    p_mping_s = read_parental_mping('RIL275.Parental_mPing_linked_100kb.table.txt')
    read_mping('../input/RIL275_RelocaTEi.summary.table', ping, p_mping, p_mping_s, excision)    

if __name__ == '__main__':
    main()

