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
python Unique_mPing.py --input RIL275_RelocaTEi.CombinedGFF.characterized.gff

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Pings	Ping_Code	RIL
#0	NA	RIL129
#1	C	RIL230
def summary_unique(unique_mping_d, ping_code, ping_number_sum, ping_single_sum):
    data = defaultdict(lambda : list())
    with open (ping_code, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Pings'):
                unit   = re.split(r'\t',line)
                ril_id = re.sub(r'RIL', '', unit[2])
                data[ril_id] = [unit[0], unit[1]]
    sum_ping_number = defaultdict(lambda : list())
    sum_ping_single = defaultdict(lambda : list())
    for ril in unique_mping_d.keys():
        if not data[ril][0] == 'NA':
            sum_ping_number[data[ril][0]].append(unique_mping_d[ril])
        if len(data[ril][1]) == 1:
            sum_ping_single[data[ril][1]].append(unique_mping_d[ril])
    ofile0 = open(ping_number_sum, 'w')
    ofile1 = open(ping_single_sum, 'w')
    for n_ping in sorted(sum_ping_number.keys(), key = int):
        values = map(int, sum_ping_number[n_ping])
        print >> ofile0, n_ping, np.mean(values), np.std(values)
    for s_ping in sorted(sum_ping_single.keys()):
        values = map(int, sum_ping_single[s_ping])
        print >> ofile1, s_ping, np.mean(values), np.std(values)
    ofile0.close()
    ofile1.close()
    return data

#Chr1    RIL231_0        transposable_element_attribute  4228091 4228092 +       .       .       ID=Chr1.4228092.spanners;Strain=RIL231_0;
def unique_mping(infile, overlap_ref_d, overlap_ril_d, output):
    data  = defaultdict(lambda : int())
    r     = re.compile(r'RIL(\d+)_\d+')
    ofile = open(output, 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'):
                unit  = re.split(r'\t', line)
                index = '%s:%s_%s_%s' %(unit[1], unit[0], unit[3], unit[4])
                if not overlap_ref_d.has_key(index) and not overlap_ril_d.has_key(index):
                    print >> ofile, line
                    ril_id = r.search(unit[1]).groups(0)[0] if r.search(unit[1]) else 'NA'
                    data[ril_id] += 1
    ofile.close()
    return data

#Chr1    RIL231_0        transposable_element_attribute  4228091 4228092 +       .       .       ID=Chr1.4228092.spanners;Strain=RIL231_0;
def parse_overlap(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'#'):
                unit = re.split(r'\t', line)
                ril  = re.sub(r'RIL', r'', unit[1])
                ril  = re.sub(r'_\d+', r'', ril)
                data[ril] += 1
    return data

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-r', '--reference')
    parser.add_argument('-p', '--project')

    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
    
    #if the reference use only mPing from HEG4, the parental mean these overlap with HEG4 mPings
    if not args.reference:
        args.reference = 'HEG4.mping.non-ref.gff'

    if not args.project:
        args.project = 'RIL275.Parental_mPing'

    overlap_ref     = 'temp.overlap_ref'
    bed_overlap_ref = 'bedtools intersect -a %s -b %s > %s' %(args.input, args.reference, overlap_ref)
    os.system(bed_overlap_ref)
    overlap_ref_d   = parse_overlap(overlap_ref)
    ofile = open('%s.table.txt' %(args.project), 'w') 
    print >> ofile, "Sample\tParental_mPing"
    for ril in sorted(overlap_ref_d.keys(), key=int):
        print >> ofile, 'RIL%s\t%s' %(ril, overlap_ref_d[ril])
    ofile.close()

if __name__ == '__main__':
    main()

