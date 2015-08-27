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
python Sum_linked_mPing_status.py --dir mPing_boundary_mPing --distance ../input/mPing_dist.100kb.list.sorted

--dir: directory of results of mPing_Boundary_Coverage.py
--distance: distance between pairs of mPing and their strand
Chr3.29404858   Chr3.29404901   43      -       +

Summary the presence/absence of mPing of these pairs within 100kb from results of mPing_Boundary_Coverage.py.
The results will be table of format (+ is covered and - is clipped):
Pairs		RILs_count	++	+-	+-	--
mping1:mping2 	N		N	N	N	N
We will manually correct results of mPing_Boundary_Coverage.py and update this results agasin to see if consistent

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Chr1:6806763-6806761,Genotype_Bin,Pseudo_mPing_status_up,Pseudo_mPing_status_down,Ref_mPing_status
#RIL1,NB,clipped,unknown,unknown
#RIL2,HEG4,covered,covered,unknown
def readcsv(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'RIL'): 
                unit = re.split(r',',line)
                if unit[1] == 'HEG4':
                    #both known and consistent
                    if unit[2] == unit[3] and unit[3] != 'unknown':
                        data[unit[0]] = unit[2]
                    #both unknonw
                    elif unit[2] == unit[3] and unit[3] == 'unknown':
                        if unit[4] != 'unknown':
                            ##use reference instead of pseudoref to determine the status
                            if unit[4] == 'covered':
                                data[unit[0]] = 'clipped'
                            elif unit[4] == 'clipped':
                                data[unit[0]] = 'covered'
                        else:
                            data[unit[0]] = 'unknown'
                    #both known and not euqal
                    elif unit[2] != 'unknown' and unit[3] != 'unknown':
                        if unit[4] != 'unknown':
                            ##use reference instead of pseudoref to determine the status
                            if unit[4] == 'covered':
                                data[unit[0]] = 'clipped'
                            elif unit[4] == 'clipped':
                                data[unit[0]] = 'covered'
                        else:
                            data[unit[0]] = 'unknown' 
                    elif unit[2] != 'unknown':
                        if unit[4] != 'unknown':
                            ##use reference instead of pseudoref to determine the status
                            if unit[4] == 'covered' and unit[2] == 'clipped':
                                data[unit[0]] = 'clipped'
                            elif unit[4] == 'clipped' and unit[2] == 'covered':
                                data[unit[0]] = 'covered'
                            else:
                                data[unit[0]] = 'unknown'
                        else:
                            data[unit[0]] = unit[2]
                    elif unit[3] != 'unknown':
                        if unit[4] != 'unknown':
                            ##use reference instead of pseudoref to determine the status
                            if unit[4] == 'covered' and unit[3] == 'clipped':
                                data[unit[0]] = 'clipped'
                            elif unit[4] == 'clipped' and unit[3] == 'covered':
                                data[unit[0]] = 'covered'
                            else:
                                data[unit[0]] = 'unknown'
                        else:
                            data[unit[0]] = unit[3] 
                    else:
                        if unit[4] != 'unknown':
                            ##use reference instead of pseudoref to determine the status
                            if unit[4] == 'covered':
                                data[unit[0]] = 'clipped'
                            elif unit[4] == 'clipped':
                                data[unit[0]] = 'covered'
                        else:
                            data[unit[0]] = 'unknown'
    return data

#RIL3,HEG4,unknown,clipped,covered,Excision,6,ABCEFG*
def readcsv_excision(infile):
    data = defaultdict(lambda : str())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'RIL'): 
                unit = re.split(r',',line)
                if unit[1] == 'HEG4':
                    data[unit[0]] = unit[5]
    return data

def summary(directory, mpings, prefix):
    ofile1 = open ('%s.status.list' %(prefix), 'w')
    ofile2 = open ('%s.table_other.txt' %(prefix), 'w')
    ofile3 = open ('%s.table_clean.txt' %(prefix), 'w')
    for pair in sorted(mpings.keys()):
        #print '%s, %s, %s, %s, %s, %s' %(pair, mpings[pair][0], mpings[pair][1], mpings[pair][2], mpings[pair][3], mpings[pair][4])
        mping1_status = readcsv_excision('%s/%s.matrix.csv' %(directory, mpings[pair][0]))
        mping2_status = readcsv_excision('%s/%s.matrix.csv' %(directory, mpings[pair][1]))
        #effect_ril, unknown, ++, +-, -+, --
        status = [0, 0, 0, 0, 0, 0]
        print >> ofile1, '%s, %s, %s, %s, %s, %s' %(pair, mpings[pair][0], mpings[pair][1], mpings[pair][2], mpings[pair][3], mpings[pair][4])
        for ril in mping1_status.keys():
            ril_status = ''
            if mping1_status[ril] != 'Check' and mping2_status[ril] != 'Check':
                status[0] += 1
                if mping1_status[ril] == 'Insertion' and mping2_status[ril] == 'Insertion':
                    status[2] += 1
                    ril_status = '++'
                elif mping1_status[ril] != 'Excision' and mping2_status[ril] == 'Excision':
                    status[3] += 1
                    ril_status = '+-'
                elif mping1_status[ril] == 'Excision' and mping2_status[ril] != 'Excision':
                    status[4] += 1
                    ril_status = '-+'
                elif mping1_status[ril] == 'Excision' and mping2_status[ril] == 'Excision':
                    status[5] += 1
                    ril_status = '--'
                else:
                    status[1] += 1
                    ril_status = 'unknown'
            else: 
                status[1] += 1
                ril_status = 'unknown'
            print >> ofile1, '%s\t%s' %(ril, ril_status)
        #print >> ofile2, '%s\t%s\t%s\t%s\t%s' %(pair, mpings[pair][2], mpings[pair][3], mpings[pair][4], '\t'.join(map(str, status)))
        if status[0] > 0:
            if float(status[2])/float(status[0]) > 0.1:
                print >> ofile3, '%s\t%s\t%s\t%s\t%s' %(pair, mpings[pair][2], mpings[pair][3], mpings[pair][4], '\t'.join(map(str, status)))
        else:
            print >> ofile2, '%s\t%s\t%s\t%s\t%s' %(pair, mpings[pair][2], mpings[pair][3], mpings[pair][4], '\t'.join(map(str, status)))
    ofile1.close()
    ofile2.close()
    ofile3.close()

#Chr3.29404858   Chr3.29404901   43      -       +
def readtable(infile):
    data = defaultdict(str)
    pair = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping1 = re.split(r'\.', unit[0])
                mping2 = re.split(r'\.', unit[1])
                mping1_idx = '%s_%s_%s' %(mping1[0], mping1[1], str(int(mping1[1]) + 2))
                mping2_idx = '%s_%s_%s' %(mping2[0], mping2[1], str(int(mping2[1]) + 2))
                #mping1_idx = '%s_%s_%s' %(mping1[0], str(int(mping1[1]) + 2), mping1[1])
                #mping2_idx = '%s_%s_%s' %(mping2[0], str(int(mping2[1]) + 2), mping2[1])
                #if data.has_key(mping1_idx):
                #    print 'present in more than 1 pairs: %s' %(mping1_idx)
                #if data.has_key(mping1_idx):
                #    print 'present in more than 1 pairs: %s' %(mping2_idx)
                data[mping1_idx] = 1
                data[mping2_idx] = 1
                index       = '%s:%s' %(unit[0], unit[1])
                pair[index] = [mping1_idx, mping2_idx, unit[2], unit[3], unit[4]]
    #print 'linked mPing: %s' %(str(len(data.keys())))
    return pair


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-dir', '--dir')
    parser.add_argument('-d', '--distance')
    parser.add_argument('-p', '--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.dir) > 0 and len(args.distance)
    except:
        usage()
        sys.exit(2)

    if not args.project:
        args.project = 'mPing_boundary.linked_100kb'

    mpings = readtable(args.distance)
    summary(args.dir, mpings, args.project)

if __name__ == '__main__':
    main()

