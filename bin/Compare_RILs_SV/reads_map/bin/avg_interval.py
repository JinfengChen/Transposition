#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser

def usage():
    test="name"
    message='''
python avg_interval.py --input mPing_boundary.linked_50Mb_debug2.table_clean.txt


Need to do for this script:
1. when deal with distance, not both mping in each are used only once. The right way is to have three table.
a, number of excision of each mPing
b, distance of closest mPing to this mPing
c, for distance we need to correct using these shared mPing in the RILs to check if the excision happend in these RILs with shared mPing. 
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


#mPing   #Excision
#Chr10_11955070_11955072 1
def read_excision(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'Chr'): 
                unit = re.split(r' |\t',line)
                mping= re.split(r'_',unit[0])
                #print mping, unit[1]
                mping.pop()
                mping = '.'.join(mping)
                #print mping
                if int(unit[1]) < 20:
                    #print mping
                    data[mping] = int(unit[1])
    return data

#Chr3.29404858   Chr3.29404901   43      -       +
def read_distance(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'Chr'): 
                unit = re.split(r'\t',line)
                #print unit[0], unit[1]
                for i in [0, 1]:
                    if data.has_key(unit[i]):
                        data[unit[i]] = int(unit[2]) if int(unit[2]) < data[unit[i]] else data[unit[i]]
                    else:
                        data[unit[i]] = int(unit[2])
    return data



def sum_distance(excision, distance, prefix):
    data1 = defaultdict(lambda : list())
    data1m = defaultdict(lambda : list())
    data2 = defaultdict(lambda : list())
    data2m = defaultdict(lambda : list())
    ofile = open('%s.avg.list' %(prefix), 'w')
    for mping in sorted(excision.keys()):
        #print mping
        if distance.has_key(mping):
            index1 = int(float(distance[mping])/100000)
            index2 = int(float(distance[mping])/10000)
            data1[index1].append(excision[mping])
            data1m[index1].append('%s:%s' %(mping, excision[mping]))
            data2[index2].append(excision[mping])
            data2m[index2].append('%s:%s' %(mping, excision[mping]))
            print >> ofile, '%s\t%s\t%s' %(mping, distance[mping], excision[mping])
    ofile1 = open('%s.sum100kb.txt' %(prefix), 'w')
    print >> ofile1, 'Distance\tNumberofObservation\tExcision_total\tExcision_Average\tExcision_std\tmPing_list\tExcision_list'
    ofile2 = open('%s.sum10kb.txt' %(prefix), 'w')
    print >> ofile2, 'Distance\tNumberofObservation\tExcision_total\tExcision_Average\tExcision_std\tmPing_list\tExcision_list'
    for i in sorted(data1.keys(), key=int):
        #print i
        print >> ofile1, '%s\t%s\t%s\t%s\t%s\t%s\t%s' %((i+1)*100000, len(data1[i]), np.sum(data1[i]), np.mean(data1[i]), np.std(data1[i]), ','.join(map(str,data1[i])), ','.join(data1m[i]))
    for i in sorted(data2.keys(), key=int):
        print >> ofile2, '%s\t%s\t%s\t%s\t%s\t%s\t%s' %((i+1)*10000, len(data2[i]), np.sum(data2[i]), np.mean(data2[i]), np.std(data2[i]), ','.join(map(str,data2[i])), ','.join(data2m[i]))
    ofile.close()
    ofile1.close()
    ofile2.close()

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
    parser.add_argument('--excision')
    parser.add_argument('--distance')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.excision) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output: 
        args.output = 'mPing_boundary.linked_50Mb_debug2.table_clean'
    
    excision = read_excision(args.excision)
    distance = read_distance(args.distance)
    sum_distance(excision, distance, args.output)

    #readtable(args.input)

if __name__ == '__main__':
    main()

