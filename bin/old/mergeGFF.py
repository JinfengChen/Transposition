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
python mergeGFF.py --input ./RelocaTE/

    '''
    print message


def getlines(infile):
    data = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if not line.startswith(r'#'):
                data.append(line)
    return data 


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
    #./RelocaTE/GN93/mping/results/GN93.mping.all_inserts.gff
    for ril in os.listdir(args.input):
        ril_dir = '%s/%s' %(args.input, ril)
        filename = '%s/%s.transposition.gff' %(ril_dir, ril)
        filename_s = '%s/%s.transposition.sort.gff' %(ril_dir, ril)
        ofile = open (filename, 'w')
        for te in os.listdir(ril_dir):
            te_dir = '%s/%s' %(ril_dir, te)
            #te_results = '%s/results' %(te_dir)
            #print te_results
            #if not os.path.isdir(te_results):
            #    ofile.close()
            #    os.system('rm %s' %(filename))
            #    continue
            te_gff = '%s/results/%s.%s.all_inserts.gff' %(te_dir, ril, te)
            if os.path.exists(te_gff):
                #print te_gff
                data = getlines(te_gff)
                if len(data) > 0:
                    gff = '\n'.join(data)
                    print >> ofile, gff
        ofile.close()
        sort = 'sort -k1,1n -k3,3n %s > %s' %(filename, filename_s)
        os.system(sort)

if __name__ == '__main__':
    main()

