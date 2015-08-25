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
from excision import bamcheck

def usage():
    test="name"
    message='''
python Excision_Landrace.py --ref_gff ../input/HEG4.hom.gff --qry_gff ../input/EG4.hom.gff --ref_bam ../input/HEG4.bam --qry_bam ../input/EG4.bam

    '''
    print message

def decode(flag):
    if int(flag) == 0:
        return 'footprint'
    elif int(flag) == 1:
        return 'clipped'
    elif int(flag) == 4:
        return 'covered'
    else:
        return 'unknown'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_gff')
    parser.add_argument('--qry_gff')
    parser.add_argument('--ref_bam')
    parser.add_argument('--qry_bam')
    parser.add_argument('--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.ref_gff) > 0 and len(args.qry_gff) > 0
    except:
        usage()
        sys.exit(2)

    if not args.project:
        args.project = 'HEG4vsEG4'
 
    #summary
    shared  = 0
    ref_sum = defaultdict(lambda : int())
    qry_sum = defaultdict(lambda : int())
    #bamcheck files
    bamcheck_file_ref = '%s.bamcheck.txt' %(os.path.splitext(os.path.split(args.ref_bam)[1])[0])
    bamcheck_file_qry = '%s.bamcheck.txt' %(os.path.splitext(os.path.split(args.qry_bam)[1])[0])
    #parse gff of mping file
    ref_gff_dict = gff_parser(args.ref_gff) 
    qry_gff_dict = gff_parser(args.qry_gff)
    #get mping from both strain into one dict
    mpings_all   = ref_gff_dict.keys()
    mpings_all.extend(qry_gff_dict.keys())
    mpings_all   = list(set(mpings_all))
    #classify shared and unique mPing, for unique mping check for footprint of excision
    ofile = open('%s.mPing.Excision.table' %(args.project), 'w')
    for mping in mpings_all:
        if ref_gff_dict.has_key(mping) and qry_gff_dict.has_key(mping):
            print >> ofile, '%s\t%s\t%s\tNA\tNA' %(mping, 1, 1)
            shared += 1
        elif ref_gff_dict.has_key(mping):
            mping_format = '%s:%s-%s' %(ref_gff_dict[mping][0][0], ref_gff_dict[mping][0][3], ref_gff_dict[mping][0][4])
            flag = bamcheck(args.qry_bam, mping_format, bamcheck_file_qry, 'qry')
            print >> ofile, '%s\t%s\t%s\t%s\t%s' %(mping, 1, 0, flag, decode(flag))
            ref_sum['unique'] += 1
            if int(flag) == 0:
                ref_sum['footprint'] += 1
        elif qry_gff_dict.has_key(mping):
            mping_format = '%s:%s-%s' %(qry_gff_dict[mping][0][0], qry_gff_dict[mping][0][3], qry_gff_dict[mping][0][4])
            flag = bamcheck(args.ref_bam, mping_format, bamcheck_file_ref, 'ref')
            print >> ofile, '%s\t%s\t%s\t%s\t%s' %(mping, 0, 1, flag, decode(flag))
            qry_sum['unique'] += 1
            if int(flag) == 0:
                qry_sum['footprint'] += 1
        else:
            #will not happen
            print >> ofile, 'NA\t0\t0\tNA\tNA'
    ofile.close()

    #output summary
    ofile = open('%s.mPing.Excision.sum' %(args.project), 'w')
    print >> ofile, 'Number of Reference mPing: %s' %(len(ref_gff_dict))
    print >> ofile, 'Number of Query mPing: %s' %(len(qry_gff_dict))
    print >> ofile, 'Number of Shared mPing: %s' %(shared)
    print >> ofile, 'Number of unique Ref mPing: %s and %s with footprint in Qry' %(ref_sum['unique'], ref_sum['footprint'])
    print >> ofile, 'Number of unique Qry mPing: %s and %s with footprint in Ref' %(qry_sum['unique'], qry_sum['footprint'])
    ofile.close()

if __name__ == '__main__':
    main()

