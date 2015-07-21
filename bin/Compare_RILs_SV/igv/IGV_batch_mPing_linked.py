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
python IGV_batch_mPing_linked.py --gff MSU_r7.Pseudo_mPing.linked_100kb.gff > mPing_linked.list


    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Chr1    PseudoGenome    Transposable_element    1132975 1133405 -       .       .       ID=Chr1_1132975_1133405;Original_ID=Chr1.1132977.spanners;TE=mping;TSD=TAA;
def read_pseudo_gff(gff):
    data = defaultdict(lambda : defaultdict(lambda : list()))
    flank= 800
    with open (gff, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                start = int(unit[3]) 
                end   = int(unit[4])
                chro  = unit[0]
                strand= unit[6]
                temp  = defaultdict(str)
                attrs = re.split(r';', unit[8])
                for attr in attrs:
                    if not attr == '':
                        idx, value = re.split(r'\=', attr)
                        temp[idx] = value
                temp['Original_ID'] = re.sub(r'.spanners', r'', temp['Original_ID'])
                temp['Original_ID'] = re.sub(r'\.', r'_', temp['Original_ID'])
                temp['Original_ID'] = re.sub(r'Chr', r'', temp['Original_ID'])
                ids = re.split(r'_', temp['Original_ID'])
                pos = '%s:%s-%s' %(chro, start-flank, end+flank)
                data[ids[0]][ids[1]] = [temp['ID'], pos]
    return data

def read_bam_list(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                bam  = os.path.split(unit[0])[1]
                ril  = re.sub(r'.bam', r'', bam)
                data[ril] = bam
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gff')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.gff) > 0
    except:
        usage()
        sys.exit(2)
    
    if not args.output:
        args.output = 'mPing_linked'

    mping_gff        = '/Users/jinfengchen/biocluster/IGV_2.3.0/igv_snapshot_batch/MSU_r7.Pseudo_mPing.gff'
    ril_bam_dir      = '/Users/jinfengchen/biocluster/IGV_2.3.0/igv_snapshot_batch/RILs_bam'
    #igv_snapshot_dir = '/Users/jinfengchen/biocluster/IGV_2.3.0/igv_snapshot_batch/RILs_SV'
    igv_snapshot_dir = '/Users/jinfengchen/biocluster/IGV_2.3.0/igv_snapshot_batch/RILs_SV/Snapshot_linked_100kb'
    mpings = read_pseudo_gff(args.gff)
    bams   = read_bam_list('bam.list')   
 
    ofiles = []
    for i in range(0, len(bams.keys())):
        #index = i/30
        #ril   = sorted(bams.keys())[i]
        #print '%s\t%s\t%s' %(i, index, ril)
        #ofile = open('%s.%s.igv' %(args.output, index), 'a')
        ril   = sorted(bams.keys())[i]
        ofile = open('%s.%s.igv' %(args.output, ril), 'w')
        ofiles.append(ofile)
        print >> ofile, 'new'
        print >> ofile, 'snapshotDirectory %s' %(igv_snapshot_dir)
        print >> ofile, 'load %s/%s.bam' %(ril_bam_dir, ril)
        print >> ofile, 'load %s' %(mping_gff)
        for chro in sorted(mpings.keys(), key=int):
            for pos in sorted(mpings[chro].keys(), key=int):
                mping = 'Chr%s_%s' %(chro, pos)
                print >> ofile, 'goto %s' %(mpings[chro][pos][1])
                print >> ofile, 'snapshot %s.%s.%s.png' %(mping, ril, mpings[chro][pos][0])
    for ofile in ofiles:
        ofile.close()

if __name__ == '__main__':
    main()

