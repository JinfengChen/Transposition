#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import glob
sys.path.append('/rhome/cjinfeng/software/ProgramPython/lib')
from utility import createdir, runjob

def usage():
    test="name"
    message='''
python IGV_mPing_Landrace.py --strain HEG4

Get subbam files of mPing flanking regions and write batch script for IGV plot.

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


def read_gff(gff):
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
                pos = '%s:%s-%s' %(chro, start-flank, end+flank)
                data[re.sub(r'Chr', r'', chro)][start] = [temp['ID'], pos]
    return data


#Chr1    1032974 1232977
def readtable(infile):
    data = []
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                idx  = '%s:%s-%s' %(unit[0], unit[1], unit[2])
                data.append(idx)
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--strain')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.strain) > 0
    except:
        usage()
        sys.exit(2)


    bedtools ='/opt/bedtools/2.17.0-25-g7b42b3b/bin/bedtools' 
    bam2fastq='/rhome/cjinfeng/BigData/software/bam2fastq/bam2fastq-1.1.0/bam2fastq'
    samtools ='/opt/linux/centos/7.x/x86_64/pkgs/samtools/1.2/bin/samtools'
    gff      = '../input/%s.hom.gff' %(args.strain)
    bams     = glob.glob('../input/%s/*.bam' %(args.strain))
    os.system('%s slop -i %s -g /rhome/cjinfeng/BigData/00.RD/seqlib/MSU7.chrlen -b 5000 > %s.mPing_5kb_flank.gff' %(bedtools, gff, args.strain))
    #os.system('bedtools merge -i ../input/Parent.ALL.mPing.100kb_flank.gff > ../input/Parent.ALL.mPing.100kb_flank.merge.table')
    

    #output directory
    outdir_bam = os.path.abspath('%s.mPing_5kb_flank_bam' %(args.strain))
    createdir(outdir_bam)
    outdir_igv = os.path.abspath('%s.mPing_5kb_flank_igv' %(args.strain))
    createdir(outdir_igv)

    #mping region
    mping_regs = os.path.abspath('%s.mPing_5kb_flank.gff' %(args.strain))
    mpings = read_gff(gff)

    #macbook path of files
    igv_snapshot_dir= '/Users/jinfengchen/biocluster/IGV_2.3.0/igv_snapshot_batch/Landrace/%s_snapshot' %(args.strain)
    igv_batch_dir   = '/Users/jinfengchen/biocluster/IGV_2.3.0/igv_snapshot_batch/Landrace/%s_igv' %(args.strain)
    igv_bam_dir     = '/Users/jinfengchen/biocluster/IGV_2.3.0/igv_snapshot_batch/Landrace/%s_bam' %(args.strain)
    mping_gff       = '/Users/jinfengchen/biocluster/IGV_2.3.0/igv_snapshot_batch/Landrace/%s.hom.gff' %(args.strain)
        
    cmd = []
    ofiles = []
    for bam in sorted(bams):
        bam      = os.path.abspath(bam)
        prefix   = os.path.split(bam)[1]
        prefix   = re.sub(r'.bam', r'.mPing_5kb_flank', prefix)
        
        #mping regions
        cmd.append('%s view -hb -L %s %s > %s/%s.bam' %(samtools, mping_regs, bam, outdir_bam, prefix))
        cmd.append('%s index %s/%s.bam' %(samtools, outdir_bam, prefix))
        #igv batch
        ofile = open('%s/%s.igv' %(outdir_igv, prefix), 'w')
        ofiles.append(ofile)
        print >> ofile, 'new'
        print >> ofile, 'snapshotDirectory %s' %(igv_snapshot_dir)
        print >> ofile, 'load %s/%s.bam' %(igv_bam_dir, prefix)
        print >> ofile, 'load %s' %(mping_gff)
        for chro in sorted(mpings.keys(), key=int):
            for pos in sorted(mpings[chro].keys(), key=int):
                mping = 'Chr%s_%s' %(chro, pos)
                print >> ofile, 'goto %s' %(mpings[chro][pos][1])
                print >> ofile, 'snapshot %s.%s.%s.png' %(mping, prefix, mpings[chro][pos][0])
    for ofile in ofiles:
        ofile.close()

    ofile = open('%s_subbam.sh' %(args.strain), 'w')
    print >> ofile, '\n'.join(cmd)
    ofile.close()

    runjob('%s_subbam.sh' %(args.strain), 2)
if __name__ == '__main__':
    main()

