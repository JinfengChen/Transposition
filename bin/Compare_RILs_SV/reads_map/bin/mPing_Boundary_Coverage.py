#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
import subprocess
sys.path.append('%s/lib' %(os.getcwd()))
from excision import bamcheck, bamcheck_simple, convert_MAP, genotyping

def usage():
    test="name"
    message='''
python mPing_Boundary_Coverage.py --bam_ref ../input/RILs_ALL_bam --bam_pseudo ../input/RILs_ALL_unmapped_mping_bam --gff_ref ../input/HEG4.ALL.mping.non-ref.gff --gff_pseudo ../input/MSU_r7.Pseudo_mPing.gff

Check read coverage at mPing boundary from bam files of read2reference and read2mping_flanking. We will create matrix
of RILs and mPing, which gives information of each end of mPing.

    '''
    print message

#Chr1    PseudoGenome    Transposable_element    1132975 1133405 -       .       .       ID=Chr1_1132975_1133405;Original_ID=Chr1.1132977.spanners;TE=mping;TSD=TAA;
#Chr1    PseudoGenome    Transposable_element    2642232 2647573 +       .       .       ID=Chr1_2642232_2647573;Original_ID=Chr1.2640500;TE=ping;TSD=TAA;
def id_mapping(infile, mping2ID_0):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
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
                temp_ids   = re.split(r'_', temp['ID'])
                temp['ID'] = '%s:%s-%s' %(temp_ids[0], temp_ids[1], temp_ids[2])
                temp_oids  = re.split(r'\.', temp['Original_ID'])
                #print temp_oids
                temp['Original_ID'] = '%s:%s-%s' %(temp_oids[0], temp_oids[1], str(int(temp_oids[1])-2))
                mping2ID_0[temp['ID']]  = temp['Original_ID']
                #print '%s\t%s' %(temp['ID'], mping2ID_0[temp['ID']])
    return data

def bamcheck_ref(bam, mping, bamck_file, ril):
    reg     = re.compile(r'(Chr\d+):(\d+)\-(\d+)')
    match   = reg.search(mping)
    chro    = match.groups(0)[0]
    l_start = int(match.groups(0)[1]) - 2
    l_end   = int(match.groups(0)[1]) + 2
    l_mping = '%s:%s-%s' %(chro, l_start, l_end)
    r_start = int(match.groups(0)[2]) - 2
    r_end   = int(match.groups(0)[2]) + 2
    r_mping = '%s:%s-%s' %(chro, r_start, r_end)
    l_flag  = bamcheck_simple(bam, l_mping, bamck_file)
    r_flag  = bamcheck_simple(bam, r_mping, bamck_file)
    return (l_flag, r_flag)

def decode(flag):
    if int(flag) == 0 or int(flag) == 4:
        return 'covered'
    elif int(flag) == 1:
        return 'clipped'
    else:
        return 'unknown'

def decode_gt(gt):
    if int(gt) == 0:
        return 'NB'
    elif int(gt) == 1:
        return 'HEG4'
    else:
        return 'NA'

def sort_mping_chr(mpings):
    data = defaultdict(lambda : defaultdict(lambda : str()))
    for mping in mpings.values():
        reg     = re.compile(r'Chr(\d+):(\d+)\-(\d+)')
        match   = reg.search(mping)
        chro    = match.groups(0)[0]
        start   = match.groups(0)[1]
        end     = match.groups(0)[2]
        data[chro][start] = mping
    sorted_mping = []
    for c in sorted(data.keys(), key=int):
        for s in sorted(data[c].keys(), key=int):
            sorted_mping.append(data[c][s])
    return sorted_mping

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_ref')
    parser.add_argument('--bam_pseudo')
    parser.add_argument('--gff_ref')
    parser.add_argument('--gff_pseudo')
    parser.add_argument('--bin_map')
    parser.add_argument('--project')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.bam_pseudo) > 0 and len(args.bam_ref) > 0
    except:
        usage()
        sys.exit(2)

    if not args.project:
        args.project = 'mPing_boundary'
    if not args.gff_ref:
        args.gff_ref = '../input/HEG4.ALL.mping.non-ref.gff'
    if not args.gff_pseudo:
        args.gff_pseudo = '../input/MSU_r7.Pseudo_mPing.gff'
    if not args.bam_ref:
        args.bam_ref = '../input/RILs_ALL_bam'
    if not args.bam_pseudo:
        args.bam_pseudo = '../input/RILs_ALL_unmapped_mping_bam'
    if not args.bin_map:
        args.bin_map = '../input/MPR.geno.bin'

    #we use mping gff from pseudogenome as key to project to everything 
    mping2ID_0  = defaultdict(lambda : str()) #ID_0 is the mping id from original call in HEG4, Chr1.1132977
    id_mapping(args.gff_pseudo, mping2ID_0)
    
    #bin map and snp genotype
    binmap = convert_MAP(args.bin_map)

    #go through RILs, for each ril determine the status of each mPing
    bamcheck_file_pseudo = '%s.bamcheck_pseudo.txt' %(args.project)
    bamcheck_file_ref    = '%s.bamcheck_ref.txt' %(args.project)
    mping_status      = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str())))
    mping_status_ref  = defaultdict(lambda : defaultdict(lambda : str()))
    mping_bin_gt  = defaultdict(lambda : defaultdict(lambda : str()))
    mping_snp_gt  = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : str())))
    rils = [1, 2, 3, 4]
    for ril in sorted(rils):
        bam_ref    = '%s/GN%s.bam' %(args.bam_ref, ril)
        bam_pseudo = '%s/RIL%s.bam' %(args.bam_pseudo, ril)
        print 'ril: %s, %s' %(ril, bam_pseudo)
        if os.path.isfile(bam_pseudo):
            for mping in sorted(mping2ID_0.keys()):
                #genotype: 0 ref, 1 non_ref, 3 unknown
                genotype = genotyping(ril, mping, binmap)
                l_flag, r_flag = bamcheck_ref(bam_pseudo, mping, bamcheck_file_pseudo, ril)
                ref_flag       = bamcheck(bam_ref, mping, bamcheck_file_ref, ril) 
                #print '%s\t%s\t%s\t%s' %(ril, mping, l_flag, r_flag)
                mping_status[ril][mping2ID_0[mping]]['up']   = decode(l_flag)
                mping_status[ril][mping2ID_0[mping]]['down'] = decode(r_flag)
                mping_status_ref[ril][mping2ID_0[mping]]     = decode(ref_flag)
                mping_bin_gt[ril][mping2ID_0[mping]] = decode_gt(genotype)
        else:
            print 'bam file not found for rils: RIL%s' %(ril)   

    #output matrix into file
    matrix_file = '%s.mping_status.matrix.txt' %(args.project)
    ofile = open(matrix_file, 'w')
    mping_ranked= sort_mping_chr(mping2ID_0)
    #mping names
    print >> ofile, 'mPing\t%s' %('\t'.join(mping_ranked))
    #mping genotype
    #mping status, matirx
    for ril in sorted(mping_status.keys(), key=int):
        inf_line = ['RIL%s' %(ril)]
        for mping in mping_ranked:
            #inf_line.append(mping_status[ril][mping]['up'])
            #inf_line.append(mping_status[ril][mping]['down'])
            status = '%s:%s' %(mping_status[ril][mping]['up'], mping_status[ril][mping]['down'])
            inf_line.append(status)
        print >> ofile, '\t'.join(inf_line)
    ofile.close()

    #output matrix for individual mPing file
    outdir = '%s_mPing' %(os.path.abspath(args.project))
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for mping in mping_ranked:
        mping_name = re.sub(r':', r'_', mping)
        mping_name = re.sub(r'-', r'_', mping_name)
        ofile = open('%s/%s.matrix.txt' %(outdir, mping_name), 'w')
        #print >> ofile, 'RILs\tDepth(X)\tGenotype_Bin\tDistance_SNP5\tGenotype_SNP5\tDistance_SNP3\tGenotype_SNP3\tmPing_status'
        print >> ofile, 'RILs\tGenotype_Bin\tmPing_status'
        for ril in sorted(mping_status.keys(), key=int):
            #status in pseudo and ref genome have different results
            #in pseudogenome, cover mean junction was covered by reads, which indicates insertion
            #in refgenome, cover mean junction was covered by reads, which indicates no insertion or excision
            status     = '%s:%s' %(mping_status[ril][mping]['up'], mping_status[ril][mping]['down'])
            status_ref = mping_status_ref[ril][mping]
            print >> ofile, 'RIL%s\t%s\t%s\t%s' %(ril, mping_bin_gt[ril][mping], status, status_ref)
        ofile.close() 

if __name__ == '__main__':
    main()

