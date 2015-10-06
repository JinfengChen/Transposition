#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
sys.path.append('%s/lib' %(os.getcwd()))

def usage():
    test="name"
    message='''
python high_excision_mPing_geneotype.py --input mping.excision.draw.highexcision


    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

def binarySearch(data, val):
    highIndex = len(data)-1
    lowIndex = 0
    while highIndex > lowIndex:
            index = (highIndex + lowIndex) / 2
            sub = int(data[index])
            #print highIndex, index, lowIndex, sub, val
            if data[lowIndex] == val:
                    return [lowIndex, lowIndex]
            elif sub == val:
                    return [index, index]
            elif data[highIndex] == val:
                    return [highIndex, highIndex]
            elif sub > val:
                    if highIndex == index:
                            return sorted([highIndex, lowIndex])
                    highIndex = index
            else:
                    if lowIndex == index:
                            return sorted([highIndex, lowIndex])
                    lowIndex = index
    return sorted([highIndex, lowIndex])

#Convert BIN MAP
'''
""      "GN1"   "GN10"  "GN100" "GN101" "GN102" "GN103" "GN104" "GN105" "GN106" "GN107" "GN108
"0100222046"    1       1       0       0       0       0       1       1       1       1
"0100500860"    1       1       0       0       0       0       1       1       1       1
'''

def convert_MAP(infile):
    rils = []
    data = defaultdict(lambda : defaultdict(lambda: defaultdict(lambda: str)))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            line = line.replace('"', '')
            if line[0:1].isdigit():
                unit = re.split(r'\t',line)
                #print '%s\t%s\t%s' %(chrs, int(unit[0][2:]), str(chr_start[chrs]))
                chrs = 'Chr%s' %(str(int(unit[0][0:2])))
                pos = int(unit[0][2:])
                for i in range(1,len(unit)):
                    #print i, rils[i], chrs, pos
                    ril = rils[i]
                    data[ril][chrs][pos] = unit[i]
            else:
                unit = re.split(r'\t',line)
                unit[0] = 'Position'
                #print unit[0], unit[1], unit[2]
                rils.extend(unit)
   
    #for t in sorted(data['GN204']['Chr10'].keys(), key=int):
    #    print t
    return data

#Convert SNP map
#""      "GN1"   "GN10"  "GN100" "GN101" "GN102" "GN103" "GN104" "GN105" "GN106" 
#"0100021547A"   NA      NA      NA      0       0       0       NA      NA      
#"0100031071A"   NA      1       0       0       0       0       1       1       
#"0100031478C"   1       1       0       0       0       0       NA      1       

def convert_MAP_SNP(infile):
    rils = []
    data = defaultdict(lambda : defaultdict(lambda: defaultdict(lambda: str)))
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            line = line.replace('"', '')
            if line[0:1].isdigit():
                unit = re.split(r'\t',line)
                #print '%s\t%s\t%s' %(chrs, int(unit[0][2:]), str(chr_start[chrs]))
                chrs = 'Chr%s' %(str(int(unit[0][0:2])))
                unit[0] = unit[0][:-1] #remove reference base
                pos = int(unit[0][2:])
                for i in range(1,len(unit)):
                    #print i, rils[i], chrs, pos, unit[i]
                    ril = rils[i]
                    data[ril][chrs][pos] = unit[i]
            else:
                unit = re.split(r'\t',line)
                unit[0] = 'Position'
                #print unit[0], unit[1], unit[2]
                rils.extend(unit)
   
    #for t in sorted(data['GN204']['Chr10'].keys(), key=int):
    #    print t
    return data

#Chr1:6806761-6806763	10,117,130,134,207,44,264,78,22
#Chr1:36270511-36270513	122,125,127,63,98,89
def read_mping_list(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = unit[0]
                rils  = re.split(r',', unit[1])
                data[mping] = rils
                #print mping, rils
    return data

def subbam(mping, ril, bam, output):
    reg = re.compile(r'(Chr\d+):(\d+)\-(\d+)')
    match = reg.search(mping)
    start = int(match.groups(0)[1])
    end   = int(match.groups(0)[2])
    chro    = match.groups(0)[0]
    region  = '%s:%s-%s' %(chro, start-100000, end+100000)
    outdir  = './%s/%s' %(output, mping)
    if not os.path.exists(output):
        os.mkdir(output)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    test_bam = '%s/%s_%s.bam' %(outdir, ril, mping)
    os.system('samtools view -hb %s %s > %s' %(bam, region, test_bam))
    os.system('samtools index %s' %(test_bam))


def snp_idx_genotype(snp_idx, snpmap, ril, chrs):
    snp_gt = []
    for i in snp_idx:
        s = snpmap[ril][chrs][i]
        snp_gt.append(s)
    return snp_gt

def find_SNP_nearby(start, binmap, snpmap, ril, chrs):
    #array = []
    #array.extend(sorted(binmap[ril][chrs].keys(), key=int))
    array_snp = []
    array_snp.extend(sorted(snpmap[ril][chrs].keys(), key=int))
    #mping after last bin on chromosome, return 0 mean genotype unknown
    if int(start) > int(array_snp[-1]):
        return ['NA']*10
    #index = binarySearch(array, int(start))
    index_snp = binarySearch(array_snp, int(start))
    #print index_snp
    snp_idx   = array_snp[(index_snp[0]-10):(index_snp[0]+10)]
    snp_gt    = snp_idx_genotype(snp_idx, snpmap, ril, chrs)
    return snp_idx, snp_gt

    #check five flanking SNP if consistent
    #pos_gt = []
    #backward_snp_idx = array_snp[(index_snp[0]-4):(index_snp[0]+1)]
    #forward_snp_idx  = array_snp[index_snp[1]:(index_snp[0]+4)]
    #gt_1 = snp_type(backward_snp_idx, snpmap, ril, chrs)
    #gt_2 = snp_type(forward_snp_idx, snpmap, ril, chrs)


    #block genotype unknown, return
    #print 't:1 %s' %(index[1])
    #print 't: %s' %(binmap[ril][chrs][array[index[1]]])
    #if str(binmap[ril][chrs][array[index[1]]]) == 'NA':
    #    pos_gt.extend([array[index[1]], 'NA'])
    #else:
    #    #block genotype known, use snp to comfirm
    #    if gt_1[0] == gt_2[0] and str(gt_1[0]) == str(binmap[ril][chrs][array[index[1]]]):
    #        #snp consistent with genotype block
    #        pos_gt.extend([array[index[1]], gt_1[0]])
    #    else:
    #        #snp inconsistent with genotype block
    #        pos_gt.extend([array[index[1]], 'NA'])
    #return pos_gt

def genotype_mping(mpings, snpmap, binmap, outfile, type):
    ofile = open(outfile, 'w')
    for mping in sorted(mpings.keys()):
        print >> ofile, '>%s' %(mping)
        #snp_header   = []
        #ril_genotype = []
        count = 0
        for ril in mpings[mping]:
            count += 1
            #print '%s\t%s' %(mping, ril)
            ril = 'GN%s' %(ril)
            p = re.compile(r'(\w+):(\d+)\-(\d+)')
            m = p.search(mping)
            chrs = ''
            start = 0
            end   = 0
            if m:
                chrs  = m.groups(0)[0]
                start = m.groups(0)[1]
                end   = m.groups(0)[2]
            #find_SNP_nearby(start, binmap, snpmap, ril, chrs)
            snp_idx, snp_gt = find_SNP_nearby(start, binmap, snpmap, ril, chrs)
            #snp_header = snp_idx
            #ril_genotype.append(snp_gt)
            ##output genotype and mPing
            if count == 1:
                snp_idx.insert(10, mping)
                print >> ofile, '%s\t%s' %('SNP', '\t'.join(map(str, snp_idx)))
            snp_gt.insert(10, type)
            print >> ofile, '%s\t%s' %(ril, '\t'.join(map(str, snp_gt))) 
    ofile.close()
            

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-t', '--type')
    parser.add_argument('--bin_map')
    parser.add_argument('--snp_map')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        args.output = 'high_excision_mping.genotype.table.txt'
    if not args.bin_map:
        args.bin_map = 'MPR.geno.bin'
    if not args.snp_map:
        args.snp_map = 'MPR.geno.data'
    if not args.type:
        args.type = 'excision'

    #bin map and snp genotype
    binmap = convert_MAP(args.bin_map)
    snpmap = convert_MAP_SNP(args.snp_map)

    #read mping and rils
    mpings = read_mping_list(args.input)

    #genotyping 
    genotype_mping(mpings, snpmap, binmap, args.output, args.type)


if __name__ == '__main__':
    main()

