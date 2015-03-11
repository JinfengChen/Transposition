#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO
import subprocess

def usage():
    test="name"
    message='''
python Excision.py --input ../input/10092013.mpings.gff
    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

'''
Chr1    RelocaTE        mPing   10903901        10903903
'''
def readmPing(nonref):
    data = defaultdict(lambda: int)
    #data.update(readfile(ref, 0))
    #print len(data.keys())
    data.update(readfile(nonref, 1))
    #print len(data.keys())
    #data.update(readfile(shared, 2)) 
    #print len(data.keys())
    return data

def readfile(infile, flag):
    data = defaultdict(lambda: int)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = '%s:%s-%s' %(unit[0], unit[3], unit[4])
                data[mping] = flag
    return data

'''
Chr1    RIL10_0 RelocaTE        1132975 1132977 .       .       .       ID=mPing_1;Strain=RIL10_0;TSD=TAA;
'''
def readtable(infile):
    data = defaultdict(lambda: defaultdict(lambda: int()))
    p = re.compile(r'Strain=RIL(\d+)\_\d+;')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = '%s:%s-%s' %(unit[0], unit[3], unit[4])
                m = p.search(unit[8])
                strain = m.groups(0)[0] if m else 'NA'
                #print strain, mping, line
                data[strain][mping] = 1
    return data


'''
Chr1    RIL10_0 RelocaTE        1132975 1132977 .       .       .       ID=mPing_1;Strain=RIL10_0;TSD=TAA;
'''
def mping_frequency(infile):
    data = defaultdict(lambda: defaultdict(lambda: int))
    rils = defaultdict(int)
    inf  = defaultdict()
    p = re.compile(r'Strain=RIL(\d+)\_\d+;')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = '%s:%s-%s' %(unit[0], unit[3], unit[4])
                m = p.search(unit[8])
                strain = m.groups(0)[0] if m else 'NA'
                data[mping][strain] = 1
                rils[strain] =1
                inf[mping] = [unit[0], unit[3], unit[4]]
    total = len(rils.keys())
    #print total
    mping_frq = defaultdict(lambda: int)
    for m in data.keys():
        count = len(data[m].keys())
        frq = float(count)/total
        #print '%s\t%s\t%s\t%s\t%s\t%s\t%s' %(inf[m][0], inf[m][1], inf[m][2], m, '+', str(count), str(frq))
        mping_frq[m] = frq
    return mping_frq


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
 

def findbin(start, binmap, ril, chrs):
    
    array = []
    array.extend(sorted(binmap[ril][chrs].keys(), key=int))
    #mping after last bin on chromosome, return 0 mean genotype unknown
    if int(start) > int(array[-1]):
        return 0
    index = binarySearch(array, int(start))
    #bin overlap with mping is a transition bin on recombination block, which mean genotype may not be precise, return 0
    #if transition[ril][chrs].has_key(array[index[1]]):
    #    return 0
    #return the bin that overlap with mping
    #if int(start) == 6566243 and ril == 'GN204':
    #    for t in sorted(binmap[ril][chrs].keys(), key=int):
    #        print t
    #print start, array[index[0]], array[index[1]]
    return array[index[1]]

def genotyping(ril, mping, binmap):
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
    #print ril, chrs, start, len(binmap[ril][chrs].keys())
    pos = findbin(start, binmap, ril, chrs)
    #print 'find bin', ril, chrs, start, pos
    #pos is 0 when genotype of bin that overlap with mping is unknown
    #we should use raw bin not filled or uniq bin here and check if the genotype is NA.
    #genotype = binmap[ril][chrs][pos] if pos > 0 else '3'
    if pos > 0 and binmap[ril][chrs][pos] != 'NA':
        genotype = binmap[ril][chrs][pos]
    else:
        genotype = '3'
    return genotype

def validmap(binmap):
    last = 0
    data = defaultdict(lambda : defaultdict(lambda: defaultdict(lambda: int))) 
    for ril in sorted(binmap.keys()):
        for chrs in sorted(binmap[ril].keys()):
            for pos in sorted(binmap[ril][chrs].keys(), key=int):
                if binmap[ril][chrs][pos] != binmap[ril][chrs][last]:
                    data[ril][chrs][last] = 1
                    data[ril][chrs][pos] = 1
                    last = pos
                    #print ril, chrs, pos
    return data

def bamcheck(bam, mping):
    reg = re.compile(r'(Chr\d+):(\d+)\-(\d+)')
    match = reg.search(mping)
    start = int(match.groups(0)[1])
    end   = int(match.groups(0)[2])
    cmd = '/usr/local/bin/samtools view %s %s' %(bam, mping)
    ofile = open('mping.excision.bamcheck', 'a') 
    print >> ofile, bam, mping
    print >> ofile, cmd
    out = subprocess.check_output(cmd, shell=True)
    lines = re.split(r'\n', out)
    total = 0
    covered = 0
    clipped = 0
    footprint = 0
    flag = 0 # flag indicate there are indel and softclips in the reads, probably solfclip only, supports for mping insertion
    if len(lines) < 2:
        return 2
    pattern = re.compile(r'([0-9]+)([A-Z]+)')
    for line in lines:
        print >> ofile, line
        total += 1
        unit = re.split(r'\t', line)
        if len(unit) < 2:
            continue
        matches = []
        indel = 0
        soft  = 0
        for (base, match) in re.findall(pattern, unit[5]):
            if match == 'I' or match == 'D':
                indel += 1
                continue
            elif match == 'S' and int(base) > 10:
                soft += 1
                matches.append([base, match])
            else:
            #print >> ofile, base, match
                matches.append([base, match])
            #print >> ofile, matches[0]
        # merge neighbor matches that have some types after removing I and D
        j = 0
        for i in range(1, len(matches)):
            i = i - j
            if matches[i][1] == matches[i-1][1]:
                matches[i-1][0] = matches[i][0] + matches[i-1][0]
                del matches[i]
                j = j + 1
        #print >> ofile, len(matches)
        flank = 20
        if indel > 0 and soft > 0:
            flag += 1
        elif indel > 0:
            footprint += 1
        if len(matches) == 1 and matches[0][1] == 'M': # perfect match
            if int(unit[3]) < start - flank and int(unit[3]) + int(matches[0][0]) > end + flank: # read cover mping insertion site
                covered += 1
        elif len(matches) == 2: # may have soft clip
            if matches[0][1] == 'S' and int(matches[0][0]) > 10 and matches[1][1] == 'M' and int(matches[1][0]) > 10: # clip at start
                if int(unit[3]) > start - flank and int(unit[3]) < end + flank: # read at mping insertion site, so the clip should be due to mping insertion
                    covered += 1
                    clipped += 1
                elif int(unit[3]) <= start - flank and int(unit[3]) + int(matches[1][0]) >= end + flank: # read cover mping insertion site
                    covered += 1
            elif matches[1][1] == 'S' and int(matches[1][0]) > 10 and matches[0][1] == 'M' and int(matches[0][0]) > 10: # clip at end
                if int(unit[3]) + int(matches[0][0]) > start - flank and int(unit[3]) + int(matches[0][0]) < end + flank: # read at mping insertion site, so the clip should be due to mping insertion
                    covered += 1
                    clipped += 1
                elif int(unit[3]) <= start - flank and int(unit[3]) + int(matches[0][0]) >= end + flank: # read cover mping insertion site
                    covered += 1
        elif len(matches) == 3 and matches[0][1] == 'S' and matches[1][1] == 'M' and matches[2][1] == 'S': # may have soft clip, but the other end of reads have clip too
            if int(unit[3]) > start - flank and int(unit[3]) < end + flank and int(matches[0][0]) > 10: # read at mping insertion site, so the clip should be on the left if due to mping
                covered += 1
                clipped += 1
            elif int(unit[3]) + int(matches[1][0]) > start - flank and int(unit[3]) + int(matches[1][0]) < end + flank: # read start before mping insertion site, but clipped at mping
                covered += 1
                clipped += 1
            elif int(unit[3]) <= start - flank and int(unit[3]) + int(matches[1][0]) >= end + flank: # read cover the mping insertion site
                covered += 1
        print >> ofile, covered, clipped
    print >> ofile, covered, clipped
    print >> ofile, total, footprint, flag, float(float(footprint)/total)
    rate = float(float(clipped)/covered) if covered > 0 else 0
    if total <= 2:
        print >> ofile, 2
        return 2 # coverage too low
    elif float(float(clipped)/total) > 0.3:
        print >> ofile, 1
        return 1 # support mping insertion with clipped reads
    elif total > 2 and clipped < 1 and covered < 1 and flag < 1 and footprint < 1:
        print >> ofile, 2 # no supporting reads for anything
        return 2 # 
    #elif float(float(footprint)/total) > 0.3 and flag == 0:
    #    print >> ofile, 0
    #    return 0 # support mping excision, indel is possibly footprint of excision
    elif float(float(covered)/total) > 0.3:
        print >> ofile, 0
        if footprint >= 3:
            return 0 # support mping excision, many read cover the breapoint suggests precise excision, footprint
        else:
            return 4 # support mping excision, many read cover the breapoint suggests precise excision, perfect
    else:
        print >> ofile, 3
        return 3 # other case
def excision(mPing_ancestor, mPing_rils, mPing_frq, num_file):
    mPing_excision = defaultdict(lambda : defaultdict(lambda : int()))
    binmap = convert_MAP('MPR.geno.bin')
    #transition = validmap(binmap)
    ofile_num = open(num_file, 'w')
    print >> ofile_num, 'mPing\tPresence\tAbsence\tUnknown\tFootprint'
    for mping in sorted(mPing_ancestor.keys()):
        #print mping
        num_present   = 0
        num_unsure    = 0
        num_absent    = 0
        num_excision  = 0
        num_footprint = 0
        ril_absent    = []
        ril_unsure    = []
        ril_footprint = []
        if not mPing_frq.has_key(mping) or mPing_frq[mping] < 0.05:
            continue
            #print 'not in rils'
        #print 'in rils'
        for ril in sorted(mPing_rils.keys()):
            #print 'TT',ril,len(binmap[ril].keys())
            if not binmap.has_key('GN%s' %(ril)):
                continue
            genotype = genotyping(ril, mping, binmap)
            bam = '/rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam/GN%s.bam' %(ril)
            flag = 2
            #print 'Check:', mping, ril, genotype
            #if mping == 'Chr1:1715117-1715119' and int(ril) == 10:
                #print 'EX', mping, ril, genotype, mPing_rils[ril][mping], mPing_ancestor[mping]
            #if (int(genotype) == 0 and int(mPing_ancestor[mping]) != 1): # ril has genotype of reference and mping has genotype of reference or shared
                #print 'S1', mPing_rils[ril].has_key(mping)
                #if os.path.isfile(bam):
                #    flag = bamcheck(bam, mping)
            #    if not mPing_rils[ril].has_key(mping):
            #        if os.path.isfile(bam):
            #            flag = bamcheck(bam, mping)
            #            if flag == 0: ## bam check showed no mping insertion in this ril
            #                mPing_excision[mping][ril] = 1
            #            elif flag == 4:
            #                mPing_excision[mping][ril] = 2
                    #print mPing_excision[mping]
            #if (int(genotype) == 1 and int(mPing_ancestor[mping]) != 0): # ril has genotype of nonref and mping has genotype of nonref or shared
                #print 'S2', mPing_rils[ril].has_key(mping)
                #if os.path.isfile(bam):
                #    flag = bamcheck(bam, mping)
            #    if not mPing_rils[ril].has_key(mping):
            #        if os.path.isfile(bam):
            #            flag = bamcheck(bam, mping)
            #            if flag == 0: ## bam check showed no mping insertion in this ril
            #                mPing_excision[mping][ril] = 1
            #            elif flag == 4:
            #                mPing_excision[mping][ril] = 2
            if (int(genotype) == 1 and int(mPing_ancestor[mping]) == 1): # ril has genotype of nonref and mping has genotype of nonref. Check excision for non_reference insertion
                #print 'S2', mPing_rils[ril].has_key(mping)
                if not mPing_rils[ril].has_key(mping):
                    if os.path.isfile(bam):
                        flag = bamcheck(bam, mping)
                        if flag == 0: ## bam check showed no mping insertion in this ril
                            mPing_excision[mping][ril] = 1
                            num_absent += 1
                            ril_absent.append(ril)
                            num_footprint += 1
                            ril_footprint.append(ril)
                        elif flag == 4:
                            mPing_excision[mping][ril] = 2
                            num_absent += 1
                            ril_absent.append(ril)
                        elif flag == 2:
                            ##coverage too low
                            num_unsure += 1
                            ril_unsure.append(ril)
                        elif flag == 3:
                            ##other case
                            num_unsure += 1
                            ril_unsure.append(ril)
                        elif flag == 1:
                            ##have insertion
                            num_present += 1
                    #print mPing_excision[mping]
        print >> ofile_num, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(mping, num_present, num_absent, num_unsure, num_footprint, ','.join(ril_absent), ','.join(ril_unsure), ','.join(ril_footprint))
    return mPing_excision

def preplot(frq, ancestor, excision):
    ofile = open('mping.excision.table', 'w')
    for m in excision.keys():
        print '>', m, len(excision[m].keys())
        print >> ofile, m, len(excision[m].keys()), frq[m]
        for i in sorted(excision[m].keys()):
            print i, excision[m][i]
    ofile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
    mPing_ancestor = readmPing('HEG4.ALL.mping.non-ref.gff')
    #mPing_ancestor = readmPing('HEG4.mping.non-ref.gff')
    mPing_rils     = readtable(args.input)
    mPing_frq      = mping_frequency(args.input)
    mPing_excision = excision(mPing_ancestor, mPing_rils, mPing_frq, 'mping.excision.number')
    preplot(mPing_frq, mPing_ancestor, mPing_excision)

if __name__ == '__main__':
    main()

