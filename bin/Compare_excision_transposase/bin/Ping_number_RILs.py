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

python Ping_number_RIL.py --input ping_code

Summary excision based on Ping number and show details of high excision mPing in each ping number categery.

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#Pings   Ping_Code       RIL
#NA      NA      RIL158
#NA      NA      RIL242
#0       NA      RIL39
def read_ping(infile):
    data = defaultdict(lambda: list())
    r = re.compile(r'RIL(\d+)')
    with open (infile, 'r') as filehd:       
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and not line.startswith(r'Ping'): 
                unit = re.split(r'\t',line)
                m = r.search(unit[2])
                ril = m.groups(0)[0] if m else 'NA'
                #print ril, unit[0] 
                data[ril] = [unit[0], unit[1]]
    return data

#read excision
#ril\texcision
def read_excision(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r' |\t',line)
                data[unit[0]] = unit[1]
    return data



def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r' |\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data


def write_ping_RILs(ping_code, excision, prefix):
    ping2rils = defaultdict(lambda : list())
    ping2excision = defaultdict(lambda : list())
    ping_number = [0, 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 2, 3, 4, 5, 6, 7, 8]
    ping_number_rils = defaultdict(lambda : int())
    
    for ril in sorted(ping_code.keys()):
        if not ping_code[ril][0] == 'NA':
            ping_number_rils[ping_code[ril][0]] += 1
        excision_n = 0 if not excision.has_key(ril) else int(excision[ril])
        if ping_code[ril][0] == '1':
            ping2rils[ping_code[ril][1]].append('RIL%s:%s' %(ril, ping_code[ril][1]))
            #excision_n = 0 if not excision.has_key(ril) else excision[ril]
            ping2excision[ping_code[ril][1]].append(excision_n)
        else:
            ping2rils[ping_code[ril][0]].append('RIL%s:%s' %(ril, ping_code[ril][1]))
            ping2excision[ping_code[ril][0]].append(excision_n)
    ofile = open('%s.Ping_code_RILs.list' %(prefix), 'w')
    print >> ofile, 'Ping\t#RIL\tRIL:Ping_Code'
    ofile1= open('%s.Ping_code_RILs_excision.list' %(prefix), 'w')
    print >> ofile1, 'Ping\t#RIL\tExcision_Mean\tExcision_std\tExcisions\tRIL:Ping_Code'
    for ping in ping_number:
        print >> ofile, '%s\t%s\t%s' %(ping, len(ping2rils[str(ping)]), ','.join(ping2rils[str(ping)]))
        if len(ping2excision[str(ping)]) > 0:
            print >> ofile1, '%s\t%s\t%s\t%s\t%s\t%s' %(ping, len(ping2rils[str(ping)]), np.mean(ping2excision[str(ping)]), np.std(ping2excision[str(ping)]), ','.join(map(str, ping2excision[str(ping)])), ','.join(ping2rils[str(ping)]))
            #print >> ofile1, '%s\t%s\t%s\t%s' %(ping, len(ping2rils[str(ping)]), ','.join(map(str, ping2excision[str(ping)])), ','.join(ping2rils[str(ping)]))
        else:
            print >> ofile1, '%s\t%s\tNA\tNA\t%s\t%s' %(ping, len(ping2rils[str(ping)]), ','.join(map(str, ping2excision[str(ping)])), ','.join(ping2rils[str(ping)]))
    ofile.close()
    ofile1.close()
    return ping_number_rils

#Chr1:36267659-36267661  101,112,95,119,66,89,127,143,151
def write_ping_high_excision(infile, ping_code):
    mping = ['total']
    ping_number_excision = defaultdict(lambda : defaultdict(lambda : int()))
    ofile =open('%s.ping_code' %(infile), 'w')
    ofile1 =open('%s.ping_code.list' %(infile), 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                rils = re.split(r',', unit[1])
                codes = []
                mping.append(unit[0])
                for ril in rils:
                    codes.append('RIL%s:%s:%s' %(ril, ping_code[ril][1], ping_code[ril][0]))
                    print >> ofile1, '%s\t%s' %(ping_code[ril][0], ping_code[ril][1])
                    ping_number_excision[ping_code[ril][0]]['total'] += 1
                    ping_number_excision[ping_code[ril][0]][unit[0]] += 1
                print >> ofile, '%s\t%s' %(unit[0], ','.join(codes))
    ofile.close()
    ofile1.close()
    return ping_number_excision, mping

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ping_code')
    parser.add_argument('--high_excision')
    parser.add_argument('--excision')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.ping_code) > 0
    except:
        usage()
        sys.exit(2)

    prefix = 'RIL275.Ping_number'


    #read excision number in ril
    excision = read_excision(args.excision)

    #read ping code
    ping_code = read_ping(args.ping_code)
    ping_number_rils = write_ping_RILs(ping_code, excision, prefix)

    #high excision ping code
    ping_number_excision, mpings = write_ping_high_excision(args.high_excision, ping_code)

    #summary table of ping number and excision
    ofile = open('%s.high_excision.sum.txt' %(prefix), 'w')
    print >> ofile, '#Ping\t#RILs\t%s' %('\t'.join(mpings))
    for ping_n in sorted(ping_number_rils.keys(), key=int):
        line = [ping_n, ping_number_rils[ping_n]]
        #print 'ping number: %s' %(ping_n)
        for mping in mpings:
            #print 'mping: %s' %(mping)
            n = 0 if not ping_number_excision[ping_n][mping] else ping_number_excision[ping_n][mping]
            line.append(n)
            #print 'excision: %s' %(n)
        print >> ofile, '\t'.join(map(str, line))
    ofile.close()

if __name__ == '__main__':
    main()

