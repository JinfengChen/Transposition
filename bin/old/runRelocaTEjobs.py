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
python runRelocaTEjobs.py --input ./RelocaTE/GN100/run_these_jobs.sh

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

#RelocaTE/GN100/run_these_jobs.sh
def splitjob(infile):
    data = defaultdict(str)
    s = re.compile(r'step\_(\d+)')
    filepath = os.path.abspath(os.path.dirname(infile))
    part1 = open(filepath + '/runpart1.sh', 'w')
    part2 = open(filepath + '/runpart2.sh', 'w')
    part3 = open(filepath + '/runpart3.sh', 'w')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            m = s.search(line)
            if m:
                step = m.groups(0)[0]
                #print step, filepath
                if int(step) < 3:
                    print >> part1, line
                elif int(step) == 6:
                    print >> part3, line
                else:
                    print >> part2, line
    data['part1'] = filepath + '/runpart1.sh'
    data['part2'] = filepath + '/runpart2.sh'
    data['part3'] = filepath + '/runpart3.sh'
    return data            

def runjob(script, lines):
    cmd = 'perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --maxjob 30 --lines %s --interval 120 --resource walltime=100:00:00 --convert no %s' %(lines, script)
    #print cmd 
    os.system(cmd)

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

    scripts = splitjob(args.input)
    runjob(scripts['part1'], 10000)
    runjob(scripts['part2'], 15)
    runjob(scripts['part3'], 10000)
    #finish_check()

if __name__ == '__main__':
    main()

