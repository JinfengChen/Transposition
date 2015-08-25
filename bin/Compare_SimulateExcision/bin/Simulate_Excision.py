#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
sys.path.append('/rhome/cjinfeng/BigData/software/ProgramPython/lib')
from utility import gff_parser, createdir
import random

def usage():
    test="name"
    message='''
python Simulate_Excision.py --input ../input/RIL275_RelocaTEi.CombinedGFF.characterized.unique_mPing.gff


    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid

##read sample mPing from gff
def readtable(infile):
    data = defaultdict(lambda : str())
    rank = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                idx  = '%s:%s-%s' %(unit[0], unit[4], unit[5])
                rank += 1
                data[int(rank)] = line
    return data

def writefile(lines, filename):
    ofile = open(filename, 'w')
    print >> ofile, lines
    ofile.close()

##plot simulation compared with sample started
def distr_gff(gff, prefix):
    os.system('bedtools intersect -a %s -b MSU_r7.all.final.full.utr.gff3 -wao > %s.intersect' %(gff, prefix))
    os.system('python mPing_position.py --input %s.intersect' %(prefix))
    os.system('python mPing_intron.py --input %s.intersect' %(prefix))
    os.system('bedtools closest -a %s -b MSU_r7.all.final.mRNA.gff -d > %s.mRNA.intersect' %(gff, prefix))
    os.system('python mPing_intergenic.py --input %s.mRNA.intersect' %(prefix))

def update_key(gff):
    rank = 0
    gff_new = defaultdict(lambda : str())
    for i in sorted(gff.keys()):
        #print 'update key'
        #print rank, i
        if gff[i]:
            gff_new[rank] = gff[i]
            rank += 1
    return gff_new

def excision_frequency(samples):
    data = defaultdict(lambda : defaultdict(lambda : int()))
    for i in samples.keys():
        unit = re.split(r'\t', samples[i])
        data[unit[0]][unit[3]] = 1
    distance = defaultdict(lambda : int())
    for chro in data.keys():
        start_chro = sorted(data[chro].keys(), key=int)
        distance['%s_%s' %(chro, start_chro[0])] = abs(int(start_chro[0]) - int(start_chro[1]))
        for i in range(1, len(start_chro)-1):
            dist1 = abs(int(start_chro[i]) - int(start_chro[i-1]))
            dist2 = abs(int(start_chro[i+1]) - int(start_chro[i]))
            dist  = dist1 if dist1 < dist2 else dist2
            distance['%s_%s' %(chro, start_chro[i])] = dist
        distance['%s_%s' %(chro, start_chro[-1])] = abs(int(start_chro[-2]) - int(start_chro[-1]))
    frq = defaultdict(lambda : float())
    for index in sorted(distance.keys()): 
        #print 'index, distance: %s, %s' %(index, distance[index])
        if distance[index] < 10000:
            frq[index] = 0.002
        elif distance[index] < 100000:
            frq[index] = 0.0004
        else:
            frq[index] = 0.0001
    return frq

##random sample "num" insertion from gff and return the sampled list
def sample_mPing(gff, num):
    insertions = defaultdict(lambda : str())
    for n in range(num):
        rn = random.randint(0, len(gff.keys())-1)
        insertions[n] = gff[rn]
    #print 'Insertion: %s' %(insertions[0])
    return update_key(insertions)

##random sample "num" insertion from gff and return the sampled list
def insert_mPing(gff, somatic, num):
    insertions = gff
    index_base = len(gff.keys())
    for n in range(num):
        rn = random.randint(0, len(somatic.keys())-1)
        index = index_base + n
        #print 'insertion'
        #print  n, rn, index
        #print somatic[rn]
        insertions[index] = somatic[rn]
        #print 'Insertion: %s' %(insertions[0])
    return update_key(insertions)

##random remove "num" insertion from gff and return the sampled list
def excise_mPing_frq(gff, gff_raw, num, frq):
    insertions = gff_raw
    excision   = 0
    for n in range(len(gff_raw.keys())):
        rn = random.randint(0, 10000)
        unit = re.split(r'\t', gff_raw[n])
        index= '%s_%s' %(unit[0], unit[3])
        prob = frq[index] * 10000
        #print index
        #print 'random: %s, %s' %(rn, prob)
        if rn < prob:
            del insertions[n]
            #print n, rn
            #print 'Excision: %s' %(insertions[rn])
            excision += 1
    #print 'number of excision: %s' %(excision)
        #del insertions[rn]
    return update_key(insertions)



##random remove "num" insertion from gff and return the sampled list
def excise_mPing1(gff, gff_raw, num, frq):
    insertions = gff_raw
    for n in range(num):
        rn = random.randint(0, len(gff_raw.keys())-1)
        print n, rn
        print 'Excision: %s' %(insertions[rn])
        del insertions[rn]
    return update_key(insertions)

def evolve(samples, somatic, generation):
    gain = 40
    loss = 1
    for i in range(generation): 
        samples = insert_mPing(samples, somatic, gain)
        frq     = excision_frequency(samples)
        samples = excise_mPing_frq(samples, samples, loss, frq)
    return samples

def valid_sample(samples):
    n_sample = len(samples.keys())
    n_values = 0
    for k in samples.keys():
        #print '%s:%s' %(k, samples[k])
        if len(samples[k]) > 10:
            n_values += 1
    #print 'number of key: %s' %(n_sample)
    #print 'number of value: %s' %(n_values)

def simulate_excision(sample):
    sim_size   = 1000 #size of subsample
    sim_run    = 10  #number of run
    sim_generation = 1000
    sample_num = len(sample.keys())
    outdir = 'simulation_samplesize%s_numofrun%s' %(sim_size, sim_run)
    createdir(outdir)
    for r in range(sim_run):
        #sample a start sample from all somatic insertions
        samples = sample_mPing(sample, sim_size)
        #for n in range(sim_size):
        #    rn = random.randint(1, sample_num)
        #    print 'run%s\tsample%s\trn:%s' %(r, n, rn)
        #    samples.append(sample[int(rn)])
        samples_e = evolve(samples, sample, sim_generation)
        valid_sample(samples_e)
        prefix = 'Simulation%s' %('%04d' %(r))
        writefile('\n'.join(samples_e.values()), '%s.gff' %(prefix))
        distr_gff('%s.gff' %(prefix), prefix)
        os.system('mv %s.* %s' %(prefix, outdir))
    os.system('python Sim_Sum.py --input %s --output %s_results' %(outdir, outdir))
    distr_file='%s_results.mRNA.5primer.distance.distr' %(outdir)
    R_cmd='''
error.bar <- function(x, y, upper, lower=upper, color,length=0.06,...){
     if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
     stop("vectors must be same length")
     arrows(x,y+upper, x, y-lower, col=color,angle=90, code=3, length=length, ...)
 }

pdf("mping_intergenic_5distance_withsim.pdf")

par(mar=c(6,4,4,2), cex=1.2)
som5 <- read.table("random.mRNA.5primer.distance.distr")
#str5 <- read.table("../mPing_distr/Strains.mRNA.5primer.distance.distr")
#ril5 <- read.table("../mPing_distr/RIL.mRNA.5primer.distance.distr")
sim5 <- read.table("%s")

som5 <- som5[-1,]
#str5 <- str5[-1,]
#ril5 <- ril5[-1,]
sim5 <- sim5[-1,]

som5 <- som5[-length(som5[,1]),]
#str5 <- str5[-length(str5[,1]),]
#ril5 <- ril5[-length(ril5[,1]),]
sim5 <- sim5[-length(sim5[,1]),]

plot(rev(som5[,4]), type='b', pch= 1,lwd = 2 , col="aquamarine3", xaxt='n', frame.plot = FALSE, ylim=c(0,0.2), ylab="Proportion", xlab="")
#lines(rev(ril5[,4]), type='b',pch= 2,lwd = 2 , col="steelblue2")
#lines(rev(str5[,4]), type='b',pch= 3,lwd = 2 , col="sandybrown")
lines(rev(sim5[,4]), type='b',pch= 20, cex=0.2,lwd = 2 , col="dim gray")
error.bar(1:length(sim5[,4]), rev(sim5[,4]), rev(sim5[,7]-sim5[,4]), rev(sim5[,7]-sim5[,4]), 'dim gray')

#yaxis <- seq(1:length(som5[,1])+0.5
axis(1,seq(1:length(som5[,1])),line=0, labels=rep("",length(som5[,1])))
text(seq(1:length(som5[,1][-1]))+0.5,rep(-0.02,7), cex=1, offset=2,labels=rev(som5[,1]*500/-1000)[-1],srt=55,xpd=TRUE)

legend('topright', bty='n', border='NA', lty= c(1,2,3,4), pch = c(1,2,3,20), cex=1 , lwd = 2 ,col=c("aquamarine3", "steelblue2", "sandybrown", "dim gray"), c("Somatic", "RIL", "Strains", "Simulation"))
mtext("Distance to TSS (kp)", side=1,cex=1.2, at=9,line=3)

dev.off()
''' %(distr_file)
    writefile(R_cmd, 'mping_intergenic_5distance_withsim.R')
    os.system('cat mping_intergenic_5distance_withsim.R | R --slave')

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

    sample_start = readtable(args.input)
    simulate_excision(sample_start)
    distr_gff(args.input, 'random')
    

if __name__ == '__main__':
    main()

