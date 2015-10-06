#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l mem=60gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR

module load samtools
PATH=$PATH:~/BigData/software/SVcaller/ROOT/bin/
start=`date +%s`

# Example speedseq commands on a small slice of chromosome 20

# 1. Align with BWA
#../bin/speedseq align \
#    -o example \
#    -M 3 \
#    -p \
#    -R "@RG\tID:NA12878\tSM:NA12878\tLB:lib1" \
#    data/human_g1k_v37_20_42220611-42542245.fasta \
#    data/NA12878.20slice.30X.fastq.gz

# 2. Detect SNVs and indels
#../bin/speedseq var \
#    -o example \
#    data/human_g1k_v37_20_42220611-42542245.fasta \
#    example.bam

# 3. Detect SVs
#/bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq/bin/speedseq sv \
#    -o example \
#    -B sample.1.fq.bam \
#    -S sample.1.fq.splitters.bam \
#    -D sample.1.fq.discordants.bam \
#    -R MSU_r7.fa \
#    -d \
#    -k \
#    -v
workdir=/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/lumpy/bin

/usr/bin/python2.7 /bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq//bin/cnvnator_wrapper.py --cnvnator /bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq//bin/cnvnator-multi -T cnvnator-temp -t 12 -w 100 -b sample.1.fq.bam -o sample.1.fq.bam.readdepth -c /bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq//annotations/cnvnator_chroms

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
