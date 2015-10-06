#!/bin/bash
#PBS -l nodes=1:ppn=24
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
#/bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq/bin/speedseq align \
#     -t $PBS_NP \
#     -o sample \
#     -R "@RG\tID:id\tSM:RIL26\tLB:lib" \
#     -M 3 \
#     MSU_r7.fa \
#     sample.1.fq \
#     sample.2.fq

# 2. Detect SNVs and indels
#/bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq/bin/speedseq var \
#    -o example \
#    data/human_g1k_v37_20_42220611-42542245.fasta \
#    example.bam

# 3. Detect SVs
/bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq/bin/speedseq sv \
    -o sample \
    -B sample.bam \
    -S sample.splitters.bam \
    -D sample.discordants.bam \
    -R MSU_r7.fa \
    -t $PBS_NP \
    -d \
    -k \
    -v

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
