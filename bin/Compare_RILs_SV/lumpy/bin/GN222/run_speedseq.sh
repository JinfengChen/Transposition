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
RIL=GN222

# Example speedseq commands on a small slice of chromosome 20

# 1. Align with BWA
#/bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq/bin/speedseq align -t 16 -o $RIL -R "@RG\tID:id\tSM:$RIL\tLB:lib"\
#   MSU_r7.fa \
#   $RIL\.1.fq \
#   $RIL\.2.fq

/bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq/bin/speedseq align \
     -t $PBS_NP \
     -o $RIL \
     -R "@RG\tID:id\tSM:$RIL\tLB:lib" \
     /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/lumpy/bin/MSU_r7.fa \
     $RIL\.1.fq \
     $RIL\.2.fq

# 2. Detect SNVs and indels
#/bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq/bin/speedseq var \
#    -o example \
#    data/human_g1k_v37_20_42220611-42542245.fasta \
#    example.bam

# 3. Detect SVs
/bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq/bin/speedseq sv \
    -o $RIL \
    -B $RIL\.bam \
    -S $RIL\.splitters.bam \
    -D $RIL\.discordants.bam \
    -R /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/lumpy/bin/MSU_r7.fa \
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
