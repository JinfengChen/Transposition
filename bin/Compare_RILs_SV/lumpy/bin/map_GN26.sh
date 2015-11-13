#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l mem=40gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`

RIL=GN26
echo "Do thing here"
~/BigData/software/SVcaller/speedseq/bin/speedseq align -t 16 -o $RIL -R "@RG\tID:id\tSM:$RIL\tLB:lib"\
   MSU_r7.fa \
   $RIL\.1.fq \
   $RIL\.2.fq

#samtools sort sample.1.fq.discordants.unsorted.bam sample.discordants
#samtools sort sample.1.fq.splitters.unsorted.bam sample.splitters

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

