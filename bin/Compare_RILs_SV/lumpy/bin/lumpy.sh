#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=40gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`

echo "Do thing here"

~/BigData/software/SVcaller/lumpy-sv/bin/lumpyexpress \
     -B sample.1.fq.bam \
     -S sample.1.fq.splitters.bam \
     -D sample.1.fq.discordants.bam \
     -o sample.1.fq.vcf

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

