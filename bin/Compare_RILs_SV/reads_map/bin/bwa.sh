#!/bin/bash
#PBS -l nodes=1:ppn=24
#PBS -l mem=60gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`

perl step1_Mapping.pl -ref ../input/MSU_r7.Pseudo_mPing.10kb_flank.fa -1 ../input/RILs_ALL_fastq/RIL1/RIL1_1.fq -2 ../input/RILs_ALL_fastq/RIL1/RIL1_2.fq -min 0 --max 500 -cpu 24 --tool bwa --project test

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

