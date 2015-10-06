#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l mem=16gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`
fq1=../input/RILs_ALL_unmapped_mping_fastq/RIL1/RIL1_1.fq
fq2=../input/RILs_ALL_unmapped_mping_fastq/RIL1/RIL1_2.fq
ref=../input/MSU_r7.Pseudo_mPing.fa
perl step1_Mapping.pl -ref $ref -1 $fq1 -2 $fq2 -min 0 --max 500 -cpu $PBS_NP --tool bwa --project test1

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

