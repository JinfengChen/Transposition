#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=10:00:00
#PBS -d ./

#cd $PBS_O_WORKDIR


for i in `ls -d RelocaTEi_* | sed 's/@//'`
do
   echo $i
   cd $i
   cd repeat
   rm -R blat_output
   rm -R flanking_seq
   rm -R te_containing_fq/
   rm -R te_only_read_portions_fa/
   rm bwa_aln/*.mates.bam*
   rm bwa_aln/*.unPaired.bam*
   rm bwa_aln/*.bwa.bam*
   rm -R fastq_split/
   rm -R fastq
   cd ../..
done

