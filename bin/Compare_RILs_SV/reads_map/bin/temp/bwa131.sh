#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l mem=30gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

#cd $PBS_O_WORKDIR


start=`date +%s`
fq1=../input/RILs_ALL_unmapped_mping_fastq/RIL131/RIL131_1.fq
fq2=../input/RILs_ALL_unmapped_mping_fastq/RIL131/RIL131_2.fq
ref=../input/bwa_0.7/MSU_r7.Pseudo_mPing.fa
#perl step1_Mapping.pl -ref $ref -1 $fq1 -2 $fq2 -min 0 --max 500 -cpu $PBS_NP --tool bwa --project RIL131_test
#/opt/linux/centos/7.x/x86_64/pkgs/bwa/0.7.12/bin/bwa aln -t $PBS_NP $ref $fq1 > RIL131.1.sai
#/opt/linux/centos/7.x/x86_64/pkgs/bwa/0.7.12/bin/bwa aln -t $PBS_NP $ref $fq2 > RIL131.2.sai
#/opt/linux/centos/7.x/x86_64/pkgs/bwa/0.7.12/bin/bwa sampe -a 500 -A $ref RIL131.1.sai RIL131.2.sai $fq1 $fq2 > RIL131_test.sam

prefix=RIL131
/opt/linux/centos/7.x/x86_64/pkgs/samtools/1.2/bin/samtools view -@ $PBS_NP -bS -o $prefix.raw.bam $prefix.sam
/opt/linux/centos/7.x/x86_64/pkgs/samtools/1.2/bin/samtools sort -@ $PBS_NP $prefix.raw.bam $prefix.sort
java -Xmx10G -jar /opt/picard/1.81/MarkDuplicates.jar ASSUME_SORTED=TRUE REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT INPUT=$prefix.sort.bam OUTPUT=$prefix.bam METRICS_FILE=$prefix.dupli


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
