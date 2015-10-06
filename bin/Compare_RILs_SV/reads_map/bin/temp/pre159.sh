#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -V
#PBS -d ./

prefix=159
mkdir RIL$prefix

#/opt/linux/centos/7.x/x86_64/pkgs/samtools/1.2/bin/samtools view -hb -f 4 -F264 /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/input/RILs_ALL_bam/GN$prefix\.bam > /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.unmapped1.bam
#/opt/linux/centos/7.x/x86_64/pkgs/samtools/1.2/bin/samtools view -hb -f 8 -F260 /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/input/RILs_ALL_bam/GN$prefix\.bam > /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.unmapped2.bam
#/opt/linux/centos/7.x/x86_64/pkgs/samtools/1.2/bin/samtools view -hb -f 12 -F256  /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/input/RILs_ALL_bam/GN$prefix\.bam > /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.unmapped3.bam
#/opt/linux/centos/7.x/x86_64/pkgs/samtools/1.2/bin/samtools merge /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.unmapped.bam /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.unmapped[123].bam
#/rhome/cjinfeng/BigData/software/bam2fastq/bam2fastq-1.1.0/bam2fastq /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.unmapped.bam -o /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.unmapped#.fq
#/opt/linux/centos/7.x/x86_64/pkgs/samtools/1.2/bin/samtools view -hb -L /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/input/Parent.ALL.mPing.100kb_flank.merge.table /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/input/RILs_ALL_bam/GN$prefix\.bam > /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.mping.bam
#/opt/linux/centos/7.x/x86_64/pkgs/samtools/1.2/bin/samtools view -hb -f 3 /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.mping.bam > /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.mping.mapped.bam
#/rhome/cjinfeng/BigData/software/bam2fastq/bam2fastq-1.1.0/bam2fastq /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.mping.mapped.bam -o /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.mping.mapped#.fq
#cat /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.unmapped_1.fq /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.mping.mapped_1.fq > /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\_1.fq
#cat /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.unmapped_2.fq /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.mping.mapped_2.fq > /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\_2.fq
rm /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.unmapped* /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix/RIL$prefix\.mping*
mv /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL$prefix ../input/RILs_ALL_unmapped_mping_fastq3

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

