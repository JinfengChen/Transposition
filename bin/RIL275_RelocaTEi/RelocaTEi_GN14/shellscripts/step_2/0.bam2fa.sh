/opt/samtools-0.1.16/samtools view -h /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam/GN14.bam | awk '$5<60' | samtools view -Shb - | samtools sort -m 500000000 -n - /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN14/repeat/fastq/GN14.subset 2> /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN14/run.std
/opt/bedtools/2.17.0-25-g7b42b3b/bin//bedtools bamtofastq -i /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN14/repeat/fastq/GN14.subset.bam -fq /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN14/repeat/fastq/GN14_1.fq -fq2 /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN14/repeat/fastq/GN14_2.fq 2> /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN14/run.std
/rhome/cjinfeng/software/tools/RelocaTE_1.0.3_i/RelocaTE/scripts/relocaTE_fq2fa.pl /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN14/repeat/fastq/GN14_1.fq /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN14/repeat/fastq/GN14_1.fa
/rhome/cjinfeng/software/tools/RelocaTE_1.0.3_i/RelocaTE/scripts/relocaTE_fq2fa.pl /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN14/repeat/fastq/GN14_2.fq /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN14/repeat/fastq/GN14_2.fa
