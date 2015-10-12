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

# Example speedseq commands on a small slice of chromosome 20

# 1. Align with BWA
#../bin/speedseq align \
#    -o example \
#    -M 3 \
#    -p \
#    -R "@RG\tID:NA12878\tSM:NA12878\tLB:lib1" \
#    data/human_g1k_v37_20_42220611-42542245.fasta \
#    data/NA12878.20slice.30X.fastq.gz

# 2. Detect SNVs and indels
#../bin/speedseq var \
#    -o example \
#    data/human_g1k_v37_20_42220611-42542245.fasta \
#    example.bam

# 3. Detect SVs
#/bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq/bin/speedseq sv \
#    -o example \
#    -B sample.1.fq.bam \
#    -S sample.1.fq.splitters.bam \
#    -D sample.1.fq.discordants.bam \
#    -R MSU_r7.fa \
#    -d \
#    -k \
#    -v
workdir=/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/lumpy/bin
#
#/usr/bin/python2.7 /bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq//bin/cnvnator_wrapper.py --cnvnator /bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq//bin/cnvnator-multi -T $workdir/example/cnvnator-temp -t 12 -w 100 -b $workdir/sample.bam -o $workdir/example/sample.bam.readdepth -c /bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq//annotations/cnvnator_chroms

#>>>EXTRACTING READ MAPPING FROM BAM/SAM FILES
echo "Step1. EXTRACTING READ MAPPING FROM BAM/SAM FILES"
cnvnator=/rhome/cjinfeng/BigData/software/SVcaller/speedseq/src/cnvnator/bin
#temp=/rhome/cjinfeng/BigData/software/SVcaller/speedseq/example/cnvnator-temp
#bam=/rhome/cjinfeng/BigData/software/SVcaller/speedseq/example/sample.bam
temp=$workdir/cnvnator-temp
bam=$workdir/sample.1.fq.bam
bin=1000
chr=Chr10
$cnvnator/cnvnator-multi -root $temp/MSU7.root -tree $bam -unique
echo "Step1 Done"
echo ""

#>>>GENERATING HISTOGRAM
echo "Step2. GENERATING HISTOGRAM"
chrs=/bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq//annotations/cnvnator_chroms
$cnvnator/cnvnator-multi -root $temp/MSU7.root -outroot $temp/MSU7.hist.root -his $bin -d $chrs
echo "Step2 Done"
echo ""

#>>>CALCULATING STATISTICS
echo "Step3. CALCULATING STATISTICS"
$cnvnator/cnvnator-multi -root $temp/MSU7.hist.root -stat $bin
echo "Step3. Done"
echo ""

#>>>RD SIGNAL PARTITIONING
echo "Step4. RD SIGNAL PARTITIONING"
$cnvnator/cnvnator-multi -root $temp/MSU7.hist.root -partition $bin -ngc
echo "Step4 Done"
echo ""

#>>>CNV CALLING
echo "Step5. CNV CALLING"
$cnvnator/cnvnator-multi -root $temp/MSU7.hist.root -call $bin -ngc > $bam\.readdepth.txt
echo "Step5 Done"
echo ""


end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"
