echo "map"
qsub -q highmem map.sh

echo "SV"
qsub -q highmem run_speedseq.sh

echo "convert chr name"
samtools view -h sample.1.fq.bam | sed 's/Chr//g' | samtools view -Sb - -o sample.bam
samtools view -h sample.1.fq.discordants.bam | sed 's/Chr//g' | samtools view -Sb - -o sample.discordants.bam
samtools view -h sample.1.fq.splitters.bam | sed 's/Chr//g' | samtools view -Sb - -o sample.splitters.bam

echo "rsw_seq"
#nb chr3
bedtools bamtobed -i test.bam | cut -f2 > test.bed
#r26 chr3
bedtools bamtobed -i chr3.r26.bam | cut -f2 > chr3.r26.bed
#rsw_seq
~/BigData/software/SVcaller/rSW_seq/a.out chr3.r26.bed test.bed 2803040 3617875 50 > test.out
