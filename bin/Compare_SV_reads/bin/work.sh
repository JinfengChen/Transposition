echo "Chr2_28724474"
~/BigData/software/seqtk-master/seqtk subseq ~/BigData/00.RD/RILs/QTL_pipe/input/fastq/Illumina_fixed_link/RIL11_0/RIL11_0_ACTGAT_FC153L5_p1.fq Chr2_28724475_RIL11.reads > Chr2_28724475_RIL11.reads.1.fq
~/BigData/software/seqtk-master/seqtk subseq ~/BigData/00.RD/RILs/QTL_pipe/input/fastq/Illumina_fixed_link/RIL11_0/RIL11_0_ACTGAT_FC153L5_p2.fq Chr2_28724475_RIL11.reads > Chr2_28724475_RIL11.reads.2.fq
blat ../input/MSU_r7.fa Chr2_28724475_RIL11.reads.new_seq.fa Chr2_28724475_RIL11.reads.blat
blat ../input/MSU_r7.fa Chr2_28724475_RIL11.sanger.fa Chr2_28724475_RIL11.sanger.blat

echo "Chr5_15210100"
samtools view ../../../../QTL_pipe/input/fastq/RILs_ALL_bam/GN222.bam Chr5:15,210,100-15,210,448 | less -S
~/BigData/software/seqtk-master/seqtk subseq ~/BigData/00.RD/RILs/QTL_pipe/input/fastq/Illumina_fixed_link/RIL222_0/RIL222_0_GTGGCC_FC197L6_p1.fq Chr5_15210100_RIL222.reads > Chr5_15210100_RIL222.reads.1.fq
~/BigData/software/seqtk-master/seqtk subseq ~/BigData/00.RD/RILs/QTL_pipe/input/fastq/Illumina_fixed_link/RIL222_0/RIL222_0_GTGGCC_FC197L6_p2.fq Chr5_15210100_RIL222.reads > Chr5_15210100_RIL222.reads.2.fq
blat -minScore=10 -tileSize=7 ../input/MSU_r7.fa Chr5_15210100_RIL222.reads.new_seq.fa Chr5_15210100_RIL222.reads.blat

echo "Chr8_24625306"
samtools view ../../../../QTL_pipe/input/fastq/RILs_ALL_bam/GN22.bam Chr8:24,634,481-24,634,811 | less -S
blat -minScore=10 -tileSize=7 ../input/MSU_r7.fa Chr8_24625306_RIL22.reads.new_seq.fa Chr8_24625306_RIL22.reads.blat
blat -minScore=10 -tileSize=7 ../input/MSU_r7.fa Chr8_24625306_RIL133.sanger.fa Chr8_24625306_RIL133.sanger.blat
blastall -p blastn -i Chr8_24625306_RIL133.sanger.fa -d ../input/MSU_r7.fa -o Chr8_24625306_RIL133.sanger.blast -e 1e-5

echo "Chr3_6441699"
blat -minScore=10 -tileSize=7 ../input/MSU_r7.fa Chr3_6441699_RIL31.sanger.fa Chr3_6441699_RIL31.sanger.blat
samtools view ../../../../QTL_pipe/input/fastq/RILs_ALL_bam/GN31.bam Chr3:6,441,661-6,441,766 | less -S
samtools view ../../../../QTL_pipe/input/fastq/RILs_ALL_bam/GN31.bam Chr3:6,442,461-6,442,513 | less -S 

#R26
blat -minScore=10 -tileSize=7 ../input/MSU_r7.fa Chr3_6441699_RIL26.sanger_right_junction.fa Chr3_6441699_RIL26.sanger_right_junction.blat

echo "Chr8.24625267, RIL155"
samtools view ../../../../QTL_pipe/input/fastq/RILs_ALL_bam/GN155.bam Chr8:24,625,267-24,625,441
blat -minScore=10 -tileSize=7 ../input/MSU_r7.fa Chr8_24625267_RIL155.reads.1.fa Chr8_24625267_RIL155.reads.1.blat

echo "Chr8.24625267, RIL22"
blat -minScore=10 -tileSize=7 ../input/MSU_r7.fa Chr8_24625306_RIL22.sanger.fa Chr8_24625306_RIL22.sanger.blat
blat -minScore=10 -tileSize=7 ../input/MSU_r7.fa Chr8_24625306_RIL22.sanger_left_junction.fa Chr8_24625306_RIL22.sanger_left_junction.blat

