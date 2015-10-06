samtools view ~/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam/GN26.bam Chr3:6,441,398-6,443,105 | less -S
blat -minScore=10 -tileSize=7 ~/BigData/00.RD/seqlib/MSU_r7.fa HP7_R26.unmapped.fa HP7_R26.unmapped.blat
