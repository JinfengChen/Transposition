samtools view ~/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam/GN34.bam Chr2:30,934,039-30,934,880 | less -S
blat -minScore=10 -tileSize=7 ~/BigData/00.RD/seqlib/MSU_r7.fa HP6_R34.unmapped.fa HP6_R34.unmapped.blat

