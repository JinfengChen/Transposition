samtools view ~/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam/GN155.bam Chr8:24625200-24625350 | less -S
blat -minScore=10 -tileSize=7 ~/BigData/00.RD/seqlib/MSU_r7.fa HP13_R155.unmapped.fa HP13_R155.unmapped.blat
