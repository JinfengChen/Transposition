blat -minScore=10 -tileSize=7 /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL91/mping/mping.fa /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_fastq/RIL91/RIL91_2.fa /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL91/mping/blat_output/RIL91_2.te_mping.blatout 1>> /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL91/mping/blat_output/blat.out
perl /rhome/cjinfeng/software/tools/RelocaTE_1.0.3/RelocaTE/scripts/relocaTE_trim.pl /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL91/mping/blat_output/RIL91_2.te_mping.blatout /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_fastq/RIL91/RIL91_2.fq 10 0 > /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL91/mping/flanking_seq/RIL91_2.te_mping.flankingReads.fq
