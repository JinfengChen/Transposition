/usr/local/bin/blat -minScore=10 -tileSize=7 /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/TE/mping.fa /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN4/repeat/fastq/GN4_1.fa /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN4/repeat/blat_output/GN4_1.te_repeat.blatout 1>>/bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN4/repeat/blat_output/blat.out 2>>/bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN4/repeat/blat_output/blat.out
python /rhome/cjinfeng/software/tools/RelocaTE_1.0.3_i/RelocaTE/scripts/relocaTE_trim.py /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN4/repeat/blat_output/GN4_1.te_repeat.blatout /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN4/repeat/fastq/GN4_1.fq 10 1 > /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN4/repeat/flanking_seq/GN4_1.te_repeat.flankingReads.fq
