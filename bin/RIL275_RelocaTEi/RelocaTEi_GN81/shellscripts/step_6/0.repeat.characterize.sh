cat /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN81/repeat/results/*.all_nonref_insert.gff > /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN81/repeat/results/ALL.all_nonref_insert.gff
cat /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN81/repeat/results/*.all_nonref_insert.txt | grep "^TE" -v > /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN81/repeat/results/ALL.all_nonref_insert.txt
perl /rhome/cjinfeng/software/tools/RelocaTE_1.0.3_i/RelocaTE/scripts/characterizer.pl -s /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTEi_GN81/repeat/results/ALL.all_nonref_insert.txt -b /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam/GN81.bam -g /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa --samtools /opt/samtools-0.1.16/samtools
