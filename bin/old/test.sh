perl /rhome/cjinfeng/software/tools/RelocaTE-master/scripts/relocaTE.pl -t ../input/TE/Rice.TE.fa -g /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa -d ../input/fq/GN100 -e GN100 -o ./RelocaTE/GN100 -r ./all_reference_repeat.txt -p 1 -a 0
#python runRelocaTEjobs.py --input ./RelocaTE/GN100/run_these_jobs.sh
#perl /rhome/cjinfeng/software/tools/RelocaTE-master/scripts/characterizer.pl -s ./RelocaTE/GN100/mping/results/GN100.mping.all_nonref.txt -b /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam/GN100.bam -g /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa -x 1
#perl /rhome/cjinfeng/software/tools/RelocaTE-master/scripts/relocaTE.pl -t ../input/TE/mPing_Ping_Pong.fa -g /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa -d ../input/fq/GN101 -e GN101 -o ./RelocaTE/GN101 -r 1 -p 1 -a 0
#python runRelocaTEjobs.py --input ./RelocaTE/GN101/run_these_jobs.sh
