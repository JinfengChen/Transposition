echo "prepare run"
perl prefastq_ln.pl --RIL /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/Illumina --fastq /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/fastq
perl prepare_run.pl --input RIL.fastq.list
echo "run"
bash scripts_11.sh > scripts_11.sh.log 2> scripts_11.sh.log2 &
bash scripts_12.sh > scripts_12.sh.log 2> scripts_12.sh.log2 &
bash scripts_13.sh > scripts_13.sh.log 2> scripts_13.sh.log2 &

