echo "test run"
python RunRelocaTEi_bam_Pop.py --input /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam > log 2> log2 &
python RunRelocaTEi_bam_Pop.py --input /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam_simulation > log 2> log2 &
echo "test with existingRun: subdata extracted"
python RunRelocaTEi_bam_Pop.py --input /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam --existingRun ./v2 > log 2> log2&
echo "ref as chr3 not whole genome, whole genome give less confident. but we still need to use whole genome"
python RunRelocaTEi_bam_Pop.py --input /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam_simulation --genome /rhome/cjinfeng/BigData/00.RD/RelocaTE_i/Simulation/Reference/MSU7.Chr3.fa
echo "test mping, ping, pong: when using similar element like ping and mping, some read will map to both, will get confused some place, less confident at breakpoint"
python RunRelocaTEi_bam_Pop.py --input /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam_simulation --repeat /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/TE/mping_all.fa > log 2> log2 &
echo "test Insert Size"
python RunInsertSize.py --input /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam &
echo "test temp"
python RunTEMP_Pop.py --input /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam_simulation > log 2> log2 &
python RunTEMP_Pop_qsub.py --input /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam_simulation > log 2> log2 &

echo "Run RelocaTEi on 275 RILs"
python RunInsertSize.py --input /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam --project RIL275 > log 2> log2 &
python RunRelocaTEi_bam_Pop.py --input /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam --project RIL275_RelocaTEi > log 2> log2 &
python RunCheckResults.py --input RIL275_RelocaTEi --tools RelocaTEi > RIL275_RelocaTEi.lowcpmping.list

echo "Run TEMP on 275 RILs"
python RunTEMP_Pop.py --input /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam --project RIL275_TEMP > log 2> log2 & 
python RunTEMP_Pop_qsub.py --input /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam --project RIL275_TEMP > log 2> log2 &
python RunCheckResults.py --input RIL275_TEMP --tools TEMP > RIL275_TEMP.lowcpmping.list
python RunTEMP_sumtable.py --input RIL275_TEMP

