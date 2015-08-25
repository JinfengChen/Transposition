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
python RunRelocaTEi_sumtable.py --input RIL275_RelocaTEi
python RunRelocaTEi_sumtable_clean.py --input RIL275_RelocaTEi
python RunRelocaTEi_CombinedGFF.py --input RIL275_RelocaTEi


echo "Run TEMP on 275 RILs"
python RunTEMP_Pop.py --input /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam --project RIL275_TEMP > log 2> log2 & 
python RunTEMP_Pop_qsub.py --input /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/bam --project RIL275_TEMP > log 2> log2 &
python RunCheckResults.py --input RIL275_TEMP --tools TEMP > RIL275_TEMP.lowcpmping.list
python RunTEMP_sumtable.py --input RIL275_TEMP


echo "Run RelocaTE on 275 RILs"
cut -f2,3 EG4.mping.all_reference.txt > Nipponbare.mPing.txt
python RunRelocaTE_Pop.py --input /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_fastq --project RIL275_RelocaTE
python RunCheckResults.py --input RIL275_RelocaTE --tools RelocaTE > RIL275_RelocaTE.lowcpmping.list

python RunRelocaTE_CombinedGFF.py --input RIL275_RelocaTE

echo "add unique mping"
python MergeTable.py --table1 Compare_fig/RIL275_RelocaTEi.summary.table --table2 RIL275_RelocaTE.sofia.unique.table2 > Compare_fig/RIL275_RelocaTEi_unique.summary.table

