echo "sliding windows RIL unique mPing"
perl splitGFF.pl ../input/RIL275_RelocaTEi.CombinedGFF.characterized.unique_mPing.gff > log 2>&1 &
perl gff2windowsbar4TE.pl --window 200000 --step 50000 --chrlen ../input/MSU7.chrlen --bar unique_mping_bar --gff ../input/RIL275_RelocaTEi.CombinedGFF.characterized.unique_mPing.gff.chr
cat unique_mping_bar/Chr* > unique_mping.density.txt

echo "sliding windows RIL excision mPing"
perl splitGFF.pl ../input/mping.excision.non_ref.gff
perl gff2windowsbar4TE.pl --window 200000 --step 50000 --chrlen ../input/MSU7.chrlen --bar excision_mping_bar --gff ../input/mping.excision.non_ref.gff.chr
cat excision_mping_bar/Chr* > excision_mping.density.txt

echo "sliding windows bedtools: results same with perl script"
python mPing_slidingWin.py --gff ../input/RIL275_RelocaTEi.CombinedGFF.characterized.unique_mPing.gff --type unique > log 2>&1 &
python mPing_slidingWin.py --gff ../input/mping.excision.non_ref.gff --type excision > log 2>&1 &
python mPing_slidingWin.py --gff ../input/HEG4.ALL.mping.non-ref.gff --type HEG4 > log 2>&1 &

echo "correlation"
paste MSU7.SlidingWin.HEG4.density MSU7.SlidingWin.excision.density | awk '$7>0' > MSU7.SlidingWin.HEG4vsExcision.table

