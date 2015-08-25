echo "unique ping, all confident calls. cleaned all unsure TSD, use perfect match for overlap"
python Clean_Calls.py --gff RIL275_RelocaTEi.CombinedGFF.characterized.gff
python Unique_mPing_clean.py --input RIL275_RelocaTEi.CombinedGFF.characterized.clean.gff
cat ping_number.R | R --slave
python Sum_unique_clean.py --input RIL275_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff > RIL275_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.summary
python Sum_class_clean.py --input RIL275_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff
python Sum_Ping_mPing.py
echo "Check"
awk '$2!=$11' RIL275_RelocaTEi.CombinedGFF.characterized.clean.overlap_ril | cut -f1,2,3,4,5,6 | grep "RIL1_0" | uniq | wc -l
grep "RIL1_0" RIL275_RelocaTEi.CombinedGFF.characterized.clean.overlap_ref -c

#echo "unique ping, all calls"
#python Unique_mPing.py --input RIL275_RelocaTEi.CombinedGFF.ALL.gff
#cat ping_number.ALL.R | R --slave
#python Sum_unique.py --input RIL275_RelocaTEi.CombinedGFF.ALL.unique_mPing.gff > RIL275_RelocaTEi.CombinedGFF.ALL.unique_mPing.summary

echo "annotate"
perl AnnoTEsite.pl --te RIL275_RelocaTEi.CombinedGFF.characterized.unique_mPing.gff --gene /rhome/cjinfeng/BigData/00.RD/seqlib/GFF/MSU7.gene.gff &
