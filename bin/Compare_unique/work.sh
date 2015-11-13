echo "unique ping, all confident calls. cleaned all unsure TSD, use perfect match for overlap"
python Clean_Calls.py --gff RIL275_RelocaTEi.CombinedGFF.characterized.gff

echo "make unique mPing gff, generate unique mPing number for ping number and single ping"
python Unique_mPing_clean.py --input RIL275_RelocaTEi.CombinedGFF.characterized.clean.gff
cat ping_number.R | R --slave

echo "get list of unique mPing number for each RILs, generate mping frquency for shared mPing in RILs, can check parental mPing, shared in RILs"
python Sum_unique_clean.py --input RIL275_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff > RIL275_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.summary
python MergePingCode.py --input RIL230_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.txt

echo "classify mPing into parental, shared in ril, or unique and then classify into hom, het, som"
python Sum_class_clean.py --input RIL275_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff

echo "different from Unique_mPing_clean.py. have het and som mPing number but from summary table which is roughly right but not clean"
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


echo "update to 230 core"
#remove 6 ping from HEG4.ALL.mping.non-ref.AF0.1.gff, in excision we used ping but can not use here. 
python Clean_Calls.py --gff RIL230_RelocaTEi.CombinedGFF.characterized.gff
python Unique_mPing_clean.py --input RIL230_RelocaTEi.CombinedGFF.characterized.clean.gff --reference HEG4.ALL.mping.non-ref.AF0.1.gff --code RIL230_RelocaTE.sofia.ping_code.table
python Sum_unique_clean.py --input RIL230_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff > RIL230_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.summary
python MergePingCode.py --input RIL230_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.txt
echo "need to rewrite this when classify hom/het/som, use large proportion of genotype to represent. this will minimum the error of het by method"
python Sum_class_clean.py --input RIL230_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff
python Sum_Ping_mPing.py --code RIL230_RelocaTE.sofia.ping_code.table --output RIL230_RelocaTEi.CombinedGFF.characterized.clean
python Sum_Ping_mPing_Narrow_range.py --code RIL230_RelocaTE.sofia.ping_code.table --output RIL230_RelocaTEi.CombinedGFF.characterized.clean.narrow_range
awk '$2>=163 && $2<=203' RIL230_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt > RIL230_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.narrow_range.txt
python Sum_Ping_mPing_High_depth.py --code RIL230_RelocaTE.sofia.ping_code.table --output RIL230_RelocaTEi.CombinedGFF.characterized.clean.high_depth
python Sum_Ping_mPing_High_depth_table.py --table RIL230_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt --output high_depth
python Sum_Ping_mPing_High_depth.py --code RIL230_RelocaTE.sofia.ping_code.table --output RIL230_RelocaTEi.CombinedGFF.characterized.clean.high_narrow --narrow
python Sum_Ping_mPing_High_depth_table.py --table RIL230_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt --output high_narrow --narrow

