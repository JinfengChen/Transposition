echo "unique ping, all confident calls"
python Unique_mPing.py --input RIL275_RelocaTEi.CombinedGFF.characterized.gff
cat ping_number.R | R --slave
python Sum_unique.py --input RIL275_RelocaTEi.CombinedGFF.characterized.unique_mPing.gff > RIL275_RelocaTEi.CombinedGFF.characterized.unique_mPing.summary

echo "unique ping, all calls"
python Unique_mPing.py --input RIL275_RelocaTEi.CombinedGFF.ALL.gff
cat ping_number.ALL.R | R --slave
python Sum_unique.py --input RIL275_RelocaTEi.CombinedGFF.ALL.unique_mPing.gff > RIL275_RelocaTEi.CombinedGFF.ALL.unique_mPing.summary

echo "annotate"
perl AnnoTEsite.pl --te RIL275_RelocaTEi.CombinedGFF.characterized.unique_mPing.gff --gene /rhome/cjinfeng/BigData/00.RD/seqlib/GFF/MSU7.gene.gff &
