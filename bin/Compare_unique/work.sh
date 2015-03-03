echo "unique ping, all confident calls"
python Unique_mPing.py --input RIL275_RelocaTEi.CombinedGFF.characterized.gff
cat ping_number.R | R --slave

echo "unique ping, all calls"
python Unique_mPing.py --input RIL275_RelocaTEi.CombinedGFF.ALL.gff
cat ping_number.ALL.R | R --slave
