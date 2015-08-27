echo "Number of Parental mPing in each RILs"
python Parental_mPing_RILs.py --input ../input/RIL275_RelocaTEi.CombinedGFF.characterized.gff --reference ../input/HEG4.ALL.mping.non-ref.AF0.1.gff --project RIL275.Parental_mPing
python Parental_mPing_RILs.py --input ../input/RIL275_RelocaTEi.CombinedGFF.characterized.gff --reference ../input/HEG4.ALL.mping.non-ref.AF0.1.linked_100kb.gff --project RIL275.Parental_mPing_linked_100kb

echo "add ping number and excision number to mping summary table"
python Excision_transposases.py > Excision_transposases.table.txt
echo "summary excision according to ping number"
python Excision_Ping_number.py
echo "Ping vs. Excision and Ping, mPing vs. Excision"
cat ping_number_excision.R | R --slave
cat ping_mping_number_excision.R | R --slave

echo "Ping code in RILs and high excision mPing"
python Ping_number_RILs.py --ping_code ../input/RIL275_RelocaTE.sofia.ping_code.table --high_excision mping.excision.draw.highexcision --exicision ../input/mping.excision.ril.counts

#for each mPing, add excision code and ping code to the end of csv matrix
python Ping_number_RILs.High_exicison.py --csv ../../Compare_RILs_SV/reads_map/bin/mPing_boundary_mPing --ping_code ../input/RIL275_RelocaTE.sofia.ping_code.table

