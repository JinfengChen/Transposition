echo "Prepare excision number for RILs"
python Excision_Number_In_RILs.py --input RIL230.sample.list

echo "Prepare mPing, Ping, excision table"
python Excision_Number_Ping.py --excision Excision_newpipe_version1.footprint.list.noPing.rils.txt --ping RIL230_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.narrow_range.txt
python Excision_Number_Ping.py --excision Excision_newpipe_version1.footprint.list.noPing.rils.txt --ping RIL230_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.ping_code.txt --output RIL230.RIL_mPing_Ping_Excision.table.txt 


echo "correlation and plot"
cat ping_number_clean_excision_cor_plot.R | R --slave

