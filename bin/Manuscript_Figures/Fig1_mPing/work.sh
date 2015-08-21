echo "default RelocaTE2"
paste RIL275_RelocaTEi.summary.table RIL275_RelocaTEi.CombinedGFF.characterized.mping.shared_unique_table.txt > RIL275_RelocaTEi.summary.unique.table
cat mPing_Call_Correlation_unique.R | R --slave

echo "clean all unsure TSD"
paste RIL275_RelocaTEi.summary_clean.table RIL275_RelocaTEi.CombinedGFF.characterized.clean.mping.shared_unique_table.txt > RIL275_RelocaTEi.summary_clean.unique.table
cat mPing_Call_Correlation_unique_clean.R | R --slave

cat mPing_class_proportion.R | R --slave
cat mPing_type_proportion.R | R --slave
cat ping_number_clean.R | R --slave
cat ping_number_clean_hom.R | R --slave
cat ping_number_clean_som.R | R --slave
