echo "parental mPing allele frequency"
grep "Parental" RIL275_RelocaTEi.CombinedGFF.characterized.clean.mping.ril.frequency > RIL275_RelocaTEi.CombinedGFF.characterized.clean.mping.ril.parental.frequency

echo "nonparental mPing alleley frequency"
grep -v "Parental" RIL275_RelocaTEi.CombinedGFF.characterized.clean.mping.ril.frequency > RIL275_RelocaTEi.CombinedGFF.characterized.clean.mping.ril.nonparental.frequency

echo "plot"
cat mping.allele_frq.ALL_shared.R | R --slave
cat mping.allele_frq.parental.R | R --slave
cat mping.allele_frq.nonparental.R | R --slave
cat mping.excision.avg.R | R --slave
