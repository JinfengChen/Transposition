#new excison pipe
#call excision
#/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/Sum_excision_distance.py
#footprint
#/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Compare_footprint_events/bin/footprint_events.py

#mping.excision.avg_frequency.py
#get mping.excision.avg.table use Excision_newpipe_version1.footprint.list.noPing.txt 
#and RIL230_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.parental.frequency 
python mping.excision.avg_frequency.py

echo "parental mPing allele frequency"
grep "Parental" RIL275_RelocaTEi.CombinedGFF.characterized.clean.mping.ril.frequency > RIL275_RelocaTEi.CombinedGFF.characterized.clean.mping.ril.parental.frequency
grep "Parental" RIL230_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.frequency > RIL230_RelocaTEi.CombinedGFF.characterized.clean.shared_mping.ril.parental.frequency

echo "nonparental mPing alleley frequency"
grep -v "Parental" RIL275_RelocaTEi.CombinedGFF.characterized.clean.mping.ril.frequency > RIL275_RelocaTEi.CombinedGFF.characterized.clean.mping.ril.nonparental.frequency

echo "plot"
cat mping.allele_frq.ALL_shared.R | R --slave
cat mping.allele_frq.parental.R | R --slave
cat mping.allele_frq.nonparental.R | R --slave
cat mping.excision.avg.R | R --slave
