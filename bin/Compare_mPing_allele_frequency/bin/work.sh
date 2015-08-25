python mPing_allele_frequency.py --frequency ../input/mping.ril.frquency --gff ../input/HEG4.ALL.mping.non-ref.gff > HEG4.ALL.mping.non-ref.allele.frq
#get subgff of allele frequency more than 0.1
python mPing_allele_frequency.py --frequency ../input/mping.ril.frquency --gff ../input/HEG4.ALL.mping.non-ref.gff --subgff HEG4.ALL.mping.non-ref.AF0.1.gff > HEG4.ALL.mping.non-ref.allele.frq
grep "NA" -v HEG4.ALL.mping.non-ref.allele.frq | awk '$2>0.1' | less -S
#552 test from HEG4 mPing/Ping
#446 present in RILs
#419 present with frequency more than 10%
