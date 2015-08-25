echo "mping between AF 0.1-0.3"
awk '$7>0.1 && $7<0.3' RIL275_RelocaTEi.CombinedGFF.characterized.clean.mping.ril.parental.frequency | sort -k1,1 -k2,2n
awk '$7>0.1 && $7<0.3' RIL275_RelocaTEi.CombinedGFF.characterized.clean.mping.ril.parental.frequency | sort -k1,1 -k2,2n > RIL275_RelocaTEi.CombinedGFF.characterized.clean.mping.ril.parental.frequency.lower

echo "distort block"
python Segragation_block.py --input HEG4vsNB.Marker.Segregation.table
bedtools intersect -wo -a ../../Compare_unique/HEG4.ALL.mping.non-ref.gff -b HEG4vsNB.Marker.Segregation.distort.bed | less -S

