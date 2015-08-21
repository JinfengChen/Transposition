grep "hom" ../../../Compare_unique/RIL275_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff > RIL.gff
grep -e "het" -e "som" ../../../Compare_unique/RIL275_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff > Somatic.gff

echo "Landrace split, sofia version"
grep "HEG4" Strains.gff > Landrace_sofia/HEG4.gff
grep "A123" Strains.gff > Landrace_sofia/A123.gff
grep "A119" Strains.gff > Landrace_sofia/A119.gff
grep -e "\,EG4" -e "\=EG4" Strains.gff > Landrace_sofia/EG4.gff

