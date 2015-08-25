grep "hom" ../../../Compare_unique/RIL275_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff > RIL.gff
grep -e "het" -e "som" ../../../Compare_unique/RIL275_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff > Somatic.gff

echo "Landrace split, sofia version"
grep "HEG4" Strains.gff > Landrace_sofia/HEG4.gff
grep "A123" Strains.gff > Landrace_sofia/A123.gff
grep "A119" Strains.gff > Landrace_sofia/A119.gff
grep -e "\,EG4" -e "\=EG4" Strains.gff > Landrace_sofia/EG4.gff

echo "Landrace RelocaTE2"
bedtools intersect -a Landrace_RelocaTE2/HEG4.hom.gff -b Strains.gff -v > Landrace_RelocaTE2/HEG4.add.gff
bedtools intersect -a Landrace_RelocaTE2/EG4.hom.gff -b Strains.gff -v > Landrace_RelocaTE2/EG4.add.gff
bedtools intersect -a Landrace_RelocaTE2/A123.hom.gff -b Strains.gff -v > Landrace_RelocaTE2/A123.add.gff
bedtools intersect -a Landrace_RelocaTE2/A119.hom.gff -b Strains.gff -v > Landrace_RelocaTE2/A119.add.gff
echo "remove duplicat in all four and add these to Strains.gff"
cd Landrace_RelocaTE2 
cat HEG4.add.gff EG4.add.specific.gff A123.add.gff A119.add.gff > Strains.add.gff
cat Strains.add.gff ../Landrace_sofia/Strains.gff > ../Strains.gff
