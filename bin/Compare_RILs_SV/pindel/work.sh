echo "pindel SV"
cp ~/BigData/00.RD/Variations/SV/pindel/bin/RIL275.pindel.raw.gff ./
cp ~/BigData/00.RD/Variations/SV/pindel/bin/RIL275.pindel.gff ./
cp ~/BigData/00.RD/Variations/SV/pindel/bin/RIL275.pindel.inv.gff ./

echo "3,3"
bedtools intersect -a RIL275.pindel.gff -b HEG4.Deletion.final.gff -v | sort -k1,1 -k4,4n > RIL275.pindel.temp.gff
bedtools intersect -a RIL275.pindel.temp.gff -b HEG4.pindel.deletion.gff -v | sort -k1,1 -k4,4n > RIL275.pindel.denovo.gff
bedtools window -w 1000 -a RIL275.pindel.denovo.gff -b HEG4.ALL.mping.non-ref.gff | less -S
echo "raw"
bedtools intersect -a RIL275.pindel.raw.gff -b HEG4.Deletion.final.gff -v | sort -k1,1 -k4,4n > RIL275.pindel.raw.temp.gff
bedtools intersect -a RIL275.pindel.raw.temp.gff -b HEG4.pindel.deletion.gff -v | sort -k1,1 -k4,4n > RIL275.pindel.raw.denovo.gff
bedtools window -w 1000 -a RIL275.pindel.raw.denovo.gff -b HEG4.ALL.mping.non-ref.gff | less -S
echo "inv"
bedtools window -w 1000 -a RIL275.pindel.inv.gff -b HEG4.ALL.mping.non-ref.gff | less -S
bedtools window -w 1000 -a RIL275.pindel.inv.gff -b Parent.ALL.mPing.gff | less -S

