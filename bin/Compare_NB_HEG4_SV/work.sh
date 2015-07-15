echo "ALL possible mPing assoicated with SV in HEG4, 10 kb range of mPing and Ping"
bedtools window -w 10000 -a Parent.ALL.mPing.gff -b HEG4.pindel.gff | grep "Deletion" | grep "Reference-only" -v > Parent.ALL.mPing.SV_pindel.gff

echo "overlap of SV and mPing: 61"
bedtools window -w 10000 -a Parent.ALL.mPing.gff -b HEG4.Deletion.final.gff | grep "Deletion" | grep "Reference-only" -v | wc -l
