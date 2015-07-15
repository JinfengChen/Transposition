echo "DHS region, not work well"
bedtools window -w 2000 -a ../input/MSU7.mRNA.gff -b ../input/GSM655033_Rice_Seedling_DHsites.MSU7.Corrected.gff > MSU7.mRNA_DHS.overlap &
python Gene_DHS.py --input MSU7.mRNA_DHS.overlap > Gene_DHS.list
python listdiff.py mPing_Number.list Gene_DHS.list | sort -k2,2rn > mPing_Number.DHS.list

echo "DHS read density"
python Create_5k_TSS.py --input ../input/MSU7.mRNA.gff > ../input/MSU7.mRNA.TSS_5k.bed

