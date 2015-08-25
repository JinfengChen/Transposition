echo "self"
#python IGV_mPing_Landrace.py --strain HEG4
#python IGV_mPing_Landrace.py --strain EG4
#python IGV_mPing_Landrace.py --strain A123
#python IGV_mPing_Landrace.py --strain A119
echo "cross"
python IGV_mPing_Landrace.py --strain HEG4 --cross
python IGV_mPing_Landrace.py --strain EG4 --cross
python IGV_mPing_Landrace.py --strain A123 --cross
python IGV_mPing_Landrace.py --strain A119 --cross

echo "unique mPing used to check cross pairs"
#bedtools intersect -a HEG4.hom.gff -b EG4.hom.gff -v > HEG4.unique.gff
#bedtools intersect -a EG4.hom.gff -b HEG4.hom.gff -v > EG4.unique.gff
#bedtools intersect -a A123.hom.gff -b A119.hom.gff -v > A123.unique.gff
#bedtools intersect -a A119.hom.gff -b A123.hom.gff -v > A119.unique.gff
