#excision found in HEG4 and EG4 are almost right with only false repeat_Chr3_16551323_16551325
python Excision_Landrace.py --ref_gff ../input/HEG4.hom.gff --qry_gff ../input/EG4.hom.gff --ref_bam ../input/HEG4.bam --qry_bam ../input/EG4.bam --project HEG4vsEG4
sort HEG4vsEG4.mPing.Excision.table | grep "footprint"
#all excision found in A123 and A119 are false negative due to existing indels in A123/A119 or shared indels.
python Excision_Landrace.py --ref_gff ../input/A123.hom.gff --qry_gff ../input/A119.hom.gff --ref_bam ../input/A123.bam --qry_bam ../input/A119.bam --project A123vsA119
sort A123vsA119.mPing.Excision.table | grep "footprint"


echo "distance"
perl mPing_dist_ref.pl --input ../input/HEG4.hom.gff | sort -k1,1 -k3,3n > HEG4.mPing.distance
perl mPing_dist_ref.pl --input ../input/EG4.hom.gff | sort -k1,1 -k3,3n > EG4.mPing.distance
perl mPing_dist_ref.pl --input ../input/A123.hom.gff | sort -k1,1 -k3,3n > A123.mPing.distance
perl mPing_dist_ref.pl --input ../input/A119.hom.gff | sort -k1,1 -k3,3n > A119.mPing.distance
cat mPing_Distance_landrace.R | R --slave

