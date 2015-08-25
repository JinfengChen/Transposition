echo "all HEG4 mPing, non-ref + shared"
cat HEG4.ALL.mping.non-ref.gff HEG4.mping.shared.gff | sort -k1,1 -k4,4n > HEG4.ALL.mPing.gff
echo "all NB mPing, ref_only + shared"
cat HEG4.mping.ref_only.gff HEG4.mping.shared.gff | sort -k1,1 -k4,4n > NB.ALL.mPing.gff
echo "all parent mPing in RILs, non-ref + ref_only + shared"
cat HEG4.ALL.mping.non-ref.gff HEG4.mping.shared.gff HEG4.mping.ref_only.gff | sort -k1,1 -k4,4n > Parent.ALL.mPing.gff

echo "distance"
perl mPing_dist_ref.pl --input NB.ALL.mPing.gff | sort -k3,3n

perl mPing_dist_ref.pl --input NB.ALL.mPing.gff | sort -k3,3n > NB.ALL.mPing.distance.txt
perl mPing_dist_ref.pl --input HEG4.ALL.mPing.gff | sort -k3,3n > HEG4.ALL.mPing.distance.txt
perl mPing_dist_ref.pl --input Parent.ALL.mPing.gff | sort -k3,3n > Parent.ALL.mPing.distance.txt
