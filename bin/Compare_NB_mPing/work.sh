echo "all HEG4 mPing, non-ref + shared"
cat HEG4.ALL.mping.non-ref.gff HEG4.mping.shared.gff | sort -k1,1 -k4,4n > HEG4.ALL.mPing.gff
echo "all NB mPing, ref_only + shared"
cat HEG4.mping.ref_only.gff HEG4.mping.shared.gff | sort -k1,1 -k4,4n > NB.ALL.mPing.gff
echo "all parent mPing in RILs, non-ref + ref_only + shared"
cat HEG4.ALL.mping.non-ref.gff HEG4.mping.shared.gff HEG4.mping.ref_only.gff | sort -k1,1 -k4,4n > Parent.ALL.mPing.gff

echo "distance"
perl mPing_dist_ref.pl --input NB.ALL.mPing.gff | sort -k3,3n

