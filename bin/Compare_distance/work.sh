perl mPing_dist.pl --input HEG4.ALL.mping.non-ref.gff
sort -k3,3n -k1,1n -k2,2n mPing_dist.100kb.list > mPing_dist.100kb.list.sorted
#fix shortest distance from both side
sort -k3,3n -k1,1n -k2,2n mPing_dist1.100kb.list | uniq > mPing_dist1.100kb.list.sorted
sort -k3,3n -k1,1n -k2,2n mPing_dist1.50Mb.list | uniq > mPing_dist1.50Mb.list.sorted
#only these mPing with AF>0.1
perl mPing_dist.pl --input HEG4.ALL.mping.non-ref.AF0.1.gff
sort -k3,3n -k1,1n -k2,2n mPing_dist2.100kb.list | uniq > mPing_dist2.100kb.list.sorted
sort -k3,3n -k1,1n -k2,2n mPing_dist2.50Mb.list | uniq > mPing_dist2.50Mb.list.sorted
#RIL mPing AF>0.1 
perl mPing_dist.pl --input RIL275_RelocaTEi.CombinedGFF.characterized.AF0.1.gff
sort -k3,3n -k1,1n -k2,2n mPing_dist_RIL_AF0.1.50Mb.list | uniq > mPing_dist_RIL_AF0.1.50Mb.list.sorted

