echo "scatterhist plot of excision vs. distance"
python Excision_Distance.py --excision1 ../input/mping.excision.table --distance ../input/mPing_dist_RIL_AF0.1.50Mb.list.sorted --gff ../input/HEG4.ALL.mping.non-ref.AF0.1.gff
python Excision_Distance.py --excision2 ../input/mPing_boundary.linked_50Mb_debug2.mping_excision.list --distance ../input/mPing_dist_RIL_AF0.1.50Mb.list.sorted --gff ../input/HEG4.ALL.mping.non-ref.AF0.1.gff
python Excision_Distance.py --excision1 ../input/Excision_newpipe_version1.footprint.list.txt --distance ../input/mPing_dist_RIL_AF0.1.50Mb.list.sorted --gff ../input/HEG4.ALL.mping.non-ref.AF0.1.gff --output Excision_distance.matrix_events
echo "run plot on mac and save pdf"
Excision_distance.matrix_events.1.R

echo "boxplot of distance among group: high excision vs. all mping vs. control"
python Excision_Distance_Boxplot.py --highexcision ../input/mping.excision.draw.highexcision --distance ../input/mPing_dist_RIL_AF0.1.50Mb.list.sorted --gff ../input/HEG4.ALL.mping.non-ref.AF0.1.gff
echo "new excision version, potential high excision"
awk '$2>5' Excision_distance.matrix.2.txt > ../input/mping.excision.draw.highexcision.2
python Excision_Distance_Boxplot.2.py --highexcision ../input/mping.excision.draw.highexcision.2 --distance ../input/mPing_dist_RIL_AF0.1.50Mb.list.sorted --gff ../input/HEG4.ALL.mping.non-ref.AF0.1.gff

echo "Events draw"
awk '$2>=5' ../input/Excision_newpipe_version1.footprint.list.txt >> ../input/Excision_newpipe_version1.footprint.high.txt
awk '$2>=4' ../input/Excision_newpipe_version1.footprint.list.txt >> ../input/Excision_newpipe_version1.footprint.high4.txt
python Excision_Distance_Boxplot.py --highexcision ../input/Excision_newpipe_version1.footprint.high.txt --distance ../input/mPing_dist_RIL_AF0.1.50Mb.list.sorted --gff ../input/HEG4.ALL.mping.non-ref.AF0.1.gff --output Excision_Events_distance.boxplot
python Excision_Distance_Boxplot.py --highexcision ../input/Excision_newpipe_version1.footprint.high4.txt --distance ../input/mPing_dist_RIL_AF0.1.50Mb.list.sorted --gff ../input/HEG4.ALL.mping.non-ref.AF0.1.gff --output Excision_Events_distance.boxplot4
echo "draw boxplot to compare high vs. all and test significant"
#cutoff 5
cat Excision_Events_distance.boxplot.R | R --slave
#cutoff 4
cat Excision_Events_distance.boxplot4.R | R --slave
