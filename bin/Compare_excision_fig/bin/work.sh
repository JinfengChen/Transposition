echo "scatterhist plot of excision vs. distance"
python Excision_Distance.py --excision1 ../input/mping.excision.table --distance ../input/mPing_dist_RIL_AF0.1.50Mb.list.sorted --gff ../input/HEG4.ALL.mping.non-ref.AF0.1.gff
python Excision_Distance.py --excision2 ../input/mPing_boundary.linked_50Mb_debug2.mping_excision.list --distance ../input/mPing_dist_RIL_AF0.1.50Mb.list.sorted --gff ../input/HEG4.ALL.mping.non-ref.AF0.1.gff

echo "boxplot of distance among group: high excision vs. all mping vs. control"
python Excision_Distance_Boxplot.py --highexcision ../input/mping.excision.draw.highexcision --distance ../input/mPing_dist_RIL_AF0.1.50Mb.list.sorted --gff ../input/HEG4.ALL.mping.non-ref.AF0.1.gff
echo "new excision version, potential high excision"
awk '$2>5' Excision_distance.matrix.2.txt > ../input/mping.excision.draw.highexcision.2
python Excision_Distance_Boxplot.2.py --highexcision ../input/mping.excision.draw.highexcision.2 --distance ../input/mPing_dist_RIL_AF0.1.50Mb.list.sorted --gff ../input/HEG4.ALL.mping.non-ref.AF0.1.gff
