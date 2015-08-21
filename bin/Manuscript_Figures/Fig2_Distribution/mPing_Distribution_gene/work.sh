#cat mping_position_breakY_withsim.R | R --slave
#cat mping_intergenic_5distance_withsim.R | R --slave
#cat mping_intergenic_3distance_withsim.R | R --slave
#cat mping_intron_distance_point_withsim.R | R --slave
#cat mping_intron_length_point_withsim.R | R --slave 

dirname=/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Figures/Fig2_Distribution/mPing_distr/
echo "Intron length"
grep "intron" /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Figures/Fig2_Distribution/mPing_distr/Somatic.intersect | awk '{print $14-$13}' > Somatic.intron.length
grep "intron" /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Figures/Fig2_Distribution/mPing_distr/Strains.intersect | awk '{print $14-$13}' > Strains.intron.length
grep "intron" /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Figures/Fig2_Distribution/mPing_distr/RIL.intersect | awk '{print $14-$13}' > RIL.intron.length
grep "intron" /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Manuscript_Figures/Fig2_Distribution/mPing_distr/Simulation.intersect | awk '{print $14-$13}' > Simulation.intron.length
cat mPing_intron_length_boxplot.R | R --slave

echo "Intron distance"
perl -e 'while(<>){@unit=split("\t");if($unit[11]=~/intron/){$dist = abs($unit[3]-$unit[12]) > abs($unit[13]-$unit[4]) ? abs($unit[3]-$unit[12]) : abs($unit[13]-$unit[4]); print "$dist\n"}}' $dirname/Strains.intersect > Strains.intron.distance
perl -e 'while(<>){@unit=split("\t");if($unit[11]=~/intron/){$dist = abs($unit[3]-$unit[12]) > abs($unit[13]-$unit[4]) ? abs($unit[3]-$unit[12]) : abs($unit[13]-$unit[4]); print "$dist\n"}}' $dirname/Somatic.intersect > Somatic.intron.distance
perl -e 'while(<>){@unit=split("\t");if($unit[11]=~/intron/){$dist = abs($unit[3]-$unit[12]) > abs($unit[13]-$unit[4]) ? abs($unit[3]-$unit[12]) : abs($unit[13]-$unit[4]); print "$dist\n"}}' $dirname/Simulation.intersect > Simulation.intron.distance
perl -e 'while(<>){@unit=split("\t");if($unit[11]=~/intron/){$dist = abs($unit[3]-$unit[12]) > abs($unit[13]-$unit[4]) ? abs($unit[3]-$unit[12]) : abs($unit[13]-$unit[4]); print "$dist\n"}}' $dirname/RIL.intersect > RIL.intron.distance
cat mPing_intron_distance_boxplot.R | R --slave
