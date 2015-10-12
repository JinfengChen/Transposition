#perl runcnvnoter_rils.pl --bam /rhome/cjinfeng/BigData/00.RD/Variations/SV/pindel/input/RILs_ALL_bam/ --project RIL_275 > log 2>&1 &
#perl runcnvnoter_rils.pl --bam /rhome/cjinfeng/BigData/00.RD/Variations/SV/pindel/input/RILs_ALL_bam_core/ --project RIL_230 > log 2>&1 &


echo "run rils"
qsub -q highmem cnvnator_rils.sh
awk '$4 < 0.05 || $4>2' GN162.readdepth.bed > GN162.filtered.bed
