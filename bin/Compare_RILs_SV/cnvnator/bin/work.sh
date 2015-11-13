#perl runcnvnoter_rils.pl --bam /rhome/cjinfeng/BigData/00.RD/Variations/SV/pindel/input/RILs_ALL_bam/ --project RIL_275 > log 2>&1 &
#perl runcnvnoter_rils.pl --bam /rhome/cjinfeng/BigData/00.RD/Variations/SV/pindel/input/RILs_ALL_bam_core/ --project RIL_230 > log 2>&1 &


echo "run rils"
qsub -q highmem cnvnator_rils.sh

echo "Filter"
python Filter_CNVnator.py --input RIL275_correction
#or 
qsub run_filter.sh

echo "prepare draw bam and validation table of manual check"
python Prepare_validation.py --input RIL230_core_filtered
sort -k1,1 -k2,2n RIL230_core_filtered.validation_table.txt > RIL230_core_filtered.validation_table.sorted.txt
python get_excision_bam.py --input RIL230_core_filtered.draw.txt  --output RIL230_core_filtered.draw

