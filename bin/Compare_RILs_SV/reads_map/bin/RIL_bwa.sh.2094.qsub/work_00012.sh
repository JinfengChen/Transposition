perl /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/step1_Mapping.pl -ref /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/input/bwa_0.7/MSU_r7.Pseudo_mPing_RILs.fa -1 /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/input/RILs_ALL_unmapped_mping_fastq/RIL109/RIL109_1.fq -2 /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/input/RILs_ALL_unmapped_mping_fastq/RIL109/RIL109_2.fq -min 0 --max 500 -cpu 12 --tool bwa --project /bigdata/stajichlab/cjinfeng/00.RD/RILs/Transpostion/bin/Compare_RILs_SV/reads_map/bin/RIL109; echo This-Work-is-Completed!
