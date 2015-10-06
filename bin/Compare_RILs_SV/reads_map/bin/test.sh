#MSU_r7.Pseudo_mPing.SNPmap.gff
#python mPing_Boundary_Coverage.py --bam_ref ../input/RILs_ALL_bam --bam_pseudo ../input/RILs_ALL_unmapped_mping_bam_HEG4_mPing --gff_ref ../input/HEG4.ALL.mping.non-ref.debug.gff --gff_pseudo ../input/MSU_r7.Pseudo_mPing.SNPmap.gff > log 2>&1

python mPing_Boundary_Coverage.py --bam_ref ../input/RILs_ALL_bam --bam_pseudo ../input/RILs_ALL_unmapped_mping_bam --gff_ref ../input/RIL275_RelocaTEi.CombinedGFF.characterized.AF0.1.gff --gff_pseudo ../input/MSU_r7.Pseudo_mPing_RILs.gff > log 2>&1
