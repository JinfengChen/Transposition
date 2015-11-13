python Ping_genotype.py --input Ping.list

echo "275"
python Ping_genotype.py --input Ping.list --snp_map RILs_275/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam --output RILs_275_ping_genotype &

echo "230"
#bam files linked from 275, should get results from RIL275
python RIL230.py RILs_275_ping_genotype.table.ping_code.list RILs_230_ping_genotype.table.ping_code.list

echo "multilib"
python Ping_genotype.py --input Ping.list --snp_map RILs_multilib/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_multi_lib --output RILs_multilib_ping_genotype &
python Ping_genotype.py --input Ping_otherMarker14.list --snp_map RILs_multilib/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_multi_lib --output RILs_multilib_ping_other_marker14_genotype &
python Ping_genotype.py --input Ping_otherMarker.list --snp_map RILs_multilib/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_multi_lib --output RILs_multilib_ping_other_marker_genotype &

echo "correction"
python Ping_genotype.py --input Ping.list --snp_map RILs_correction/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_correction --output RILs_correction_ping_genotype &

echo "fix_ID"
python Ping_genotype.py --input Ping.list --snp_map RILs_fixID/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_fixID --output RILs_fixID_ping_genotype &


echo "14 marker"
python Ping_genotype.py --input Ping_otherMarker14.list --snp_map RILs_275/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam --output RILs_275_ping_other_marker14_genotype &
python RIL230.py RILs_275_ping_other_marker14_genotype.table.ping_code.list RILs_230_ping_other_marker14_genotype.table.ping_code.list
echo "Ping and other Marker to genotype RILs"
python Ping_genotype.py --input Ping_otherMarker.list --snp_map RILs_275/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam --output RILs_275_ping_other_marker_genotype &
python RIL230.py RILs_275_ping_other_marker_genotype.table.ping_code.list RILs_230_ping_other_marker_genotype.table.ping_code.list

echo "ping8 and mping13 marker"
python Ping_genotype.py --input Ping_mPing13.list --snp_map RILs_275/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam --output RILs_275_ping_mping13_genotype
python RIL230.py RILs_275_ping_mping13_genotype.table.ping_code.list RILs_230_ping_mping13_genotype.table.ping_code.list
python Ping_genotype.py --input Ping_mPing13.list --snp_map RILs_multilib/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_multi_lib --output RILs_multilib_ping_mping13_genotype &
python Ping_genotype.py --input Ping_mPing13.list --snp_map RILs_correction/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_correction --output RILs_correction_ping_mping13_genotype &
python Ping_genotype.py --input Ping_mPing13.list --snp_map RILs_fixID/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_fixID --output RILs_fixID_ping_mping13_genotype &


