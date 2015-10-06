python Ping_genotype.py --input Ping.list

echo "275"
python Ping_genotype.py --input Ping.list --snp_map RILs_275/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam --output RILs_275_ping_genotype &

echo "multilib"
python Ping_genotype.py --input Ping.list --snp_map RILs_multilib/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_multi_lib --output RILs_multilib_ping_genotype &

echo "correction"
python Ping_genotype.py --input Ping.list --snp_map RILs_correction/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_correction --output RILs_correction_ping_genotype &

echo "fix_ID"
python Ping_genotype.py --input Ping.list --snp_map RILs_fixID/MPR.geno.data --bam /rhome/cjinfeng/BigData/00.RD/RILs/QTL_pipe/input/fastq/RILs_ALL_bam_fixID --output RILs_fixID_ping_genotype &
