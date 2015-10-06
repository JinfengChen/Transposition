echo "snp genotype around high excision mPing in RILs excised"
python high_excision_mPing_geneotype.py --input mping.excision.draw.highexcision.1 --type excision --output high_excision_mping.genotype_excision.table.txt

echo "snp genotype around pairs of high excision mPing in RILs excised"
python high_excision_mPing_geneotype.py --input mping.excision.draw.highexcision.pairs --type insertion --output high_excision_mping.genotype_insertion.table.txt
