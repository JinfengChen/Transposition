python mPing_intergenic_gene_list.py --input A119.mRNA.intersect | awk '$NF <= -1500 && $NF >= -2000' > A119.mRNA.intersect.1.5_2.0kb.gene.list
python mPing_intergenic_gene_list.py --input A119.mRNA.intersect | awk '$NF <= -1000 && $NF >= -1500' > A119.mRNA.intersect.1.0_1.5kb.gene.list
python mPing_intergenic_gene_list.py --input HEG4.mRNA.intersect | awk '$NF <= -1500 && $NF >= -2000' > HEG4.mRNA.intersect.1.5_2.0kb.gene.list
python mPing_intergenic_gene_list.py --input HEG4.mRNA.intersect | awk '$NF <= -1000 && $NF >= -1500' > HEG4.mRNA.intersect.1.0_1.5kb.gene.list

python mPing_intergenic_subsample_HEG4.py
paste ../RIL.mRNA.5primer.distance.distr ../Strains.mRNA.5primer.distance.distr HEG4.mRNA.5primer.distance.distr A119.mRNA.5primer.distance.distr > Table.mRNA.5primer.distance.distr.txt

echo "expression"
python Get_Gene_Expression.py --gff HEG4.mRNA.intersect.1.5_2.0kb.gene.list --e ~/BigData/02.Transcription/Transcriptome/bin/DiffExpression/NB_HEG4_EG4.table
python Get_Gene_Expression.py --gff HEG4.mRNA.intersect.1.0_1.5kb.gene.list --e ~/BigData/02.Transcription/Transcriptome/bin/DiffExpression/NB_HEG4_EG4.table

