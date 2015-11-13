echo "DHS region, not work well"
bedtools window -w 2000 -a ../input/MSU7.mRNA.gff -b ../input/GSM655033_Rice_Seedling_DHsites.MSU7.Corrected.gff > MSU7.mRNA_DHS.overlap &
python Gene_DHS.py --input MSU7.mRNA_DHS.overlap > Gene_DHS.list
python listdiff.py mPing_Number.list Gene_DHS.list | sort -k2,2rn > mPing_Number.DHS.list

echo "DHS read density"
python Create_5k_TSS.py --input ../input/MSU7.mRNA.gff > ../input/MSU7.mRNA.TSS_5k.bed
DHS_count.sh
cut -f4,7 ../input/MSU7.mRNA.TSS_5k.bed.count > Gene_DHS.count
python listdiff.py mPing_Number.list Gene_DHS.list Gene_DHS.count | sort -k2,2rn > mPing_Number.DHS.list

echo "merge"
python merge_DHS_EXP.py mPing_Number.DHS.list mPing_Number.Exp.list |sort -k2,2rn > mPing_Number.DHS_Exp.list

echo "fpkm"
cut -f1,8,9,10 ~/BigData/02.Transcription/Transcriptome/bin/DiffExpression_Jason/bin/Landrace/cuffnorm/genes.fpkm_table > gene.fpkm.list
python merge_DHS_EXP.py mPing_Number.DHS.list gene.fpkm.list | sort -k5,5rn > mPing_Number.DHS_FPKM.list

python Gene_DHS_promoter.py --input MSU7.mRNA_DHS.overlap > Gene_DHS_promoter.list
python listdiff.py mPing_Number.list Gene_DHS.list Gene_DHS_promoter.list | sort -k2,2rn > mPing_Number.DHS_promoter.list
python merge_DHS_EXP.py mPing_Number.DHS_promoter.list gene.fpkm.list | sort -k5,5rn > mPing_Number.DHS_pro_FPKM.list


echo "overlap of mPing with DHS"
python Distance2DHS.py

python TSDprofile_DHS.py ../input/RIL230_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff test.profile > test.sum &

