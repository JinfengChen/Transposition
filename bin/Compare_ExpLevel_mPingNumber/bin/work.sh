echo "mPing number asscoiated with gene"
python mPing_Number.py --input ../input/RIL275_RelocaTEi.CombinedGFF.characterized.unique_mPing.gff.mPing.annotation | sort -k2,2n > mPing_Number.list
echo "merged with gene expression level"
python listdiff.py mPing_Number.list ../input/all_sample.expr_value.gene_exp.list | sort -k2,2rn > mPing_Number.Exp.list

