echo "Annotation insertion, unique"
perl AnnoTEsite.pl --te ../input/RIL275_RelocaTEi.CombinedGFF.characterized.clean.unique_mPing.gff --gene /rhome/cjinfeng/HEG4_cjinfeng/seqlib/GFF/MSU7.gene.gff > log 2>&1 &

echo "simple summary"
#1847 gene have more than 1 allele that cauased by mPing (within 3 kb)
awk '$7>-3000 && $7<3000' RIL275_RelocaTEi.mPing.annotation | cut -f6 | sort | uniq -d | wc -l
#3083 gene have mPing within 3kb
awk '$7>-3000 && $7<3000' RIL275_RelocaTEi.mPing.annotation | cut -f6 | sort | uniq | wc -l

echo "summary single gene table"

