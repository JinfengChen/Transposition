cat MSU7.centromere.txt MSU7.chr_start_end.txt > MSU7.black.bed
awk -F"\t" '{print "Chr"$1"\t"$2"\t"$3}' NB_P.high_low_500bp.bedgraph > NB_P.high_low_500bp_chr.bedgraph
awk -F"\t" '{print "Chr"$1"\t"$2"\t"$3}' HEG4_P.high_low_500bp.bedgraph > HEG4_P.high_low_500bp_chr.bedgraph
