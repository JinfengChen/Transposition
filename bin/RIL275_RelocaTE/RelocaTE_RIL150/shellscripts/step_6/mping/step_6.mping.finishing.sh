
`mkdir -p /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/all_files`

#combine confident insertions to one file
echo "TE	TSD	Exper	insertion_site	strand	left_flanking_read_count	right_flanking_read_count	left_flanking_seq	right_flanking_seq" > /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/temp
for i in `ls /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/*.mping.confident_nonref_insert.txt` ; do grep -v flanking_read_count $i >> /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/temp ; done
mv /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/temp /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/RIL150.mping.confident_nonref.txt
mv /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/*.mping.confident_nonref_insert.txt /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/all_files

#combine all insertions to one file
echo "TE	TSD	Exper	chromosome	insertion_site	strand	combined_read_count	right_flanking_read_count	left_flanking_read_count" > /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/temp2
for i in `ls /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/*.mping.all_nonref_insert.txt` ; do grep -v total $i | grep -v Note >> /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/temp2 ; done
mv /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/temp2 /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/RIL150.mping.all_nonref.txt
mv /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/*.mping.all_nonref_insert.txt /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/all_files

#combine confident insertions ref seqs to one file
for i in `ls /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/*.mping.confident_nonref_genomeflank.fa` ; do cat $i  >> /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/temp3 ; done
mv /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/temp3 /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/RIL150.mping.confident_nonref_genomeflanks.fa
mv /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/*.mping.confident_nonref_genomeflank.fa /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/all_files

#combine confident insertions gff to one file
echo "##gff-version 3" > /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/temp4
for i in `ls /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/*.mping.all_insert.gff` ; do grep -v gff $i  >> /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/temp4 ; done        
mv /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/temp4 /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/RIL150.mping.all_inserts.gff
mv /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/*.mping.all_insert.gff /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/all_files

#combine confident insertions reads to one file
for i in `ls /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/*.mping.confident_nonref_insert_reads_list.txt` ; do cat $i  >> /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/temp5 ; done
mv /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/temp5 /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/RIL150.mping.confident_nonref_reads_list.txt
mv /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/*.mping.confident_nonref_insert_reads_list.txt /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/results/all_files

# move other outfiles somewhere else
if [ -e /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/bowtie-build.out ] ; then 
  mv /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/bowtie-build.out /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/bowtie_aln/.
fi
mv /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/existingTE.blat.stdout /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/blat_output/.
mv /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/existingTE.blatout /bigdata/cjinfeng/00.RD/RILs/Transpostion/bin/RelocaTE_RIL150/mping/blat_output/.

