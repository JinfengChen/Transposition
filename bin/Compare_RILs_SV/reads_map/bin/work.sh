echo "Map reads to mPing and 10kb flanking regions, non_ref mPings were inserted into genome"
cd ../input
samtools faidx MSU_r7.Pseudo_mPing.fa
perl ~/software/bin/fastaDeal.pl --attr id:len MSU_r7.Pseudo_mPing.fa > MSU_r7.Pseudo_mPing.chrlen
#bedtools slop -i MSU_r7.Pseudo_mPing.gff -g MSU_r7.Pseudo_mPing.chrlen -b 10000 > MSU_r7.Pseudo_mPing.10kb_flank.bed
#bedtools getfasta -fi MSU_r7.Pseudo_mPing.fa -bed MSU_r7.Pseudo_mPing.10kb_flank.bed -fo MSU_r7.Pseudo_mPing.10kb_flank.fa

#get linked mPing within 100kb, deal with first. add 2 to ping position in MSU_r7.Pseudo_mPing.gff
python Get_linked_mPing_gff.py --gff ../input/MSU_r7.Pseudo_mPing.gff --distance ../input/mPing_dist.100kb.list.sorted
python Get_linked_mPing_gff_ref.py --gff ../input/HEG4.ALL.mping.non-ref.gff --distance ../input/mPing_dist2.100kb.list.sorted
python Get_linked_mPing_gff_ref.py --gff ../input/HEG4.ALL.mping.non-ref.AF0.1.gff --distance ../input/mPing_dist2.100kb.list.sorted
#bedtools slop -i MSU_r7.Pseudo_mPing.linked_100kb.gff -g MSU_r7.Pseudo_mPing.chrlen -b 10000 > MSU_r7.Pseudo_mPing.linked_100kb.10kb_flank.bed
#100kb linked
#bedtools getfasta -fi MSU_r7.Pseudo_mPing.fa -bed MSU_r7.Pseudo_mPing.linked_100kb.10kb_flank.bed -fo MSU_r7.Pseudo_mPing.linked_100kb.10kb_flank_temp.fa
#sed 's/:/_/g' MSU_r7.Pseudo_mPing.linked_100kb.10kb_flank_temp.fa | sed 's/-/_/g' > MSU_r7.Pseudo_mPing.linked_100kb.10kb_flank.fa
#rm MSU_r7.Pseudo_mPing.linked_100kb.10kb_flank_temp.fa
#python runRIL_bwa.py --input ../input/RILs_ALL_fastq --ref ../input/MSU_r7.Pseudo_mPing.linked_100kb.10kb_flank.fa > log 2>&1 &

#100kb unlinked (other)
#bedtools slop -i MSU_r7.Pseudo_mPing.other.gff -g MSU_r7.Pseudo_mPing.chrlen -b 10000 > MSU_r7.Pseudo_mPing.other.10kb_flank.bed
#bedtools getfasta -fi MSU_r7.Pseudo_mPing.fa -bed MSU_r7.Pseudo_mPing.other.10kb_flank.bed -fo MSU_r7.Pseudo_mPing.other.10kb_flank_temp.fa
#sed 's/:/_/g' MSU_r7.Pseudo_mPing.other.10kb_flank_temp.fa | sed 's/-/_/g' > MSU_r7.Pseudo_mPing.other.10kb_flank.fa
#rm MSU_r7.Pseudo_mPing.other.10kb_flank_temp.fa


#bwa index MSU_r7.Pseudo_mPing.10kb_flank.fa
#python runRIL_bwa.py --input ../input/RILs_ALL_fastq

echo "get unmapped reads and mPing regions reads, map to Pseudogenome"
echo "Prepare unmapped and mPing region reads"
python PrepareFastq.py --bam ../input/RILs_ALL_bam > log 2>&1 &
echo "Run bwa mapping to pseudogenome"
python runRIL_bwa.py --input ../input/RILs_ALL_unmapped_mping_fastq > log 2>&1 &
python runRIL_bwa.py --input ../input/RILs_ALL_unmapped_mping_fastq2 > log 2>&1 &
python runRIL_bwa.py --input ../input/RILs_ALL_unmapped_mping_fastq3 > log 2>&1 &
perl /rhome/cjinfeng/BigData/software/bin/qsub-pbs.pl --maxjob 30 --lines 1 --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no RIL_bwa.sh > log1 2>&1 &
echo "Redo mapping using all RILs mPing AF>0.1, 20150827:RILs_ALL_unmapped_mping_bam"

echo "chech finished"
ls RILs_ALL_bam/*.bam | sed 's/RILs_ALL_bam\///' | sed 's/.bam//' | sed 's/GN/RIL/' > rils.list0
ls RILs_ALL_unmapped_mping_bam/*.bam | sed 's/RILs_ALL_unmapped_mping_bam\///' | sed 's/.bam//' > rils.list1
ls RILs_ALL_unmapped_mping_fastq/*/*_1.fq | sed 's/RILs_ALL_unmapped_mping_fastq\/RIL.*\///' | sed 's/_1.fq//' > rils.list2

echo "For every mPing check both boundary if read cover the junction: cover, clipped reads (not cover), unsure"
python mPing_Boundary_Coverage.py --bam_ref ../input/RILs_ALL_bam --bam_pseudo ../input/RILs_ALL_unmapped_mping_bam --gff_ref ../input/HEG4.ALL.mping.non-ref.gff --gff_pseudo ../input/MSU_r7.Pseudo_mPing.gff > log 2>&1 &
python mPing_Boundary_Coverage.py --bam_ref ../input/RILs_ALL_bam --bam_pseudo ../input/RILs_ALL_unmapped_mping_bam --gff_ref ../input/HEG4.ALL.mping.non-ref.100.gff --gff_pseudo ../input/MSU_r7.Pseudo_mPing.100.gff > log 2>&1 &
python mPing_Boundary_Coverage.py --bam_ref ../input/RILs_ALL_bam --bam_pseudo ../input/RILs_ALL_unmapped_mping_bam --gff_ref ../input/HEG4.ALL.mping.non-ref.high.gff --gff_pseudo ../input/MSU_r7.Pseudo_mPing.high.gff > log 2>&1 &
python mPing_Boundary_Coverage.py --bam_ref ../input/RILs_ALL_bam --bam_pseudo ../input/RILs_ALL_unmapped_mping_bam --gff_ref ../input/HEG4.ALL.mping.non-ref.debug.gff --gff_pseudo ../input/MSU_r7.Pseudo_mPing.debug.gff > log 2>&1 &

echo "add SNPmap to comfirm binmap"
#test run
python mPing_Boundary_Coverage.py --bam_ref ../input/RILs_ALL_bam --bam_pseudo ../input/RILs_ALL_unmapped_mping_bam_HEG4_mPing --gff_ref ../input/HEG4.ALL.mping.non-ref.debug.gff --gff_pseudo ../input/MSU_r7.Pseudo_mPing.SNPmap.gff > log 2>&1
#update to use mPing in RIL with AF>0.1
python mPing_Boundary_Coverage.py --bam_ref ../input/RILs_ALL_bam --bam_pseudo ../input/RILs_ALL_unmapped_mping_bam --gff_ref ../input/RIL275_RelocaTEi.CombinedGFF.characterized.AF0.1.gff --gff_pseudo ../input/MSU_r7.Pseudo_mPing_RILs.gff > log 2>&1 
#check process
grep "^../input/RILs_ALL_bam/GN" mPing_boundary.bamcheck_ref.txt | cut -d" " -f1 | uniq | sort | uniq | wc -l

echo "summarize results from mPing_Boundary_Coverage.py"
python Sum_linked_mPing_status.py --dir mPing_boundary_mPing --distance ../input/mPing_dist.100kb.list.sorted | sort -k2,2n > mPing_boundary.linked_100kb.table.txt
python Sum_linked_mPing_status.py --dir mPing_boundary_mPing_manual --distance ../input/mPing_dist.100kb.list.sorted.1 --project mPing_boundary.linked_100kb_manual
python Sum_linked_mPing_status.py --dir mPing_boundary_mPing --distance ../input/mPing_dist.100kb.list.sorted.1 --project mPing_boundary.linked_100kb_debug
python Sum_linked_mPing_status.py --dir mPing_boundary_mPing --distance ../input/mPing_dist.100kb.list.sorted.2 --project mPing_boundary.linked_100kb_debug
python Sum_linked_mPing_status.py --dir mPing_boundary_mPing --distance ../input/mPing_dist.100kb.list.sorted --project mPing_boundary.linked_100kb_debug
#mPing_dist1: modified to count both side to shortest distance, only one change within 100kb
python Sum_linked_mPing_status.py --dir mPing_boundary_mPing --distance ../input/mPing_dist1.100kb.list.sorted --project mPing_boundary.linked_100kb_debug1
python Sum_linked_mPing_status.py --dir mPing_boundary_mPing --distance ../input/mPing_dist1.50Mb.list.sorted --project mPing_boundary.linked_50Mb_debug1
#mPing_dist2: removed these AF<0.1
python Sum_linked_mPing_status.py --dir mPing_boundary_mPing --distance ../input/mPing_dist2.100kb.list.sorted --project mPing_boundary.linked_100kb_debug2
python Sum_linked_mPing_status.py --dir mPing_boundary_mPing --distance ../input/mPing_dist2.50Mb.list.sorted --project mPing_boundary.linked_50Mb_debug2

echo "ave in 100 kb interval"
python avg_interval.py --input mPing_boundary.linked_50Mb_debug2.table_clean.txt > mPing_boundary.linked_50Mb_debug2.table_clean.sum
python avg_interval.py --input mPing_boundary.linked_50Mb_debug1.table_clean.txt > mPing_boundary.linked_50Mb_debug1.table_clean.sum
cat avg_interval_test.R | R --slave

echo "update to use only mPing, distance and excision, not pairs which have many duplicate"
python avg_interval.py --excision mPing_boundary_mPing_debug2_results/mPing_boundary.linked_50Mb_debug2.mping_excision.list --distance ../input/mPing_dist2.50Mb.list.sorted
cat avg_interval_test.R | R --slave | grep "p-value"

echo "plot point view of excision and distance"
python Sum_linked_mPing_statusV2.py --dir ../../../Compare_excision_transposase/bin/High_excision_csv_Ping --distance ../input/mPing_dist_RIL_AF0.1.50Mb.list.sorted --project mPing_boundary.linked_50Mb_debug2
awk '$8<20 && $9<20 && $10<20' mPing_boundary.linked_50Mb_debug2.table_clean.txt > mPing_boundary.linked_50Mb_debug2.table_clean.rough_removed_early_insertion.txt
cat mPing_boundary.linked_50Mb_debug2.table_clean.point_plot.R | R --slave
#411 mping and 881 excision including early excision for this pipeline, compare details with former
python Sum_excision_distance.py --dir ../../../Compare_excision_transposase/bin/High_excision_csv_Ping --distance ../input/mPing_dist_RIL_AF0.1.50Mb.list.sorted --project mPing_boundary.linked_50Mb_debug2
cat mPing_boundary.linked_50Mb_debug2.table_clean.point_view.R | R --slave
cat mPing_boundary.linked_50Mb_debug2.table_clean.proportion_accumulation.R | R --slave
#ks test p-value=0.0001, difference between excision and control
cat mPing_boundary.linked_50Mb_debug2.table_clean.proportion_accumulation.test.R | R --slave
