echo "Map reads to mPing and 10kb flanking regions, non_ref mPings were inserted into genome"
cd ../input
samtools faidx MSU_r7.Pseudo_mPing.fa
perl ~/software/bin/fastaDeal.pl --attr id:len MSU_r7.Pseudo_mPing.fa > MSU_r7.Pseudo_mPing.chrlen
#bedtools slop -i MSU_r7.Pseudo_mPing.gff -g MSU_r7.Pseudo_mPing.chrlen -b 10000 > MSU_r7.Pseudo_mPing.10kb_flank.bed
#bedtools getfasta -fi MSU_r7.Pseudo_mPing.fa -bed MSU_r7.Pseudo_mPing.10kb_flank.bed -fo MSU_r7.Pseudo_mPing.10kb_flank.fa

#get linked mPing within 100kb, deal with first. add 2 to ping position in MSU_r7.Pseudo_mPing.gff
#python Get_linked_mPing_gff.py --gff ../input/MSU_r7.Pseudo_mPing.gff --distance ../input/mPing_dist.100kb.list.sorted
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
perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --maxjob 30 --lines 1 --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=20G --convert no RIL_bwa.sh > log1 2>&1 &

echo "chech finished"
ls RILs_ALL_bam/*.bam | sed 's/RILs_ALL_bam\///' | sed 's/.bam//' | sed 's/GN/RIL/' > rils.list0
ls RILs_ALL_unmapped_mping_bam/*.bam | sed 's/RILs_ALL_unmapped_mping_bam\///' | sed 's/.bam//' > rils.list1
ls RILs_ALL_unmapped_mping_fastq/*/*_1.fq | sed 's/RILs_ALL_unmapped_mping_fastq\/RIL.*\///' | sed 's/_1.fq//' > rils.list2


echo "For every mPing check both boundary if read cover the junction: cover, clipped reads (not cover), unsure"

