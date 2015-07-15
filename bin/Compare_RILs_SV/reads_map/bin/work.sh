echo "Map reads to mPing and 10kb flanking regions, non_ref mPings were inserted into genome"
cd ../input
samtools faidx MSU_r7.Pseudo_mPing.fa
perl ~/software/bin/fastaDeal.pl --attr id:len MSU_r7.Pseudo_mPing.fa > MSU_r7.Pseudo_mPing.chrlen
#bedtools slop -i MSU_r7.Pseudo_mPing.gff -g MSU_r7.Pseudo_mPing.chrlen -b 10000 > MSU_r7.Pseudo_mPing.10kb_flank.bed
#bedtools getfasta -fi MSU_r7.Pseudo_mPing.fa -bed MSU_r7.Pseudo_mPing.10kb_flank.bed -fo MSU_r7.Pseudo_mPing.10kb_flank.fa

#get linked mPing within 100kb, deal with first. add 2 to ping position in MSU_r7.Pseudo_mPing.gff
python Get_linked_mPing_gff.py --gff ../input/MSU_r7.Pseudo_mPing.gff --distance ../input/mPing_dist.100kb.list.sorted
bedtools slop -i MSU_r7.Pseudo_mPing.linked_100kb.gff -g MSU_r7.Pseudo_mPing.chrlen -b 10000 > MSU_r7.Pseudo_mPing.linked_100kb.10kb_flank.bed
#100kb linked
bedtools getfasta -fi MSU_r7.Pseudo_mPing.fa -bed MSU_r7.Pseudo_mPing.linked_100kb.10kb_flank.bed -fo MSU_r7.Pseudo_mPing.linked_100kb.10kb_flank_temp.fa
sed 's/:/_/g' MSU_r7.Pseudo_mPing.linked_100kb.10kb_flank_temp.fa | sed 's/-/_/g' > MSU_r7.Pseudo_mPing.linked_100kb.10kb_flank.fa
rm MSU_r7.Pseudo_mPing.linked_100kb.10kb_flank_temp.fa
#100kb unlinked (other)
bedtools slop -i MSU_r7.Pseudo_mPing.other.gff -g MSU_r7.Pseudo_mPing.chrlen -b 10000 > MSU_r7.Pseudo_mPing.other.10kb_flank.bed
bedtools getfasta -fi MSU_r7.Pseudo_mPing.fa -bed MSU_r7.Pseudo_mPing.other.10kb_flank.bed -fo MSU_r7.Pseudo_mPing.other.10kb_flank_temp.fa
sed 's/:/_/g' MSU_r7.Pseudo_mPing.other.10kb_flank_temp.fa | sed 's/-/_/g' > MSU_r7.Pseudo_mPing.other.10kb_flank.fa
rm MSU_r7.Pseudo_mPing.other.10kb_flank_temp.fa


bwa index MSU_r7.Pseudo_mPing.10kb_flank.fa
python runRIL_bwa.py --input ../input/RILs_ALL_fastq




echo "For every mPing check both boundary if read cover the junction: cover, clipped reads (not cover), unsure"

