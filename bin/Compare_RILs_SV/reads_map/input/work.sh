cd ../input
samtools faidx MSU_r7.Pseudo_mPing.fa
perl ~/software/bin/fastaDeal.pl --attr id:len MSU_r7.Pseudo_mPing.fa > MSU_r7.Pseudo_mPing.chrlen
bedtools slop -i MSU_r7.Pseudo_mPing.gff -g MSU_r7.Pseudo_mPing.chrlen -b 10000 > MSU_r7.Pseudo_mPing.10kb_flank.bed
bedtools getfasta -fi MSU_r7.Pseudo_mPing.fa -bed MSU_r7.Pseudo_mPing.10kb_flank.bed -fo MSU_r7.Pseudo_mPing.10kb_flank.fa
bwa index MSU_r7.Pseudo_mPing.10kb_flank.fa

echo "debug site"
grep "68235" MSU_r7.Pseudo_mPing.gff > MSU_r7.Pseudo_mPing.high.gff 
grep "68336" MSU_r7.Pseudo_mPing.gff >> MSU_r7.Pseudo_mPing.high.gff 
grep "363090" MSU_r7.Pseudo_mPing.gff >> MSU_r7.Pseudo_mPing.high.gff
grep "363123" MSU_r7.Pseudo_mPing.gff >> MSU_r7.Pseudo_mPing.high.gff

grep "875060" HEG4.ALL.mping.non-ref.gff > HEG4.ALL.mping.non-ref.debug.gff 
grep "875074" HEG4.ALL.mping.non-ref.gff >> HEG4.ALL.mping.non-ref.debug.gff 
grep "260838" HEG4.ALL.mping.non-ref.gff >> HEG4.ALL.mping.non-ref.debug.gff 
grep "260848" HEG4.ALL.mping.non-ref.gff >> HEG4.ALL.mping.non-ref.debug.gff 

grep "875060" MSU_r7.Pseudo_mPing.gff > MSU_r7.Pseudo_mPing.debug.gff
grep "875074" MSU_r7.Pseudo_mPing.gff >> MSU_r7.Pseudo_mPing.debug.gff
grep "260838" MSU_r7.Pseudo_mPing.gff >> MSU_r7.Pseudo_mPing.debug.gff
grep "260848" MSU_r7.Pseudo_mPing.gff >> MSU_r7.Pseudo_mPing.debug.gff

grep "87506" mPing_dist.100kb.list.sorted > mPing_dist.100kb.list.sorted.2
grep "260838" mPing_dist.100kb.list.sorted >> mPing_dist.100kb.list.sorted.2
