#echo "long intron gene with mping"
#grep "intron" RIL.intersect | awk '$14-$13 > 2000' | perl -e 'while(<>){$line=$_;if($line=~/Parent=(\w+\.\d+)/){print "$1\n"}}'

echo "RelocaTE2"
echo "Somatic"
bedtools intersect -a ../mPing_gff/Somatic.gff -b ../mPing_gff/MSU_r7.all.final.full.utr.gff3 -wao > Somatic.intersect
bedtools closest -a ../mPing_gff/Somatic.gff -b ../mPing_gff/MSU_r7.all.final.mRNA.gff -d > Somatic.mRNA.intersect
python mPing_position.py --input Somatic.intersect --mrna Somatic.mRNA.intersect
python mPing_intron.py --input Somatic.intersect
python mPing_intergenic.py --input Somatic.mRNA.intersect


echo "RIL"
bedtools intersect -a ../mPing_gff/RIL.gff -b ../mPing_gff/MSU_r7.all.final.full.utr.gff3 -wao > RIL.intersect
bedtools closest -a ../mPing_gff/RIL.gff -b ../mPing_gff/MSU_r7.all.final.mRNA.gff -d > RIL.mRNA.intersect
python mPing_position.py --input RIL.intersect --mrna RIL.mRNA.intersect
python mPing_intron.py --input RIL.intersect
python mPing_intergenic.py --input RIL.mRNA.intersect

echo "Strains"
bedtools intersect -a ../mPing_gff/Strains.gff -b ../mPing_gff/MSU_r7.all.final.full.utr.gff3 -wao > Strains.intersect
bedtools closest -a ../mPing_gff/Strains.gff -b ../mPing_gff/MSU_r7.all.final.mRNA.gff -d > Strains.mRNA.intersect
python mPing_position.py --input Strains.intersect --mrna Strains.mRNA.intersect
python mPing_intron.py --input Strains.intersect 
python mPing_intergenic.py --input Strains.mRNA.intersect 


echo "Simulation, use summary not one"
bedtools intersect -a ../mPing_gff/Simulate0001.gff -b ../mPing_gff/MSU_r7.all.final.full.utr.gff3 -wao > Simulation.intersect
bedtools closest -a ../mPing_gff/Simulate0001.gff -b ../mPing_gff/MSU_r7.all.final.mRNA.gff -d > Simulation.mRNA.intersect
python mPing_position.py --input Simulation.intersect --mrna Simulation.mRNA.intersect
python mPing_intron.py --input Simulation.intersect 
python mPing_intergenic.py --input Simulation.mRNA.intersect

echo "Landrace split, sofia version"
bedtools intersect -a ../mPing_gff/Landrace_RelocaTE2/HEG4.gff -b ../mPing_gff/MSU_r7.all.final.full.utr.gff3 -wao > HEG4.intersect
bedtools intersect -a ../mPing_gff/Landrace_RelocaTE2/EG4.gff -b ../mPing_gff/MSU_r7.all.final.full.utr.gff3 -wao > EG4.intersect
bedtools intersect -a ../mPing_gff/Landrace_RelocaTE2/A123.gff -b ../mPing_gff/MSU_r7.all.final.full.utr.gff3 -wao > A123.intersect
bedtools intersect -a ../mPing_gff/Landrace_RelocaTE2/A119.gff -b ../mPing_gff/MSU_r7.all.final.full.utr.gff3 -wao > A119.intersect

python mPing_position.py --input HEG4.intersect
python mPing_position.py --input EG4.intersect
python mPing_position.py --input A123.intersect
python mPing_position.py --input A119.intersect

python mPing_intron.py --input HEG4.intersect
python mPing_intron.py --input EG4.intersect
python mPing_intron.py --input A123.intersect
python mPing_intron.py --input A119.intersect

bedtools closest -a ../mPing_gff/Landrace_RelocaTE2/HEG4.gff -b ../mPing_gff/MSU_r7.all.final.mRNA.gff -d > HEG4.mRNA.intersect
bedtools closest -a ../mPing_gff/Landrace_RelocaTE2/EG4.gff -b ../mPing_gff/MSU_r7.all.final.mRNA.gff -d > EG4.mRNA.intersect
bedtools closest -a ../mPing_gff/Landrace_RelocaTE2/A123.gff -b ../mPing_gff/MSU_r7.all.final.mRNA.gff -d > A123.mRNA.intersect
bedtools closest -a ../mPing_gff/Landrace_RelocaTE2/A119.gff -b ../mPing_gff/MSU_r7.all.final.mRNA.gff -d > A119.mRNA.intersect

python mPing_intergenic.py --input HEG4.mRNA.intersect
python mPing_intergenic.py --input EG4.mRNA.intersect
python mPing_intergenic.py --input A123.mRNA.intersect
python mPing_intergenic.py --input A119.mRNA.intersect


