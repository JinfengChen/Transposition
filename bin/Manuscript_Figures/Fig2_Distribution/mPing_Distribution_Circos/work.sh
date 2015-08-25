echo "Based on sofia results"
python /rhome/cjinfeng/BigData/00.RD/Circos/bin/scripts/distri_data_pre_gff.py --head MSU7.circos.head --input OpenChromatin.gff.list
python /rhome/cjinfeng/BigData/00.RD/Circos/bin/scripts/distri_data_pre_gff.py --head MSU7.circos.head --input mPing.gff.list
source ~/.bashrc
python /rhome/cjinfeng/BigData/00.RD/Circos/bin/CircosConf.py --input circos.mPing.config --output mPing_intron.conf
python /rhome/cjinfeng/BigData/00.RD/Circos/bin/CircosConf.py --input circos.mPing_chromatin.config --output mPing_chromatin.conf

echo "RelocaTE2"
python /rhome/cjinfeng/BigData/00.RD/Circos/bin/scripts/distri_data_pre_gff.py --head MSU7.circos.head --input mPing.gff.RelocaTE2.list
source ~/.bashrc
python /rhome/cjinfeng/BigData/00.RD/Circos/bin/CircosConf.py --input circos.mPing.config --output mPing_intron.conf
python /rhome/cjinfeng/BigData/00.RD/Circos/bin/CircosConf.py --input circos.mPing_chromatin.config --output mPing_chromatin.conf

echo "test difference on chr4"
grep "Chr4" GFF.Somatic.histogram.txt | awk '$2>17800000' | cut -d" " -f4 | perl ~/software/bin/numberStat.pl
grep "Chr4" GFF.Somatic.histogram.txt | awk '$2<17800000' | cut -d" " -f4 | perl ~/software/bin/numberStat.pl
grep "Chr4" GFF.Strain.histogram.txt | awk '$2<17800000' | cut -d" " -f4 | perl ~/software/bin/numberStat.pl
grep "Chr4" GFF.Strain.histogram.txt | awk '$2>17800000' | cut -d" " -f4 | perl ~/software/bin/numberStat.pl
grep "Chr4" GFF.RIL.histogram.txt | awk '$2>17800000' | cut -d" " -f4 | perl ~/software/bin/numberStat.pl
grep "Chr4" GFF.RIL.histogram.txt | awk '$2<17800000' | cut -d" " -f4 | perl ~/software/bin/numberStat.pl
cat fisher.test.R | R --slave
#strains and RIL showed significant different but not with Somatic
#The sample size of strain are small we might just say overall similar?
