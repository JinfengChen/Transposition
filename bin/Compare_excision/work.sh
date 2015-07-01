echo "analyze exision in RILs, codes were modified from ~/BigData/00.RD/RILs/Figures/mPing_stability/bin"


echo "allele frequncy of mPing in RILs, count all insertions, including hom, het and somatic"
python AlleleFrq.py --input RIL275_RelocaTEi.CombinedGFF.characterized.gff > mping.ril.frquency

echo "analyze allele frequency of all insertion site in RILs"
export R_LIBS=$R_LIBS:"/rhome/cjinfeng/software/tools/R-2.15.3/library/"
#from this distribution, we know that these insertion with frequency from 0.1-0.5 may have abnormal frequency due to excision or other factor

cat mping.ril_breakY.R | R --slave
#singleton
awk '$6==1' mping.ril.frquency | wc -l
#shared in >1
awk '$6>1' mping.ril.frquency | wc -l
#all insertion site
wc -l mping.ril.frequency


echo "Analyze excision in RILs"
python Excision.py --input RIL275_RelocaTEi.CombinedGFF.characterized.gff > log 2> log2 &
qsub -q js runEx.sh

python Excision_Ref.py --input RIL275_RelocaTE.CombinedGFF.Shared.gff --mode Shared
python Excision_Ref.py --input RIL275_RelocaTE.CombinedGFF.Ref_only.gff --mode Ref
python Excision_Ref.py --input RIL275_RelocaTEi.CombinedGFF.characterized.gff --mode Non_Ref

python Excision_test.py --input RIL275_RelocaTEi.CombinedGFF.characterized.gff > log3 2> log4 &

echo "stat"
#207 of 555 mping have excision events
wc -l mping.excision.table
#374 of 905 excision events have footprint
grep ">" -v mping.excision.table.log | wc -l
grep ">" -v mping.excision.table.log | awk '$2==1' | wc -l

#14 excision event occured in more than 10 RILs
awk '$2>=10' mping.excision.non_ref.table | wc -l
awk '$2>=10' mping.excision.non_ref.table > mping.excision.non_ref.10more.table
#24 excision event occured in 5-10 RILs
awk '$2>=5 && $2<10' mping.excision.non_ref.table | wc -l
awk '$2>=5 && $2<10' mping.excision.non_ref.table > mping.excision.non_ref.5-10.table

echo "Draw allele frequency of insertions sites that have excision events in RILs. This is partial data of allele frequency of all insertion sites"
cat mping.excision.frq.R | R --slave

echo "Average number of excision events in each interval of allele frequency. This tell us that 0.1-0.4 frequency alleles are parital caused by excisions. How about other not caused by excision? Hotspot?"
echo "Or maybe the insertion site is around the boundary of recombination bin, so we do not know if there is excision or not"
python mping.excision.avg.py
cat mping.excision.avg.R | R --slave

echo "correlation of excision number with ping, unique mping, mping number in RILs"
python mping.excision.ril.py --input bamcheck.log
cat mping.excision.ril.R | R --slave


echo "test difference between non_ref and ref_pnly excision"
cat fisher.test.R | R --slave

echo "get sub bam for draw example"
python get_excision_bam.py --input mping.excision.draw.example

echo "excision mPing gff"
python Excision_GFF.py --input mping.excision.non_ref.table --gff HEG4.ALL.mping.non-ref.gff
perl AnnoTEsite.pl --te mping.excision.non_ref.gff --gene /rhome/cjinfeng/HEG4_cjinfeng/seqlib/GFF/MSU7.gene.gff &

