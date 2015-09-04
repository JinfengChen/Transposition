echo "testing on highexcision, almost right, early events are not included in input file, blacklist remove cases also"
python footprint.py --input ../input/mping.excision.draw.highexcision > log 2>&1 &
python footprint_events.py --input ../input/mping.excision.draw.highexcision --blacklist ../input/Bam.Core.blacklist > log 2>&1 &
echo "run on 419 mPing with 891 excision"
python footprint_events.py --input ../input/mPing_boundary.linked_50Mb_debug2.mping_excision.list --blacklist ../input/Bam.Core.blacklist --output Excision_newpipe_version1 > log 2>&1 &

echo "distribution of excision events"
cat mping.excision_events.distr.R | R --slave
echo "test for significant: binom.test  p=0.0003 for 5, p=0.003 for 4"
cat mping.excision_events.binomial_test.R | R --slave


echo "summary excision number"
python sum_excision.py --input Excision_newpipe_version1.footprint.list.noPing.txt

echo "manual footprint"
#copy mPing footprint out into mping.txt then use print_fp.py print certain range of footprint. open in mac using sublime and snapshot by cmd+shift+4 
python print_fp.py --input t.txt --start 40 --end 80 > t.txt.1
python print_fp_helper.py --input Excision_newpipe_version1.footprint.high_and_pairs.txt

