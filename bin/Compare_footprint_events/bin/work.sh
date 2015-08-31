echo "testing on highexcision, almost right, early events are not included in input file, blacklist remove cases also"
python footprint.py --input ../input/mping.excision.draw.highexcision > log 2>&1 &
python footprint_events.py --input ../input/mping.excision.draw.highexcision --blacklist ../input/Bam.Core.blacklist > log 2>&1 &
echo "run on 419 mPing with 891 excision"
python footprint_events.py --input ../input/mPing_boundary.linked_50Mb_debug2.mping_excision.list --blacklist ../input/Bam.Core.blacklist --output Excision_newpipe_version1 > log 2>&1 &

echo "distribution of excision events"
cat mping.excision_events.distr.R | R --slave
echo "test for significant: binom.test  p=0.0003 for 5, p=0.003 for 4"
cat mping.excision_events.binomial_test.R | R --slave

