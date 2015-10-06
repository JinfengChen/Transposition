echo "High excision RILs check list"
python mPing_unknown_gt_RILs.py --input ../input/mping.excision.draw.highexcision.1 --matrix ../input/High_excision_csv_Ping_HEG4_mPing_only/ --output mping.excision.draw.highexcision
python get_excision_bam.py --input mping.excision.draw.highexcision.draw --output mping.excision.draw.highexcision.draw_bam
tar -zcvf mping.excision.draw.highexcision.draw_bam.tar.gz mping.excision.draw.highexcision.draw_bam/**

echo "High excision pairs RILs check list"
python mPing_unknown_gt_RILs.py --input ../input/mping.excision.draw.highexcision.pairs --matrix ../input/High_excision_csv_Ping --output mping.excision.draw.highexcision.pairs
python get_excision_bam.py --input mping.excision.draw.highexcision.pairs.draw --output mping.excision.draw.highexcision.pairs.draw_bam
tar -zcvf mping.excision.draw.highexcision.pairs.draw_bam.tar.gz mping.excision.draw.highexcision.pairs.draw_bam/**

