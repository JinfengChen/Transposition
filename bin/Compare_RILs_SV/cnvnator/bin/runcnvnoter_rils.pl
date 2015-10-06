#!/usr/bin/perl
use Getopt::Long;
use File::Spec;

GetOptions (\%opt,"bam:s","project:s","help");


my $help=<<USAGE;
perl $0 --bam
Prepare shell for pindel run on biocluster using pbs. After this script, run "qsub -q js pindel.sh".
--bam: dir of bam file to run, bam should be indexed using samtools
Note that BP is breakpoints file for which one breakpoints is not known and LI is long insertion which two breakpoint are determined but the length and insertion sequence was not known. So these kind of thing can not be usd in drawing statistic figure, but need to investigate in detail.
USAGE


if ($opt{help} or keys %opt < 0){
    print "$help\n";
    exit();
}

$opt{project} ||= "HEG4";
$opt{project} = File::Spec->rel2abs($opt{project});
$opt{bam} ||= "../input";
config() unless (-e "$opt{project}.pindel.config");
#`perl /rhome/cjinfeng/software/bin/qsub-pbs.pl --maxjob 10 --lines 3 --interval 120 --resource nodes=1:ppn=12,walltime=100:00:00,mem=12g --convert no $opt{project}.pindel.sh`;


sub shell
{
my ($bam, $ril)=@_;
my $prefix = $1 if ($bam=~/(.*)\.bam/);
###-c ALLmeans to run all chromosome
###-T 12 means using 12 cpus
###-x 5  means detect SV < 32,000 bp,6=129,472, see pindel options
###-v 50 means report inversion > 50 bp, see pindel options
my $Shell=<<CMD;
#!/bin/sh
#PBS -l nodes=1:ppn=12
#PBS -l mem=10gb
#PBS -l walltime=100:00:00
#PBS -d ./

date

module load samtools
PATH=\$PATH:~/BigData/software/SVcaller/ROOT/bin/

readn=`head -n 1 $prefix.flagstat | cut -d" " -f1`
if [ ! -e $ril.readdepth.bed ] && [ \$readn -gt 18600000 ]; then

   samtools view -h $bam | sed 's/Chr//g' | samtools view -Sb - -o $ril.bam
   samtools index $ril.bam 
   /usr/bin/python2.7 /bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq//bin/cnvnator_wrapper.py --cnvnator /bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq//bin/cnvnator-multi -T $ril-cnvnator-temp -t \$PBS_NP -w 100 -b $ril.bam -o $ril.readdepth -c /bigdata/stajichlab/cjinfeng/software/SVcaller/speedseq/annotations/cnvnator_chroms
   rm $ril.bam $ril.bam.bai
else
   echo "Read coverage low than 5X: \$readn"
fi

date
echo "Done"

CMD


open OUT, ">$opt{project}.$ril.cnvnator.sh" or die "$!";
     print OUT "$Shell\n";
close OUT;
}


sub config
{
my @bam=glob("$opt{bam}/*.bam");
#open SH, ">$opt{project}.cnvnoter.sh" or die "$!";
for (my $i=0;$i<@bam;$i++){
    #../input/RILs_ALL_bam/GN1.bam
    my $ril = $1 if ($bam[$i]=~/(GN\d+)\.bam/); 
    #open OUT, ">$opt{project}.$ril.pindel.config" or die "$!";
    #print OUT "$bam[$i]\t250\t$opt{project}.$ril\_250\n";
    #close OUT;
    #print SH "/opt/linux/centos/7.x/x86_64/pkgs/pindel/0.2.5a7/bin/pindel -T 12 -c ALL -f /rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa -o $opt{project}.$ril.pindel -i $opt{project}.$ril.pindel.config -x 6 -v 50\n";
    #print SH "/rhome/cjinfeng/BigData/software/SVcaller/pindel/pindel2vcf -P $opt{project}.$ril.pindel -r /rhome/cjinfeng/HEG4_cjinfeng/seqlib/MSU_r7.fa -R MSU7 -d 2013 -v $opt{project}.$ril.pindel.vcf\n";
    #print SH "/rhome/cjinfeng/HEG4_cjinfeng/Variations/SV/pindel/bin/pindel2GFF.sh $opt{project}.$ril.pindel $opt{project}.$ril.pindel.gff\n";
    shell($bam[$i], $ril)
}
#close SH;
} 


