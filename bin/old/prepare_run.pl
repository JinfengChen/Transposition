#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"input:s","help");


my $help=<<USAGE;
prepare_run.pl --input RIL.fastq.list

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}


my $list = readtable($opt{input});
my @rils = sort keys %$list;
my $lines = 20; # split by 10 rils in each bash
my %hash;
for(my $i=0; $i<@rils; $i++){
   my $index = int($i/$lines);
   #print "$index\t$rils[$i]\t$list->{$rils[$i]}\n";
   my $cmd1 = "perl /rhome/cjinfeng/software/tools/RelocaTE-master/scripts/relocaTE.pl -t ../input/TE/Rice.TE.fa -g /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa -d $list->{$rils[$i]} -e $rils[$i] -o ./RelocaTE/$rils[$i] -r ./all_reference_repeat.txt -p 1 -a 0";
   my $cmd2 = "python runRelocaTEjobs.py --input ./RelocaTE/$rils[$i]/run_these_jobs.sh";
   push @{$hash{$index}}, $cmd1;
   push @{$hash{$index}}, $cmd2;
}

foreach my $index (sort {$a <=> $b} keys %hash){
   #print "$index\n";
   open OUT, ">scripts_$index.sh" or die "$!";
   my $cmd = join("\n", @{$hash{$index}});
   print OUT "$cmd\n";
   close OUT;
}

#perl /rhome/cjinfeng/software/tools/RelocaTE-master/scripts/relocaTE.pl -t ../input/TE/Rice.TE.fa -g /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/ref/MSU_r7.fa -d ../input/fq/GN100 -e GN100 -o ./RelocaTE/GN100 -r ./all_reference_repeat.txt -p 1 -a 0
#python runRelocaTEjobs.py --input ./RelocaTE/GN100/run_these_jobs.sh


#GN100	/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/fastq/GN100
sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}
 
