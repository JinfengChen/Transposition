#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"input:s","help");


my $help=<<USAGE;

perl mPing_dist.pl --input HEG4_2.3.mPing.20X.mping.all_inserts.gff > mPing_dist.txt

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

readtable($opt{input});

sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    push @{$hash{$unit[0]}}, $unit[3]
}
close IN;

foreach my $c (keys %hash){
    my @pos = sort {$a <=> $b} @{$hash{$c}};
    for (my $i=0; $i<@pos-1; $i++){
        my $dist = $pos[$i+1] - $pos[$i];
        print "$c\t$pos[$i]\t$dist\n"; 
    }
}
}
