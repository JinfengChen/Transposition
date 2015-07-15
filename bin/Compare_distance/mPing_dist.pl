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
my %strand;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    push @{$hash{$unit[0]}}, $unit[3];
    $strand{"$unit[0]\.$unit[3]"} = $unit[5];
}
close IN;

open OUT1, ">mPing_dist.txt" or die "$!";
open OUT2, ">mPing_dist.100kb.list" or die "$!";
foreach my $c (keys %hash){
    my @pos = sort {$a <=> $b} @{$hash{$c}};
    for (my $i=0; $i<@pos-1; $i++){
        my $dist = $pos[$i+1] - $pos[$i];
        print OUT1 "$c\t$pos[$i]\t$dist\n";
        if ($dist <= 100000){
            my $mping1 = "$c\.$pos[$i]";
            my $mping2 = "$c\.$pos[$i+1]";
            print OUT2 "$mping1\t$mping2\t$dist\t$strand{$mping1}\t$strand{$mping2}\n";
        } 
    }
}
close OUT1;
close OUT2;
}
