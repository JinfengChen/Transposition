#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"te:s","gene:s","help");


my $help=<<USAGE;
perl AnnoTEsite.pl --te gff/HEG4.mping.all_inserts.gff --gene /rhome/cjinfeng/HEG4_cjinfeng/seqlib/GFF/MSU7.gene.gff
perl AnnoTEsite.pl --te gff/HEG4.mping.all_inserts.gff,gff/EG4.mping.all_inserts.gff,gff/A119.mping.all_inserts.gff,gff/A123.mping.all_inserts.gff --gene /rhome/cjinfeng/HEG4_cjinfeng/seqlib/GFF/MSU7.gene.gff
perl AnnoTEsite.pl --te gff/NB.mping.all_inserts.Non_reference.gff,gff/HEG4.mping.all_inserts.Non_reference.gff,gff/EG4.mping.all_inserts.Non_reference.gff,gff/A119.mping.all_inserts.Non_reference.gff,gff/A123.mping.all_inserts.Non_reference.gff --gene /rhome/cjinfeng/HEG4_cjinfeng/seqlib/GFF/MSU7.gene.gff

Annotate the insertion sites of transposons in genome. 
1. Compare the position of transposon to mRNA position, if -1500 <= dis < 0 then 5 UTR; if 0 < dis <= 1500 then 3 UTR; if > 1500 or < -1500 then intergenic; 
2. Compare the position of transposon to CDS position, when dis compare with mRNA == 0, use CDS dis. if dis == 0 then exon; else intron.
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $bedtools="/opt/bedtools/2.17.0-25-g7b42b3b/bin/bedtools";

my $geneprefix=basename($opt{gene},".gene.gff");
`grep "mRNA" $opt{gene} > $geneprefix.mRNA.gff` unless (-e "$geneprefix.mRNA.gff");
`grep "CDS" $opt{gene} > $geneprefix.CDS.gff` unless (-e "$geneprefix.CDS.gff");
my @te=split(",",$opt{te});

foreach(@te){
   #my $teprefix=basename($_,".mping.all_inserts.gff"); #
   my $teprefix=basename($_,".CombinedGFF.characterized.clean.unique_mPing.gff");
   `$bedtools closest -d -D b -a $_ -b  $geneprefix.mRNA.gff > $teprefix.mRNA.closest` unless (-e "$teprefix.mRNA.closest");
   `$bedtools closest -d -D b -a $_ -b  $geneprefix.CDS.gff > $teprefix.CDS.closest` unless (-e "$teprefix.CDS.closest");
   my $refmrna=closestBED("$teprefix.mRNA.closest");
   my $refcds =closestBED("$teprefix.CDS.closest");
   Annotation($teprefix,$refmrna,$refcds);
   `rm $teprefix.mRNA.closest $teprefix.CDS.closest`;
}


###########################################
sub Annotation
{
my ($prefix,$mrna,$cds)=@_;
my %hash;
foreach(keys %$mrna){
   if ($mrna->{$_}->[1] == 0){  ###in genic region
      if ($cds->{$_}->[1] == 0){###in exon
         $hash{$_}=[$mrna->{$_}->[2],$mrna->{$_}->[3],$mrna->{$_}->[4],"Exon",$mrna->{$_}->[0],$mrna->{$_}->[1]];
      }else{ ###in intron
         $hash{$_}=[$mrna->{$_}->[2],$mrna->{$_}->[3],$mrna->{$_}->[4],"Intron",$mrna->{$_}->[0],$mrna->{$_}->[1]];
      }
   }elsif($mrna->{$_}->[1] > 0 and $mrna->{$_}->[1] <= 3000){   ### 3 UTR
      $hash{$_}=[$mrna->{$_}->[2],$mrna->{$_}->[3],$mrna->{$_}->[4],"Three_primer_UTR",$mrna->{$_}->[0],$mrna->{$_}->[1]]; ### TE id->[annotation, gene, distance]
   }elsif($mrna->{$_}->[1] < 0 and $mrna->{$_}->[1] >= -3000){   ### 5 UTR
      $hash{$_}=[$mrna->{$_}->[2],$mrna->{$_}->[3],$mrna->{$_}->[4],"Five_primer_UTR",$mrna->{$_}->[0],$mrna->{$_}->[1]];
   }else{ ### intergenic
      $hash{$_}=[$mrna->{$_}->[2],$mrna->{$_}->[3],$mrna->{$_}->[4],"Intergenic",$mrna->{$_}->[0],$mrna->{$_}->[1]];
   }
}
open OUT, ">$prefix.mPing.annotation.unsort" or die "$!";
foreach(keys %hash){
   my $line=join("\t",@{$hash{$_}});
   print OUT "$_\t$line\n";
}
close OUT;
`sort -k2,2 -k3,3n $prefix.mPing.annotation.unsort > $prefix.mPing.annotation`;
`rm $prefix.mPing.annotation.unsort`;
}

sub closestBED
{
#### get the distance from TE to nearest gene
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    my $distance;
    my $gene;
    if ($unit[17]=~/ID=(.*?);/){
       $gene=$1;
    }
    my $id;
    if ($unit[8]=~/ID=(.*?);/){
       $id=$1;
    }
    $hash{$id}=[$gene,$unit[18],$unit[0],$unit[3],$unit[4]];
}
close IN;
return \%hash;
}



sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=1;
}
close IN;
return \%hash;
}
 
