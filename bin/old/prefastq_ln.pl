#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
GetOptions (\%opt,"RIL:s","fastq:s","help");


my $help=<<USAGE;
perl $0 --RIL /rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/Illumina --fastq /rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/fastq 

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{RIL} ||= "/rhome/cjinfeng/HEG4_cjinfeng/RILs/QTL_pipe/input/fastq/Illumina";
$opt{fastq} ||= "/rhome/cjinfeng/BigData/00.RD/RILs/Transpostion/input/fastq";
my $prefix = read_bam_stat("./RIL.bam.unique.stat");

my @fq=glob("$opt{RIL}/RIL*/*.fq");
my %rils;
open OUT, ">RIL.fastq.list" or die "$!";
for (my $i=0;$i<@fq;$i++){
   #print "$fq[$i]\n";
   next unless -e "$fq[$i]";
   if ($fq[$i]=~/(.*\/(RIL(\d+)\_.*))_p1\.fq/){
      #print "$2\n";
      my $name = $3;
      my $index = $2;
      next unless $prefix->{$index};
      my $fq2 =$1."_p2.fq";
      my $dirname = "$opt{fastq}/"."GN".$name;
      `mkdir $dirname` unless -e "$dirname";
      my $file1="$dirname/"."GN".$name."_1".".fq";
      my $file2="$dirname/"."GN".$name."_2".".fq";
      print OUT "GN$name\t$dirname\n";
      unless (-e "$file1" and -e "$file2"){
         #print "$fq[$i]\t$file1\n";
         #print "$fq2\t$file2\n";
         `ln -s $fq[$i] $file1`;
         `ln -s $fq2 $file2`;
      }
   }
}
close OUT;

#RIL.bam.unique.stat
#../../input/fastq/Bam/RIL1_0_CGTACG_FC153L5.recal.bam
sub read_bam_stat
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_!~/^GN/);
    my @unit=split("\t",$_);
    my $prefix = basename($unit[8], '.recal.bam');
    #print "$prefix\n"; 
    $hash{$prefix}=1;
}
close IN;
return \%hash;
}

