#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);
use strict;

my %opt;

GetOptions (\%opt,"table1:s", "table2:s", "help");


my $help=<<USAGE;
perl $0 --table1 NB_HEG4.count --table2 NB_EG4.count

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $len = getfastalen("../input/MSU7.all.cds");
my $rpkm1 = readtable($opt{table1}, $len, 'HEG4');
my $rpkm2 = readtable($opt{table2}, $len, 'EG4');

my @gene = readgene("../input/high_excision.list");
barplot(\@gene, $rpkm1, $rpkm2);

#######
sub barplot{
my ($gene, $rpkm1, $rpkm2) = @_;
for(my $i=0; $i<@$gene; $i++){
   my ($nbavg,$nbse) = ci($rpkm1->{$gene->[$i]}->{'NB'});
   #print "$nbavg,$nbse\n";
   my ($heg4avg, $heg4se) = ci($rpkm1->{$gene->[$i]}->{'HEG4'});
   #print "$heg4avg,$heg4se\n";
   my ($eg4avg, $eg4se) = ci($rpkm2->{$gene->[$i]}->{'EG4'});
   #print "$eg4avg,$eg4se\n";
   my @data = ($gene->[$i], $nbavg, $heg4avg, $eg4avg, $nbse, $heg4se, $eg4se);
   #print "$nbavg, $heg4avg, $eg4avg, $nbse, $heg4se, $eg4se\n";
   plot(@data);
}
}

sub plot{
my (@data) = @_;
my $line = join("\n", @data);
open OUT, ">$data[0].rpkm" or die "$!";
   print OUT "$line\n";
close OUT;

my $cmd=<<R;
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

pvalue <- function(x1, y1, x2, y2, top, p){
     star <- '*'
     if (p > 0.05) {star <- 'n.s.'}
     if (p < 0.001){ star <- '**'}
     if (p < 0.0001){ star <- '***'}
     segments(x1,y1,x1,top)
     segments(x1,top,x2,top)
     segments(x2,top,x2,y2)
     #segments(x1-0.2,y1,x1+0.2,y1)
     #segments(x2-0.2,y2,x2+0.2,y2)
     xt <- min(x1,x2)+abs(x2-x1)/2
     yt <- top*1.1
     text(xt,yt,star, cex=1.5)
} 

pdf("$data[0].pdf", width=5, height=7)
#png("$data[0].png", width=500, height=700)
par(mar=c(5,5,4,2))
data =read.table("$data[0].rpkm", skip=1)
expr = data[,1][1:3]
std = data[,1][4:6]
barx <- barplot(expr, space=1, col=c("black"), border=F, ylim=c(0,(max(expr)+max(std))*1.5), cex.axis = 1.5, cex.name = 1.5, cex.lab = 1.5, axis.lty=1, ylab="Expression (RPKM)")
error.bar(barx, expr, std)
axis(1,c(0.2,max(barx)+0.6),line=0,labels=c("",""))
text(barx,rep(-max(expr)*0.25,3),offset=2, cex=1.5, labels=c("NB", "HEG4", "EG4"),srt=0,xpd=TRUE)
#legend("topright",c("HEG4","Nipponbare"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("blue","orange"))

p1 = 0.1
y11 = (data[,1][1] + data[,1][4])*1.1 
x11 = barx[1]
y12 = (data[,1][2] + data[,1][5])*1.1
x12 = barx[2]-0.1
top1 = max(y11, y12)*1.1
p2 = 0.1
y21 = (data[,1][2] + data[,1][5])*1.1
x21 = barx[2] + 0.1
y22 = (data[,1][3] + data[,1][6])*1.1
x22 = barx[3]
top2 = max(y21, y22)*1.1

pvalue(x11,y11,x12,y12,top1, p1)
pvalue(x21,y21,x22,y22,top2, p2)
dev.off()

R

open OUT, ">$data[0].R" or die "$!";
print OUT "$cmd\n";
close OUT;

system("cat $data[0].R | R --slave");
}



sub readrpkm
{
my ($file)=@_;
my %hash;
my %total;
open IN, "$file" or die "$!";
my $header = <IN>;
chomp $header;
my @samples=split("\t", $header);
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    for(my $i=1; $i<@unit; $i++){
        #print "$unit[0]\t$samples[$i-1]\t$unit[$i]\n";
        $hash{$unit[0]}{$samples[$i-1]} = $unit[$i];
    }
}
close IN;
return \%hash;
}


#######
#NB1     NB2     NB3     HEG41   HEG42   HEG43
#LOC_Os01g01010  2560    1247    1549    861     1612    949
#LOC_Os01g01019  5       18      14      6       9       4
sub readtable
{
my ($file, $len, $strain)=@_;
my %hash;
my %total;
open IN, "$file" or die "$!";
my $header = <IN>;
chomp $header;
my @samples=split("\t", $header);
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    for(my $i=1; $i<@unit; $i++){
        #print "$unit[0]\t$samples[$i-1]\t$unit[$i]\n";
        $hash{$unit[0]}{$samples[$i-1]} = $unit[$i];
        $total{$samples[$i-1]} += $unit[$i];
    }
}
close IN;

my %data;
my $prefix = basename($file, '.count');
open OUT, ">$prefix.rpkm" or die "$!";
print OUT join("\t", @samples),"\n";
foreach my $g (sort keys %hash){
    my @exps;
    push @exps, $g;
    for(my $i=0; $i<@samples; $i++){
        my $count = $hash{$g}{$samples[$i]};
        my $rpkm  = (1000000000*$count)/($total{$samples[$i]}*$len->{$g});
        #print "G:$g\tL:$len->{$g}\tC:$count\tS:$total{$samples[$i]}\t$rpkm\n";
        push @exps, $rpkm;
    }
    my $line = join("\t", @exps);
    print OUT "$line\n";
    my @nb = ($exps[1], $exps[2], $exps[3]);
    my @heg4 = ($exps[4], $exps[5], $exps[6]);
    $data{$g}{'NB'} = \@nb;
    $data{$g}{$strain} = \@heg4;
}
close OUT;

return \%data;
}

#Chr3	1269856	1271783	LOC_Os03g03070	OsMADS50
#Chr6	2234119	2239162	LOC_Os06g05060	OsELF3
sub readgene
{
my ($file)=@_;
my @gene;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    push @gene, $unit[0];
}
close IN;
return @gene;
}



sub getfastalen
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$1 if $temp1[0]=~/(.*?)\.\d+/;
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $hash{$head}= length $seq;
}
close IN;
$/="\n";
return \%hash;
}


sub ci
{
my ($num)=@_;
my $loop=0;
my $total=0;
my $add_square;
foreach  (@$num) {
        next if ($_ eq "NA");
        #my $temp=log($_);
        #print "$_\n";
        my $temp=$_;
        $total+=$temp;
        $add_square+=$temp*$temp;
        $loop++;
}


my $number=$loop;
return (0,0,0) if ($number < 2);
my $mean=$total/$number;
#print "M: $mean\n";
my $SD=sqrt( ($add_square-$total*$total/$number)/ ($number-1) );
my $se=1.96*$SD/sqrt($number);
return $mean, $se;
}

 
