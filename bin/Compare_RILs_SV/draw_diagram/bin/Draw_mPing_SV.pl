#!/usr/bin/perl
use strict;
use FindBin qw ($Bin);
use Getopt::Long;
use SVG;

my %opt;
GetOptions(\%opt,"table:s","help:s");

my $help=<<USAGE;
perl $0 --table ../input/QTL.regions.list
USAGE

if(defined $opt{help} or keys %opt < 1){
        die  $help ;
}



my $svg=SVG->new(width=>1200,height=>800);
my $startw=50; 
my $endw  =1100;
my $starth=100;
my $endh  =500;


plot_deletion_two_mping($svg, $startw, $starth, 'Type I: Deletions between two nearby mPing insertions', 't.txt');
#plot_deletion_one_mping($svg, $startw+500, $starth, 'Type 2: Deletions from one mPing insertions to flanking sequence');
#plot_deletion_one_mping($svg, $startw, $starth+350, 'Type 3: Deletions accompanied by rearrangement of flanking sequence');
#plot_deletion_one_mping()
#plot_deletion_rearrangement()


writesvg("mPing_SV_diagram.svg", $svg);

###sub functions
sub plot_deletion_two_mping
{
my ($svg, $startx, $starty, $title, $plot_inf) = @_;
my $plot_data = read_inf($plot_inf);
plot_text($svg, $startx+15, $starty-50, 'start', 15, $title);

foreach my $rank (keys %$plot_data){
    my $x = $startx;
    my $y = $starty + 50 + ($rank-1)*150;
    plot_diagram_simple($svg, $plot_data->{$rank}, $x, $y);
}
}


#sub plot_deletion_one_mping
#{
#my ($svg, $startx, $starty, $title) = @_;
#
#my $subtitle =$svg->text(
#              x=>$startx+15, y=>$starty-20,
#              style=>{
#                   fontsize=>'4','text-anchor'=>'start','stroke-width'=>1
#              }
#    )->cdata($title);

#}

################################### sub for plot diagram of one SV#####################
sub plot_diagram_simple
{
my ($svg, $data, $x, $y) = @_;

#upper plot: HEG4
plot_text($svg, $x, $y, 'start', '15', $data->{'up_header'});
#plot_chromosome_line()
plot_zoom_in_mping($svg, $x, $y, 'up', $data) if exists $data->{'up_zoom_in_mping'};
plot_zoom_in($svg, $x, $y, 'up', $data) if exists $data->{'up_zoom_in'};

#lower plot: RILs
plot_text($svg, $x, $y+70, 'start', '15', $data->{'down_header'});
#plot_chromosome_line() 
plot_zoom_in_mping($svg, $x, $y+70, 'down', $data) if exists $data->{'down_zoom_in_mping'};
plot_zoom_in($svg, $x, $y+70, 'down', $data) if exists $data->{'down_zoom_in'}; 

}

################################### sub for plot_strain name ##########################
sub plot_text
{
my ($svg, $x, $y, $anchor, $size, $strain) =@_;
my $strain_name=$svg->text(
               x=>$x, y=>$y,
               style=>{
                    'font-size'=>$size,'text-anchor'=>$anchor,'stroke-width'=>1
               }
    )->cdata($strain);
}

################################### sub for plot zoom in mping #######################
sub plot_zoom_in_mping
{
my ($svg, $x, $y, $pos, $data) =@_;
my $scale  = 500/5000;
my $start = 0;
my $end   = 0;
my $upcut   = 0;
my $upcut_len= 0;
my $downcut = 0;
my $downcut_len = 0;
my $zoom_in_y = $y;
my @zoom_in;
if ($pos =~ /up/){
   $start = $data->{'up_start'};
   $end   = $data->{'up_end'};
   if (exists $data->{'up_cut'}){
       my @temp = split(',', $data->{'up_cut'});
       $upcut = $temp[0]; 
       $upcut_len = $temp[1];
   }
   $zoom_in_y = $zoom_in_y - 45;
   @zoom_in = split(";", $data->{'up_zoom_in_mping'});
}elsif($pos =~ /down/){
   $start = $data->{'down_start'};
   $end   = $data->{'down_end'};
   if (exists $data->{'down_cut'}){
       my @temp = split(',', $data->{'down_cut'});
       $downcut = $temp[0];
       $downcut_len = $temp[1];
   }
   $zoom_in_y = $zoom_in_y - 45;
   @zoom_in = split(";", $data->{'down_zoom_in_mping'});
}

for (my $i=0; $i<@zoom_in; $i++){
   my @zoom_in_inf = split(',', $zoom_in[$i]);
   my $zoom_in_x = ($zoom_in_inf[0] - $start + 1)*$scale + $x;
   if ($pos =~ /up/ and $zoom_in_inf[0] > $upcut and $upcut > 0){
       $zoom_in_x = ($zoom_in_inf[0] - $start - $upcut_len + 1)*$scale + $x;
   }elsif($pos =~ /down/ and $zoom_in_inf[0] > $downcut and $downcut > 0){
       $zoom_in_x = ($zoom_in_inf[0] - $start - $downcut_len + 1)*$scale + $x;
   }
   my $zoom_in_left_seq  = $zoom_in_inf[2];
   my $zoom_in_right_seq = $zoom_in_inf[3];
   my $zoom_in_strand    = $zoom_in_inf[1];
   print "$zoom_in_inf[0]\t$start\t$x\t$zoom_in_x\t$zoom_in_y\n";
   plot_text($svg, $zoom_in_x-50, $zoom_in_y, 'middle', '10', $zoom_in_left_seq);
   plot_text($svg, $zoom_in_x+50, $zoom_in_y, 'middle', '10', $zoom_in_right_seq);
   plot_text($svg, $zoom_in_x, $zoom_in_y, 'middle', '10', 'mPing');
   plot_zoom_in_lines($svg, $zoom_in_x, $zoom_in_y)
}

}

################################### sub for plot_zoom in breakpoint signiature #######
sub plot_zoom_in
{
my ($svg, $x, $y, $pos, $data) =@_;
my $scale  = 500/5000;
my $start = 0;
my $end   = 0;
my $upcut   = 0;
my $upcut_len= 0;
my $downcut = 0;
my $downcut_len = 0;
my $zoom_in_y = $y;
my @zoom_in;
if ($pos =~ /up/){
   $start = $data->{'up_start'};
   $end   = $data->{'up_end'};
   if (exists $data->{'up_cut'}){
       my @temp = split(',', $data->{'up_cut'});
       $upcut = $temp[0];
       $upcut_len = $temp[1];
   }

   $zoom_in_y = $zoom_in_y - 45;
   @zoom_in = split(";", $data->{'up_zoom_in'});
}elsif($pos =~ /down/){
   $start = $data->{'down_start'};
   $end   = $data->{'down_end'};
   if (exists $data->{'down_cut'}){ 
       my @temp = split(',', $data->{'down_cut'});
       $downcut = $temp[0];
       $downcut_len = $temp[1];
   }
   $zoom_in_y = $zoom_in_y - 45;
   @zoom_in = split(";", $data->{'down_zoom_in'});
}
for (my $i=0; $i<@zoom_in; $i++){
   #print $zoom_in[$i],"\n";
   my @zoom_in_inf = split(',', $zoom_in[$i]);
   #print $zoom_in_inf[0],"\n";
   my $zoom_in_x = ($zoom_in_inf[0] - $start + 1)*$scale + $x;
   if ($pos =~ /up/ and $zoom_in_inf[0] > $upcut and $upcut > 0){
       $zoom_in_x = ($zoom_in_inf[0] - $start - $upcut_len + 1)*$scale + $x;
   }elsif($pos =~ /down/ and $zoom_in_inf[0] > $downcut and $downcut > 0){
       $zoom_in_x = ($zoom_in_inf[0] - $start - $downcut_len + 1)*$scale + $x; 
   }

   my $zoom_in_left_seq = $zoom_in_inf[1];
   my $zoom_in_right_seq = $zoom_in_inf[3];
   my $zoom_in_filler_seq = $zoom_in_inf[2];
   my $zoom_in_string = $zoom_in_left_seq.' '.$zoom_in_filler_seq.' '.$zoom_in_right_seq;
   print "$zoom_in_inf[0]\t$start\t$x\t$zoom_in_x\t$zoom_in_y\t$zoom_in_string\n"; 
   plot_text($svg, $zoom_in_x, $zoom_in_y, 'middle', '10', $zoom_in_string);
   plot_zoom_in_lines($svg, $zoom_in_x, $zoom_in_y)
} 
}

################################### sub for plot zoom in lines
sub plot_zoom_in_lines
{
my ($svg, $x, $y) = @_;
my $top_left = $x-80;
my $top_right= $x+80;
my $middle_y= $y + 30;
my $bottom_y= $y + 45; 
my $baseline=$svg->line(
     x1=>$x,y1=>$middle_y,
     x2=>$x,y2=>$bottom_y,
     style=>{stroke=>'black'}
   );
my $lefline=$svg->line(
     x1=>$x,y1=>$middle_y,
     x2=>$top_left,y2=>$y,
     style=>{stroke=>'black'}
   );   
my $rightline=$svg->line(
     x1=>$x,y1=>$middle_y,
     x2=>$top_right,y2=>$y,
     style=>{stroke=>'black'}
   );   

}

################################### sub for read ploting infromation###################
sub read_inf
{
my ($file) = @_;
my %hash;
my $rank=0;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_ =~/^#/);
    if ($_=~/^>(\d+)/){
        $rank = $1;
        next;
    }else{
        my @unit=split("\=", $_);
        #print $unit[0], "\t", $unit[1], "\n";
        $hash{$rank}{$unit[0]}=$unit[1];
    }
}
close IN;
return \%hash;
}


################################### sub for write svg to file##########################
sub writesvg {
my ($file, $svg)=@_;
open OUT, ">$file" or die "can not open my file";
    print OUT $svg->xmlify;
close OUT;
system("/rhome/cjinfeng/BigData/software/draw/svg2xxx_release/svg2xxx $file -t pdf");
}
