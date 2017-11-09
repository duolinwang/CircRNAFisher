#!/usr/bin/env perl

package find_match_splice_site;
#use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(find_match_splice_site);
our @EXPORT_OK = qw(find_match_splice_site);

sub find_match_splice_site
{
my ($fragment,$readlength,$flariat)=@_;
#my $windowsize=$fragment-$readlength;
#my $max_fragmentsize=$fragment-$readlength;
if($flariat==1){$in1="./output/re_spliceplot_withstrand_detetelariat.txt"}else{$in1="./output/re_spliceplot_withstand.txt";}
print "open $in1\n";
open(IN1,$in1)or die "Can't open file \"$in1\" : $!";
#open(IN1,"./output/merged_backsplice.txt")or die "Can't open file \"./output/merged_backsplice.txt\" : $!";
open(IN2,"./output/discordant.bed")or die "Can't open file \"./output/disorder.bed\" : $!";
open(OUT,">./output/pairedreads_match_splicesite_within_fragment")or die "$!";
open(OUT2,">./output/valid_disorder.bed")or die "$!";
open(OUT3,">./output/unmatched_disorder.bed")or die "$!";
while($line1=<IN1>)
{
@a=split(/\t/,$line1);
@b=split(/,/,$a[0]);#chr pos1 pos2
if(!exists($splice_pos1{$b[0]}))
{
$splice_pos1{$b[0]}[0]="$b[1]";
$splice_pos2{$b[0]}[0]="$b[2]";
#print "$a[1]\n";
}
else
{
$num=scalar(@{$splice_pos2{$b[0]}});
$splice_pos1{$b[0]}[$num]="$b[1]";
$splice_pos2{$b[0]}[$num]="$b[2]";
}
}
while($line2=<IN2>)
{
$find=0;# flag indicate whether this disorder paired has one junction.
########read1##############
@a=split(/\t/,$line2);
$read1=$a[3];
$chr=$a[0];
$posread1s=$a[1];
$posread1d=$a[2];
$read1pos="$a[0],$a[1],$a[2]";
#print "$readlength\n";
$strand=$a[5];
chomp($strand);
$read1line=$line2;
###########read2############
$line2=<IN2>;
@a=split(/\t/,$line2);
$read2=$a[3];
$posread2s=$a[1];
$posread2d=$a[2];
$read2pos="$a[0],$a[1],$a[2]";
$read2line=$line2;
############################
$splice_num1=0;
$splices1="";
$splice_num2=0;
$splices2="";
$time=0;
 if(exists($splice_pos1{$chr}))
 {
  $num=scalar(@{$splice_pos1{$chr}});
  if($strand eq "+")
  {
    $start=$posread2s;
    $end=$posread1d;
  }
  else#strand = -
  {
    $start=$posread1s;
    $end=$posread2d;  
  } 
  #print "$read1pos $read2pos $start   $end\n";  
  for($i=0;$i<$num;$i++)
  {
     $d_start=$start-$splice_pos1{$chr}[$i];
     $d_end=$splice_pos2{$chr}[$i]-$end;
     if($d_start>=0 and $d_end>=0)
     {
        
	#print "$d_start-$d_end\n";
        if(($d_start+$d_end)<=($fragment-2*$readlength))
        {
	 if($time==0)
         {
         $time=1;
	 $find=1;
         print OUT2 $read1line;
         print OUT2 $read2line;
         }	
         $splice_num1+=1;
         if($splice_num1==1)
         {
          $splices1=$chr.",".$splice_pos1{$chr}[$i].",".$splice_pos2{$chr}[$i];
         }
         else{$splices1=$splices1.":".$chr.",".$splice_pos1{$chr}[$i].",".$splice_pos2{$chr}[$i];}
        }
     }
  }
 }#exists chr
 if($splice_num1>0)
 {
   print OUT "$read1\t$splice_num1\t$splices1\n";
   print OUT "$read2\t$splice_num1\t$splices1\n";
 }
 if($find==0)
 {
  print OUT3 $read1line;
  print OUT3 $read2line;

 }
}
close IN1;
close IN2;
close OUT2;
close OUT;
close OUT3;
}
1;
