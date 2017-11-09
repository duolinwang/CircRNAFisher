#!/usr/bin/env perl

package change2bed;
use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(get_backsplice_bed changemapped2bed);
our @EXPORT_OK = qw(get_backsplice_bed changemapped2bed);

sub get_backsplice_bed
{
my ($flariat)=@_;
my $in1;
if($flariat==1){$in1="./output/re_spliceplot_withstrand_detetelariat.txt"}else{$in1="./output/re_spliceplot_withstand.txt";}
open(IN,$in1)or die "Can't open file \"$in1\" : $!";
print "open $in1\n";
#open(IN,"./output/merged_backsplice.txt") or die "Can't open file \"./output/merged_backsplice.txt\":$!";  
open(OUT,">./output/backsplice_bed.bed") or die "$!";
my $line;
while($line=<IN>)
{
 chomp($line);
 my @a=split(/\t/,$line);
 my @b=split(/;/,$a[2]);
 my @c=split(/,/,$a[0]);
 my $chr=$c[0];
 my $i;
 for($i=0;$i<=$#b;$i++)
 {
  $b[$i]=~/(\S+):([\+|-]):(\d+),(\d+)/;
  my $readid=$1;
  my $strand=$2;
  my $left=$3;
  my $right=$4;
  my $pos1s=$c[1];
  my $pos1e=$c[1]+$left;  
  my $pos2s=$c[2]-$right;  
  my $pos2e=$c[2];  
  print OUT "$chr\t$pos1s\t$pos1e\t$readid\t255\t$strand\n";
  print OUT "$chr\t$pos2s\t$pos2e\t$readid\t255\t$strand\n";
 }
}
close IN;
close OUT;
}

sub changemapped2bed{
my ($unique)=@_;
my $infile="";
if($unique ==1 ){$infile="./output/unique_merge_hg19_trans";}
else{$infile="./output/merge_hg19_trans";}
open(IN,$infile)or die "Can't open file $infile : $!";
open(OUT,">./output/mapped_hg19.bed") or die "$!";
my (@a,$readid,$strand,$start,$score,$end,$line);
while($line=<IN>)
{
  chomp($line);
  @a=split(/\t/,$line);
  if($a[2] =~ /^(chr\d+|chrX|chrY|chrM)$/)
  {
   $readid=$a[0];
   $strand=$a[1];
   $start=$a[3];
   $end=$a[4];
   $score=$a[5];
   print OUT "$a[2]\t$start\t$end\t$readid\t$score\t$strand\n"
 }
}
close IN;
close OUT;
}

1;
