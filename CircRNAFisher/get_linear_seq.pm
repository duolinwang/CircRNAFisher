#!/usr/bin/env perl

package get_linear_seq;

use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(get_linear_seq);
our @EXPORT_OK = qw(get_linear_seq);

our %hg19=();
sub readhg19
{
my ($chrid,$refseq)=@_;

my $filename=$refseq.$chrid.".fa";

print "open $filename\n";

open(IN2,$filename)or die "can't open $filename\n";
my $l=<IN2>;###read first line#######
$hg19{$chrid}='';###global####3
while($l=<IN2>)
{
chomp($l);
$hg19{$chrid}.=$l;
}
close IN2;
}

sub get_linear_seq{
my ($refseq,$anchorl,$readlength)=@_;
my $infile="";

$infile="./output/mapped_sample_linear";
open(IN1,$infile)or die "Can't open file $infile : $!";
open(OUT,">./output/ref_linear_seq.fa")or die "$!";
my $line1;
my %linear=();
while($line1=<IN1>)
{
 chomp($line1);
 my @a=split(/\t/,$line1);
 my $chr=$a[2];
 my $pos1=$a[3];
 my $pos2=$a[4];
 if(!exists($linear{"$chr,$pos1,$pos2"}))
 {
  $linear{"$chr,$pos1,$pos2"}=1; 
  if(!exists($hg19{$chr})){&readhg19($chr,$refseq);}
  my $cur="";
  my $length=$readlength-$anchorl;
  $cur.=substr($hg19{$chr},$pos1-$length,$length);
  $cur.=substr($hg19{$chr},$pos2,$length);
  print OUT ">$chr,$pos1,$pos2\n$cur\n";
 }
}
close IN1;
close OUT;
%hg19=();
}

1;

