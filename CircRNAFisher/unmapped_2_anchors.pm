#!/usr/bin/env perl

package unmapped_2_anchors;
use strict;

sub unmapped_2_anchors{
my $anchol=$_[0];
open(IN,"./output/unmapped_sample.sam")or die "can't open file \"./output/unmapped_sample.sam\" : $!";
open(OUT,">./output/sample_anchors.fastq");
while(my $line=<IN>)
{
  chomp($line);
  my @a=split(/\t/,"$line");
  my $headname=$a[0];
  my $seq=$a[9];
  if(length($seq)>2*$anchol)
  {
  my $qual=$a[10];
  my $seq_a=substr($seq,0,$anchol);
  my $seq_b=substr($seq,-$anchol,$anchol);
  my $qual_a=substr($qual,0,$anchol);
  my $qual_b=substr($qual,-$anchol,$anchol);
  print OUT "@"."$headname"."_A_"."$seq\n";
  print OUT "$seq_a\n";
  print OUT "+\n";
  print OUT "$qual_a\n";
  print OUT "@"."$headname"."_B\n";
  print OUT "$seq_b\n";
  print OUT "+\n";
  print OUT "$qual_b\n";
  }
  else{return 2;}
}#while
close IN;
close OUT;
return 1;
}

1;
