#!/usr/bin/env perl

package order_linear_read;
use strict;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(order_linear_read);
our @EXPORT_OK = qw(order_linear_read);

sub order_linear_read
{
open(IN,"./output/mapped_sample_linear_score")or die "can't open file \"./output/mapped_sample_linear_score\" : $!";
open(OUT,">./output/my_readcount_linear") or die "$!";
my $line;
my %mread=();
my %mread_splice_site=();
while($line=<IN>)
{
  chomp($line);
  my @a=split(/\t/,$line);
  my $splice_site=join(",",@a[1..4]);
  my $readid=$a[0];
  my $score=$a[7];
  if(!exists($mread{$readid}))
  {
     $mread{$readid}=1;
     $mread_splice_site{$readid}=$splice_site.",".$score;
  }
  else
   {
     $mread{$readid}+=1;
     $mread_splice_site{$readid}.=":".$splice_site.",".$score;
   }
 } 
foreach my $r (sort{$mread{$b}<=>$mread{$a}}keys %mread)
{
 print OUT "$r\t$mread{$r}\t$mread_splice_site{$r}\n";
}

close IN;
close OUT;
}
#&order_linear_read();

1;

