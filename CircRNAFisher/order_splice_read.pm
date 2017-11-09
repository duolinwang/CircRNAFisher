#!/usr/bin/env perl

package order_splice_read;
use strict;
sub order_splice_read
{
my $reflag=$_[0];
open(IN,"./output/".$reflag."splice_sites_read_mismatch")or die "can't open file \"./output/".$reflag."splice_sites_read_mismatch\" : $!";
open(OUT,">./output/".$reflag."my_readcount_id_mismatch") or die "$!";
my $line;
my %mread=();
my %mread_splice_site=();
while($line=<IN>)
{
 if($line!~/rang/)
 {
  chomp($line);
  my @a=split(/\t/,$line);
  my $splice_site=$a[1];
  my @b=split(/:/,$a[3]);
  my $strand=$b[0];
  my $readid=$a[0];
  my $score=$b[7];
  if(!exists($mread{$readid}))
  {
     $mread{$readid}=1;
     $mread_splice_site{$readid}=$strand.",".$splice_site.",".$score;
  }
  else
   {
     $mread{$readid}+=1;
     $mread_splice_site{$readid}.=":".$strand.",".$splice_site.",".$score;
   }
 } 
}
foreach my $r (sort{$mread{$b}<=>$mread{$a}}keys %mread)
{
 print OUT "$r\t$mread{$r}\t$mread_splice_site{$r}\n";
}

close IN;
close OUT;
}

1;
