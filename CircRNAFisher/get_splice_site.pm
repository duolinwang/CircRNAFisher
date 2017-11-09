#!/usr/bin/env perl

package get_splice_site;

#use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(get_splice_site);
our @EXPORT_OK = qw(get_splice_site);

sub get_splice_site
{
my $in1;
$in1="./output/my_readcount_id_mismatch";
open(IN1,$in1)or die "Can't open file $in1 : $!";
open(OUT,">./output/spliceplot.txt")or die "$!";

while($line1=<IN1>)
{
chomp($line1);
@a=split(/\t/,$line1);
@b=split(/:/,$a[2]);
for($i=0;$i<=$#b;$i++)
{
 @c=split(",",$b[$i]);	
 $key=join(",",@c[1..3]);
 if(!exists($passsplit{$key}))
 {
 $passsplit{$key}=1;
 print OUT "$key\n";
 }
}
}
close IN1;
close OUT;
}

1;
