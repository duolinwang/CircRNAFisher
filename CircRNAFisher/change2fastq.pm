#!/usr/bin/env perl

package change2fastq;
use strict;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(change2fastq);
our @EXPORT_OK = qw(change2fastq);
sub change2fastq
{
print "changing unmapped reads into fastq format (changsam2fastq)...\n";	
open(IN,"./output/unmapped_sample.sam")or die "can't open file \"./output/unmapped_sample.sam\" : $!";
open(OUT,">./output/unmapped.fastq") or die "$!";
my $line;
while($line=<IN>)
{
chomp($line);
my @a=split(/\t/,$line);
print OUT "\@$a[0]\n";
print OUT "$a[9]\n";
print OUT "+\n";
print OUT "$a[10]\n";
}
close IN;
close OUT;
}

1;
