#!/usr/bin/env perl

package get_backsplice_seq;

use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(get_backsplice_seq);
our @EXPORT_OK = qw(get_backsplice_seq);

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

sub get_backsplice_seq{
my ($refseq,$anchorl,$readlength)=@_;
my $infile="";

#$infile="./output/spliceplot_withstand.txt";
$infile="./output/spliceplot.txt";
#$infile="./output/merged_backsplice.txt";
open(IN1,$infile)or die "Can't open file $infile : $!";
open(OUT,">./output/ref_backsplice_seq.fa")or die "$!";
my $line1;
while($line1=<IN1>)
{
 chomp($line1);
 my @a=split(/\t/,$line1);
 my @b=split(/,/,$a[0]);
 my $chr=$b[0];
 my $pos1=$b[1];
 my $pos2=$b[2];
 if(!exists($hg19{$chr})){&readhg19($chr,$refseq);}
 my $cur="";
 my $length=$readlength-$anchorl;
 $cur.=substr($hg19{$chr},$pos2-$length,$length);
 $cur.=substr($hg19{$chr},$pos1,$length);
 print OUT ">$a[0]\n$cur\n";
}

close IN1;
close OUT;
%hg19=();
}

1;
