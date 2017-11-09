#!/usr/bin/env perl

package get_unique_reads;

#use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(get_unique_reads);
our @EXPORT_OK = qw(get_unique_reads);

sub get_unique_reads{
my ($reflag)=@_;
print "Select uniquely mapped reads...\n";
######test the read count in order_splice_read and also in paper read#################
open(IN1,"./output/".$reflag."pass_paired_end_read.txt")or die "can't open file \"./output/".$reflag."pass_paired_end_read.txt\": $!";
open(OUT1,">./output/".$reflag."pass_paired_end_read_unique.txt");
open(OUT2,">./output/".$reflag."pass_paired_end_readfilter1");
#####################################################################################
my $all_read=0;
my $rc=0;
my $unique=0;
my $count=0;
while($line1=<IN1>)
{
  chomp($line1);
  @a=split(/\t/,$line1); 
  if($a[1]==1){print OUT1 "$line1\n";$unique+=1;}
  if($a[1]>1){print OUT2 "$line1\n";$count+=1;}
  $all_read+=1;
} 
print "unique read num=$unique\tmultiple read num=$count\n";

close IN1;
close OUT1;
close OUT2;
}

1;

