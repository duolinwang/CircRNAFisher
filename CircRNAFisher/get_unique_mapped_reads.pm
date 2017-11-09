#!/usr/bin/env perl

package get_unique_mapped_reads;
use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(get_unique_mapped_reads);
our @EXPORT_OK = qw(get_unique_mapped_reads);

sub get_unique_mapped_reads
{
open(IN,"./output/mapped_sample_k20.sam")or die "Can't open file \"./output/mapped_sample_k20.sam\" : $!"; #sorted by id
open(OUT,">./output/uniquemapped_hg19.sam")or die "$!";
my ($line,$idr,%readid,%readnum,$id);
if($line=<IN>)
{
$line=~/(\S+_\d)/;
$idr=$1;
$readid{$idr}=$line;
$readnum{$idr}=1;
}
while($line=<IN>)
{
  $line=~/(\S+_\d)/;
  $id=$1;
  if(exists($readid{$id}))
  {
    $readnum{$id}+=1;
  }
  else#idr！=id
  {
   if($readnum{$idr}==1)
   { 
     print OUT $readid{$idr};
   }
   %readid=();
   %readnum=(); 
   $idr=$id;
   $readid{$idr}=$line;
   $readnum{$idr}=1;
  }
}

if($readnum{$idr}==1)
{ 
  print OUT $readid{$idr};
}

close IN;
close OUT;
}

1;
