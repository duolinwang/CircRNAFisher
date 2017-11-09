#!/usr/bin/env perl

package unmapped_perl_utils;
use strict;
use print_parameter;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(add_num runcommand merge_unmapped_reads);
our @EXPORT_OK = qw(add_num runcommand merge_unmapped_reads);
our $readlength=0;
our $librarysize;

sub add_num
{
print "merging paired end reads...\n";
my ($read1, $read2) = @_;
open(IN1,$read1)or die "ERROR: can't open input file $read1 : $!";
open(IN2,$read2)or die "ERROR: can't open input file $read2 : $!";
open(OUT,">./output/merged_seq.fastq");
my $idnum=0;
my $line1;
my $line2;
my $id;
my $r_length=0;
while($line1=<IN1>)
{ 
  if($line1=~/@/)
  {$idnum+=1;
   $line2=<IN2>; 
   $id="\@SRR.".$idnum."_2\n";
   print OUT $id;
   $line2=<IN2>;
   $r_length=length($line2)-1;
   if($readlength<(length($line2)-1)){$readlength=length($line2)-1;}#record max readlength
   print OUT $line2;
   $line2=<IN2>;
   print OUT $line2;
   $line2=<IN2>;
   chomp($line2);
   if($r_length ne length($line2)){my $ln1=$idnum*4-2; my $ln2=$idnum*4; die "ERROR: the length of read sequence and read quality are different at line $ln1 and $ln2 ! \n";} 
   print OUT "$line2\n";
##################################### 
   $id="\@SRR.".$idnum."_1\n";
   print OUT $id;
   $line1=<IN1>;
   $r_length=length($line1)-1;   
   if($readlength<(length($line1)-1)){$readlength=length($line1)-1;}#record max readlength
   print OUT $line1;
   $line1=<IN1>;
   print OUT $line1;
   $line1=<IN1>;
   chomp($line1);
   if($r_length ne length($line1)){my $ln1=$idnum*4-2; my $ln2=$idnum*4; die "ERROR: the length of read sequence and read quality are different at line $ln1 and $ln2 ! \n";} 
   print OUT "$line1\n";
  } 
}
&print_parameter("readlength",$readlength);
&print_parameter("library_size",$idnum*2);# later change it into uniquely mapped reads!
close IN1;
close IN2;
close OUT;
return 1;
}

sub runcommand {
    print $_[0]."\n";
    my $state=system($_[0]);
    if ($state != 0) {
	print "failed to execute:$_[0]\n";
	exit(-1);
    }
}

sub merge_unmapped_reads {
my ($unmapped1,$unmapped2)=@_;
print "merge unmapped reads from $unmapped1 and $unmapped2...\n";
open(IN1,$unmapped1)or die "ERROR: can't open $unmapped1 : $!";#has been sorted by keys
open(IN2,$unmapped2)or die "ERROR: can't open $unmapped2 : $!";#has been sorted by keys
open(OUT,">./output/unmapped_sample.sam");
my $inpos=tell(IN2);
while(my $line1=<IN1>)
{
  my @a=split(/\t/,$line1);
  $a[0]=~s/_//;
  $a[0]=~/\S+\.(\d+)/;
  my $id1=$1;
  seek(IN2,$inpos,0);
  while(my $line2=<IN2>)
  {
   my @b=split(/\t/,$line2);
   $b[0]=~s/_//;
   $b[0]=~/\S+\.(\d+)/;
   my $id2=$1;
   if($id1<$id2){last;}
   if($id1>$id2){$inpos=tell(IN2);}
   if($id1 eq $id2)
   {
    print OUT $line1;
    $inpos=tell(IN2);#find and update $inpos
    last;
   }
  }
}#while
return 1;
close IN1;close IN2;close OUT;
}

1;
