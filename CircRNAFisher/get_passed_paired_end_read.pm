#!/usr/bin/env perl

package get_passed_paired_end_read;

#use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(get_passed_paired_end_read);
our @EXPORT_OK = qw(get_passed_paired_end_read);

sub get_passed_paired_end_read
{
my ($reflag)=@_;
if($reflag eq "re_")
{
open(IN1,"./output/passed_reads_junctions.txt")or die "Can't open file \"./output/re_record_mate_nullscore\" : $!";
open(IN2,"./output/passed_junction.txt")or die "Can't open file \"./output/re_junction_meanscore.txt\" : $! ";
open(OUT,">./output/re_pass_paired_end_read.txt")or die "$!";

while($line2=<IN2>)
{
chomp($line2);
@a=split(/\t/,$line2);
$back{$a[0]}=1;
}

while($line1=<IN1>)
{
	chomp($line1);
	@a=split(/\t/,$line1);
	@b=split(/:/,$a[2]);
	$backnum=0;
	@backprint=();
	for($i=0;$i<=$#b;$i++)
	{
 		@c=split(/,/,$b[$i]);	
 		$key=join(",",@c[0..3]);
 		if(exists($back{$key}))
 		{
  		$backprint[$backnum]=$key;	 
  		$backnum+=1;
 		}
	}

	if($backnum>0)
	{
		print OUT "$a[0]\t$backnum\t";
		for($i=0;$i<$backnum-1;$i++)
		{
 		print OUT "$backprint[$i]:"; 
		}
 		print OUT "$backprint[$i]\n";
	}
}

	close IN1;
	close IN2;
	close OUT;
}#if reflag
else
{
 open(IN1,"./output/record_mate_nullscore")or die "Can't open file \"./output/re_record_mate_nullscore\" : $!";
 open(OUT,">./output/pass_paired_end_read.txt")or die "$!";	
 while($line1=<IN1>)
 {
	chomp($line1);
	@a=split(/\t/,$line1);
	@b=split(/:/,$a[2]);
	$backnum=0;
	@backprint=();
	for($i=0;$i<=$#b;$i++)
	{
 		@c=split(/,/,$b[$i]);	
 		$key=join(",",@c[0..3]);
               	$backprint[$backnum]=$key;	 
  		$backnum+=1;
	}
		print OUT "$a[0]\t$backnum\t";
		for($i=0;$i<$backnum-1;$i++)
		{
 		print OUT "$backprint[$i]:"; 
		}
 		print OUT "$backprint[$i]\n";

 }
 close IN1;
 close OUT;
}

}



1;
