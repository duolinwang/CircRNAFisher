#!/usr/bin/env perl

package filter_candidate_junction_bymaq;

#use strict;
use List::Util qw/sum min/; 

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(filter_candidate_junction_bymaq);
our @EXPORT_OK = qw(filter_candidate_junction_bymaq);

sub filter_candidate_junction_bymaq
{
my ($maq_thred)=@_;
open(IN,"./output/record_realscore_nullmodel")or die "Can't open file \"./output/update_realscore_nullmodel : $! \"";
open(OUT,">./output/pass_maq_read_junctions.txt")or die "$!";
open(TEST,">./output/test_junctions_read_maq.txt")or die "$!";

while($line=<IN>)
{
chomp($line);
@a=split(/\t/,$line);
@real=split(/:/,$a[1]);
$line2=<IN>;
chomp($line2);
@null=split(/\t/,$line2);
$null_model=0;
for($i=0;$i<=$#null;$i++)
{
    @read=split(/,/,$null[$i]);	  
    $null_model+=10**(-($read[0]+$read[1])/10);
}
@real1score=();
@real2score=();
@real2num=();
$pass_junction_num=0;
%pass_junction_read=();
print OUT "$a[0]\t";
for($i=0;$i<=$#real;$i++)
{
 @c=split(/,/,$real[$i]);
 $back=join(",",@c[0..3]);
 $real1score=$c[4];
 $real2num=scalar(@c)-5;
 if($real2num==1){$real2score=$c[5];}
 else{$real2score=sum(@c[5..$#c])/$real2num;}
 $real_model=10**(-($real1score+$real2score)/10);
 $wp=1-$real_model/$null_model;
 if($wp==0){$maq=255;}
 else{$maq=-10*log($wp);}
 if($maq>255){$maq=255;}
 if(!exists($junction_num{$back}))
 {
  $junction_num{$back}=1;
  $junction_score{$back}[0]=$maq;
 }
 else
 {
  $junction_num{$back}+=1;
  $junction_score{$back}[$junction_num{$back}-1]=$maq;
 }
 if($maq>=$maq_thred){$pass_junction_num+=1;$pass_junction_read{$back}=sprintf("%0.2f",$maq);}
 print TEST "$back\t($line\t$line2)\t$maq\n";
}#for i
print OUT "$pass_junction_num\t";
foreach $key (keys %pass_junction_read)
{
 print OUT "$key,$pass_junction_read{$key}:";
}
print OUT "\n";
}#while
close IN;
close OUT;
close TEST;

%pass_junction=();
foreach $key (keys %junction_num)
{
 $mean=sum(@{$junction_score{$key}})/$junction_num{$key};
 if($mean>=$maq_thred){
 $pass_junction{$key}=(sprintf "%0.2f",$mean);
 }
}
%junction_score=();
%junction_num=();

open(IN,"./output/pass_maq_read_junctions.txt")or die "$!";
open(OUT,">./output/re_pass_paired_end_read.txt")or die "$!";

while($line=<IN>)
{
 	chomp($line);
	@a=split(/\t/,$line);
	@b=split(/:/,$a[2]);
	$backnum=0;
	@backprint=();
	for($i=0;$i<$a[1];$i++)
	{
	        @c=split(/,/,$b[$i]);
 		$key=join(",",@c[0..3]);
		$score=$c[4];
 		if(exists($pass_junction{$key}))
 		{
  		$backprint[$backnum]="$key,$score,$pass_junction{$key}";	 
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

close IN;
close OUT;
}

1;

