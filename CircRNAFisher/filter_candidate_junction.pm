#!/usr/bin/env perl

package filter_candidate_junction;

#use strict;
use List::Util qw/sum min/; 

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(filter_candidate_junction);
our @EXPORT_OK = qw(filter_candidate_junction);

sub filter_candidate_junction
{
my ($FDR)=@_;
open(IN0,"./output/re_my_readcount_id_mismatch")or die "Can't open file \"./output/re_my_readcount_id_mismatch : $! \"";
open(IN,"./output/re_update_record_mate_nullscore")or die "Can't open file \"./output/re_update_record_mate_nullscore\" : $!";
open(OUT,">./output/passed_junction.txt")or die "$!";
open(OUT2,">./output/passed_reads_junctions.txt")or die "$!"; #need to support corresponding junctions

while($line=<IN0>)
{
chomp($line);
@a=split(/\t/,$line);
@b=split(/:/,$a[2]);
$read1score{$a[0]}=0;
for($i=0;$i<=$#b;$i++)
{
 @c=split(/,/,$b[$i]);
 $read1score{$a[0]}+=$c[4];
}
$read1num{$a[0]}=scalar(@b);
}

while($line=<IN>)
{
chomp($line);
@a=split(/\t/,$line);
@b=split(/:/,$a[2]);
@read2score=split(/,/,$a[1]);
$read2num=scalar(@read2score);
@back=();
@real1score=();
@real2score=();
@real2num=();
for($i=0;$i<=$#b;$i++)
{
 @c=split(/,/,$b[$i]);
 $back[$i]=join(",",@c[0..3]);
 $real1score[$i]=$c[4];
 $real2num[$i]=scalar(@c)-5;
 if($real2num[$i]==1){$real2score[$i]=$c[5];}
 else{$real2score[$i]=sum(@c[5..$#c])/$real2num[$i];}
}

$passbacknum=0;
@passback=();
for($i=0;$i<=$#b;$i++)
{
 if(($real2num[$i] >= $read2num) and ($#b==0)){$add=0;}
 else
 {
  $real=$real1score[$i]+$real2score[$i];#num=1
  #print "back=$back[$i]\n";
  #$t=scalar(@read1score);
  #print "backnum=$t\n";
  #$t=scalar(@read2score);  
  #print "matenum=$t\n";
  $null=($read1score{$a[0]}*$read2num+sum(@read2score)*$read1num{$a[0]})/($read1num{$a[0]}*$read2num);
  if($null==0){$add=1;}
  else{$add=min($real/$null,1);}
 }
 if($add<=$FDR)
 {
 $passback[$passbacknum]=$back[$i];$passbacknum+=1;
 }
 if(!exists($junction_num{$back[$i]}))
 {
  $junction_num{$back[$i]}=1;
  $junction_score{$back[$i]}[0]=$add;
 }
 else
 {
  $junction_num{$back[$i]}+=1;
  $junction_score{$back[$i]}[$junction_num{$back[$i]}-1]=$add;
 }
}#for
if($passbacknum>0)
{
 print OUT2 "$a[0]\t$passbacknum\t";
 for($j=0;$j<$passbacknum-1;$j++)
 {
   print OUT2 "$passback[$j]:";
 }
 print OUT2 "$passback[$j]\n";
}

}#while

foreach $key (keys %junction_num)
{
 $mean=0;
 $N=scalar(@{$junction_score{$key}});
 $sum=sum(@{$junction_score{$key}});
 if($sum!=0 and $N !=1)
 {
    for($i=0;$i<$N;$i++)
    { 
     $mean+=(1-$junction_score{$key}[$i]/$sum)*$junction_score{$key}[$i]/($N-1);#sum((1-xi/sum)/(n-1)*xi)
    }
 }
 elsif($sum == 0)
 {
  $mean=0;
 }
 else
 {
  $mean=$junction_score{$key}[0];
 }
 if($mean<=$FDR)
 {
 print OUT "$key\t$mean\n";
 }
 
}

close IN;
close OUT;
close OUT2;
}

1;
