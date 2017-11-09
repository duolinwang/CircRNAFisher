#!/usr/bin/env perl

package merge_backsplice;
use List::Util qw/max/;
#use strict;
sub merge_backsplice
{
my $range=$_[0];
open(IN,"./output/re_spliceplot_withstand.txt")or die "!!!";
open(OUT,">./output/merged_backsplice.txt");
open(TEST,">./output/record_cluster.txt");
print "Merging multiple backsplice...\n";
while($line=<IN>)
{
chomp($line);
@a=split(/\t/,$line);
@b=split(/,/,$a[0]);
if(!exists($back_num{$b[0]}))
{
$back_num{$b[0]}[0]=$a[1];#????????
$back_start{$b[0]}[0]=$b[1];
$back_end{$b[0]}[0]=$b[2];
$back{$b[0]}[0]=$a[0];
$read{$b[0]}[0]=$a[2];
}else
{
$num=scalar(@{$back_num{$b[0]}});
$back_num{$b[0]}[$num]=$a[1];
$back_start{$b[0]}[$num]=$b[1];
$back_end{$b[0]}[$num]=$b[2];
$back{$b[0]}[$num]=$a[0];
$read{$b[0]}[$num]=$a[2];
}
}

foreach $chr (keys %back_num)
{
 @{$state{$chr}}=(0) x scalar(@{$back_num{$chr}}); 
 #print "$state{$chr}[0]-$state{$chr}[1]\n";
 for($i=0;$i<scalar(@{$back_num{$chr}});$i++)
 {
   #if($i%10==0){print "$i\n";}		  	        	 
   if($state{$chr}[$i]!=1)
   {
	   # print "state $i=$state{$chr}[$i]\n";	         
     for($j=$i;$j<scalar(@{$back_num{$chr}});$j++)#including i
     {
	   #print "j=$j\n";
       $dist=max(abs($back_start{$chr}[$i]-$back_start{$chr}[$j]),abs($back_end{$chr}[$i]-$back_end{$chr}[$j]));
       #print "dist $i $j =$dist\n";
       if($dist<=$range)
       {
        if(!exists($cluster{$chr}{$i}))
        {
         $cluster{$chr}{$i}=$j;
        }else{$cluster{$chr}{$i}.=",".$j;}
         $state{$chr}[$j]=1;
        }
      } 
      #print "cluster=$cluster{$chr}{$i}\n";
   }#j 
 }#i
}

foreach $chr (keys %cluster)
{
 foreach $i (keys %{$cluster{$chr}})
 {
  @index=split(/,/,$cluster{$chr}{$i});
  $select=0;
  $selectnum=0;
  $cluster_num=0;
  $readlist=();
  for($j=0;$j<=$#index;$j++)
  {  
     if($back_num{$chr}[$index[$j]]>$selectnum)
     {
      $selectnum=$back_num{$chr}[$index[$j]];
      $select=$index[$j];
     }
     $cluster_num+=$back_num{$chr}[$index[$j]];
     if($j==0){
     $readlist=$read{$chr}[$index[$j]];     
     }else{
       $readlist.=";".$read{$chr}[$index[$j]];
     }
     print TEST "$back{$chr}[$index[$j]];";
  }
  print TEST "\n";
  print OUT "$back{$chr}[$select]\t$cluster_num\t$readlist\n";
 }
}
print "Merge over.\n";
}

1;

