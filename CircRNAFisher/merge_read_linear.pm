#!/usr/bin/env perl

package merge_read_linear;
use List::Util qw/max/;
#use strict;
sub merge_read_linear
{
my $range=$_[0];
open(IN,"./output/re_mapped_linear.sam")or die "Can't open file: ./output/re_mapped.sam $!";
open(OUT,">./output/merged_remapped_linear.txt");
print "Merging multiple junctions from extention process...\n";
my %back_readc=();
while($line=<IN>)
{
  @a=split(/\t/,$line);
  if(!exists($back_readc{$a[2]}))
  {$back_readc{$a[2]}=1;}
  else{$back_readc{$a[2]}+=1;}
}
seek(IN,0,0);
#readfirst
$line=<IN>;
@a=split(/\t/,$line);
%id_cur=();
$id_cur{$a[0]}=$line;
@read_bk=();
$backnum=0;
$read_bk[$backnum]=$line;
$backnum+=1;

while($line=<IN>)
{
 @a=split(/\t/,$line);
 if(exists($id_cur{$a[0]}))
 {
  $read_bk[$backnum]=$line;
  $backnum+=1;
 }
 else{
  if($backnum>1)
  {
     for($i=0;$i<$backnum;$i++)
     {
	     #print OUT "backnum=$backnum\n";       	     
       @al=split(/\t/,$read_bk[$i]);
       @pos=split(/,/,$al[2]);
       $read_bk[$i]=~/AS:i:(\S+)/;
       if(!exists($back_num{$pos[0]}))
       {
        $back_num{$pos[0]}[0]=$back_readc{$al[2]};
        $back_start{$pos[0]}[0]=$pos[1];
        $back_end{$pos[0]}[0]=$pos[2];
        $read{$pos[0]}[0]=$read_bk[$i];
        $score{$pos[0]}[0]=$1;
       }else
       {
        $num=scalar(@{$back_start{$pos[0]}});
        $back_num{$pos[0]}[$num]=$back_readc{$al[2]};
        $back_start{$pos[0]}[$num]=$pos[1];
        $back_end{$pos[0]}[$num]=$pos[2];
        $read{$pos[0]}[$num]=$read_bk[$i];
        $score{$pos[0]}[$num]=$1;
       }
      }
     foreach $chr (keys %back_num)
     {
     @{$state{$chr}}=(0) x scalar(@{$back_num{$chr}}); 
     for($i=0;$i<scalar(@{$back_num{$chr}});$i++)
     {
   	  	        	 
       if($state{$chr}[$i]!=1)
       {
        for($j=$i;$j<scalar(@{$back_num{$chr}});$j++)#including i
        {
        $dist=max(abs($back_start{$chr}[$i]-$back_start{$chr}[$j]),abs($back_end{$chr}[$i]-$back_end{$chr}[$j]));
        if($dist<=$range)
        {
	if(abs($back_start{$chr}[$i]-$back_start{$chr}[$j])==abs($back_end{$chr}[$i]-$back_end{$chr}[$j]))
	 {
            if(!exists($cluster{$chr}{$i}))
            {
             $cluster{$chr}{$i}=$j;
            }else{$cluster{$chr}{$i}.=",".$j;}
            $state{$chr}[$j]=1;
         }
        }
        }#j 
       }#if 
     }#i
     }#foreach chr
     foreach $chr (keys %cluster)
     {
      foreach $i (keys %{$cluster{$chr}})
      {
       @index=split(/,/,$cluster{$chr}{$i});
       $select=0;
       $selectnum=0;
       $selectscore=0;
       #print OUT "cluster##########\n";
       
       for($j=0;$j<=$#index;$j++)
       { 
        if($back_num{$chr}[$index[$j]]>$selectnum) 
        {
         $selectnum=$back_num{$chr}[$index[$j]];
         $selectscore=$score{$chr}[$index[$j]];
	 $select=$index[$j];
        }elsif($back_num{$chr}[$index[$j]]==$selectnum and $score{$chr}[$index[$j]]>$selectscore)
        {
	  $selectnum=$back_num{$chr}[$index[$j]];
          $selectscore=$score{$chr}[$index[$j]];
	  $select=$index[$j];
	}
       }
       print OUT "$read{$chr}[$select]";
      }
     }

  }#if backnum>1
  else
  {
	  print OUT $read_bk[0];
  }  
######for next###########	 
  %id_cur=();
  $id_cur{$a[0]}=$line;
  @read_bk=();
  $backnum=0;
  $read_bk[$backnum]=$line;
  $backnum+=1;
  %back_num=();
  %back_start=();
  %back_end=();
  %read=();
  %score=();
  %cluster=();
  %state=();
 } 
}
#####when to the end the last has to do one time
 if($backnum>1)
  {
     for($i=0;$i<$backnum;$i++)
     {
       @al=split(/\t/,$read_bk[$i]);
       @pos=split(/,/,$al[2]);
       $read_bk[$i]=~/AS:i:(\S+)/;
       if(!exists($back_num{$pos[0]}))
       {
        $back_num{$pos[0]}[0]=$back_readc{$al[2]};
        $back_start{$pos[0]}[0]=$pos[1];
        $back_end{$pos[0]}[0]=$pos[2];
        $read{$pos[0]}[0]=$read_bk[$i];
        $score{$pos[0]}[0]=$1;
       }else
       {
        $num=scalar(@{$back_start{$pos[0]}});
        $back_num{$pos[0]}[$num]=$back_readc{$al[2]};
        $back_start{$pos[0]}[$num]=$pos[1];
        $back_end{$pos[0]}[$num]=$pos[2];
        $read{$pos[0]}[$num]=$read_bk[$i];
        $score{$pos[0]}[$num]=$1;
       }
      }
     foreach $chr (keys %back_num)
     {
     @{$state{$chr}}=(0) x scalar(@{$back_num{$chr}}); 
     for($i=0;$i<scalar(@{$back_num{$chr}});$i++)
     {
   	  	        	 
       if($state{$chr}[$i]!=1)
       {
        for($j=$i;$j<scalar(@{$back_num{$chr}});$j++)#including i
        {
        $dist=max(abs($back_start{$chr}[$i]-$back_start{$chr}[$j]),abs($back_end{$chr}[$i]-$back_end{$chr}[$j]));
         if($dist<=$range)
         {
	  if(abs($back_start{$chr}[$i]-$back_start{$chr}[$j])==abs($back_end{$chr}[$i]-$back_end{$chr}[$j])){
            if(!exists($cluster{$chr}{$i}))
            {
             $cluster{$chr}{$i}=$j;
            }else{$cluster{$chr}{$i}.=",".$j;}
            $state{$chr}[$j]=1;
         }
         }
        }#j 
       }#if 
     }#i
     }#foreach chr
     foreach $chr (keys %cluster)
     {
      foreach $i (keys %{$cluster{$chr}})
      {
       @index=split(/,/,$cluster{$chr}{$i});
       $select=0;
       $selectnum=0;
       $selectscore=0;
       for($j=0;$j<=$#index;$j++)
       {  
	if($back_num{$chr}[$index[$j]]>$selectnum) 
        {
         $selectnum=$back_num{$chr}[$index[$j]];
         $selectscore=$score{$chr}[$index[$j]];
	 $select=$index[$j];
        }elsif($back_num{$chr}[$index[$j]]==$selectnum and $score{$chr}[$index[$j]]>$selectscore)
        {
	  $selectnum=$back_num{$chr}[$index[$j]];
          $selectscore=$score{$chr}[$index[$j]];
	  $select=$index[$j];
	}
       }
       print OUT "$read{$chr}[$select]";
      }
     }

  }#if backnum>1
  else
  {
    print OUT $read_bk[0];
  }

}
 &merge_read_linear(15);
1;


