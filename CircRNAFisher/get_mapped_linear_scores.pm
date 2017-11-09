#!/usr/bin/env perl

package get_mapped_linear_scores;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(get_mapped_linear_scores);
our @EXPORT_OK = qw(get_mapped_linear_scores);

use List::Util qw/min/;
sub get_mapped_linear_scores
{
my ($readlength,$anchorl)=@_;
open(IN,"./output/re_mapped_linear.sam")or die "Can't open file \"./output/re_mapped.sam\" : $!"; #no gaps
open(OUT,">./output/mapped_sample_linear_score")or die "$!";
my $seqlength=$readlength-$anchorl;
my $line;
my $strand;
my $chr;
my $pos1;
my $pos2;
my $start,$end;
my $mismatchnum=0;
my $readid="";
my $leftmost=0;
my $readlength=0;#update
while($line=<IN>)
{
 my @a=split(/\t/,$line);
 @b=split(/,/,$a[2]);
 $pos1=$b[1];
 $pos2=$b[2];
 $chr=$b[0];
 $readlength=length($a[9]);
 if(($a[1]&0x0010) == 0){$strand="+";}
 else{$strand="-";}
 $readid=$a[0];
 $line=~/XM:i:(\S+)/;
 $mismatchnum=$1;
 $line=~/AS:i:(\S+)/;
 $score=$1;
 $leftmost=$a[3];
 $Q=0;
 if($mismatchnum !=0 )#must have mismatch
 {
   	my $flag=0;
        @samphred=split(//,$a[10]);
	$line=~/MD:Z:(\S+)/;
        @miscode=split(/\D/,$1);
        @misca=split(/\d+/,$1);
        @mispos=();
   	$mispos[0]=$miscode[0]+1;
   	for(my $i=1;$i<=$#miscode-1;$i++)
   	{
    	$mispos[$i]=$mispos[$i-1]+$miscode[$i]+1;
   	}
	for(my $i=0;$i<=$#mispos;$i++)
   	{
   	 if($misca[$i+1] eq "N")#misca from 1
         {
           $Q+=2;#
         }else{
	   $Q+=ord($samphred[$mispos[$i]-1])-33;	
              } 
        }
  }#has mismatch 
           $start=$pos1-($seqlength-$leftmost+1)+1;
           $end=$pos2+$readlength-($seqlength-$leftmost+1);         
           print OUT "$readid\t$strand\t$chr\t$start\t$end\t$score\t$mismatchnum\t$Q\n"; 
}#while
close IN;
close OUT;
}

1;


