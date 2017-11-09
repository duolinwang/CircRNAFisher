#!/usr/bin/env perl

package change_transpos2hg19pos;

#use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(change2hg19);
our @EXPORT_OK = qw(change2hg19);
sub get_Q
{
 	my ($MD,@phred)=@_;
	my $Q=0;
	my @miscode=split(/\D/,$MD);
	my @misca=split(/\d+/,$MD);
	my @mispos=();
   	$mispos[0]=$miscode[0]+1;
   	for(my $i=1;$i<=$#miscode-1;$i++)
   	{
    	$mispos[$i]=$mispos[$i-1]+$miscode[$i]+1;
   	}
        for($i=0;$i<=$#mispos;$i++)
        {
 		if($misca[$i+1] eq "N")#misca from 1
 		{
 		$Q+=2;
 		}else{
 		$Q+=ord($phred[$mispos[$i]-1])-33;
 		}
	}	
	return $Q;
}

sub change2hg19
{
my ($known_exons)=@_;
my $infile1="./output/mapped_sample_transcript.sam";
open(IN1,$infile1)or die "Can't open file $infile1 : $!";
open(IN2,"$known_exons")or die "Can't open file $known_exons : $!";
open(OUT,">./output/mapped_sample_trans_hg19")or die "$!";
<IN2>;
while($line2=<IN2>)
{
chomp($line2);
@uc=split(/\t/,$line2);
$uc[0]=~s/\.\d+//;
$ucsc_chr{$uc[0]}=$uc[1];
$ucsc_strand{$uc[0]}=$uc[2];
$ucsc_exonstarts{$uc[0]}=$uc[8];
$ucsc_exonends{$uc[0]}=$uc[9];
}

while($line1=<IN1>)
{
if(($line1=~/XG:i:0/))
{
	@a=();
	chomp($line1);
	@a=split(/\t/,$line1);
	$transid=$a[2];
	$transid=~s/\.\d+//;
	$leftpos=$a[3];
	$strand=$a[1];
	$readlength=length($a[9]);
	$transexon_l=0;
	#$leftsam=join("\t",@a[4..$#a]);## not exactly right!!
	if(!exists($ucsc_chr{$transid})){die "no $transid in ucsc_refseq_exonpos\n";}
	@exonstarts=split(/,/,$ucsc_exonstarts{$transid});
	@exonends=split(/,/,$ucsc_exonends{$transid});
	if($ucsc_strand{$transid} eq '-')
	{
   	 for($i=0;$i<=$#exonends;$i++)
  	 {
         $transexon_l+=$exonends[$i]-$exonstarts[$i];
         }
         $leftpos=$transexon_l-$leftpos+1; 
         if(($strand&0x0010)==0)
         {
          #$strand+=16;
          $strand="-";
         }
         else
        {
         #$strand-=16;
         $strand="+"
        }
         }
       else
       {
          if(($strand&0x0010)==0)
          {
           $strand="+";
          }
          else
          {
           $strand="-"
          }
}

    $exonsum=0;
    for($i=0;$i<=$#exonstarts;$i++)
    {
      $exonsum+=$exonends[$i]-$exonstarts[$i];
      if($leftpos<=$exonsum)
      {
        $left_leftpos=$leftpos-$exonsum+$exonends[$i]-$exonstarts[$i];
        if($ucsc_strand{$transid} eq '+'){$newpos1=$exonstarts[$i]+$left_leftpos-1;}
        if($ucsc_strand{$transid} eq '-'){$newpos2=$exonstarts[$i]+$left_leftpos;}
        last;
      }
    }

   if($ucsc_strand{$transid} eq '+')
   {
    $leftposend=$leftpos-1+$readlength;
   }
   else
   {
    $leftposend=$leftpos-$readlength+1; 
   }

   $exonsum=0;
   for($i=0;$i<=$#exonstarts;$i++)
   {
    $exonsum+=$exonends[$i]-$exonstarts[$i];
    if($leftposend<=$exonsum)
    {
     $left_leftpos=$leftposend-$exonsum+$exonends[$i]-$exonstarts[$i];
     if($ucsc_strand{$transid} eq '+'){$newpos2=$exonstarts[$i]+$left_leftpos;}
     if($ucsc_strand{$transid} eq '-'){$newpos1=$exonstarts[$i]+$left_leftpos-1;}
     last;
    }
}

  $line1=~/XM:i:(\d+)/;
  $mis_num=$1;
  if($mis_num!=0)
  {
   $line1=~/AS:i:(\S+)/;
   $score=$1;
   @phred=split(//,$a[10]);
   $line1=~/MD:Z:(\S+)/;
   $Q=&get_Q($1,@phred);
  }else{$Q=0;$score=0;} 
 print OUT "$a[0]\t$strand\t$ucsc_chr{$transid}\t$newpos1\t$newpos2\t$score\t$mis_num\t$Q\n";
}#if gap
}#while
close IN1;
close IN2;
close OUT;
}#change2hg19 

1;
