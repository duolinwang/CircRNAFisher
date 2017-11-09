#!/usr/bin/env perl

package pre_paired_filter;
#use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(changsam2fastq get_non_lariat_read merge_hg19_trans);
our @EXPORT_OK = qw(changsam2fastq get_lariat_read merge_hg19_trans);

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

sub changsam2fastq
{
print "changing unmapped reads into fastq format (changsam2fastq)...\n";	
open(IN0,"./output/re_mapped_linear.sam")or die "!!!";
open(IN,"./output/unmapped_sample.sam")or die "can't open file \"./output/unmapped_sample.sam\" : $!";
open(OUT,">./output/unmapped_back.fastq") or die "$!";
my $line;
my %lineread=();
while($line=<IN0>)
{
chomp($line);
my @a=split(/\t/,$line);
$lineread{$a[0]}=1;
}

while($line=<IN>)
{
chomp($line);
my @a=split(/\t/,$line);
if(!exists($lineread{$a[0]}))
{
print OUT "\@$a[0]\n";
print OUT "$a[9]\n";
print OUT "+\n";
print OUT "$a[10]\n";
}
}
close IN;
close OUT;
close IN0;
}

sub get_lariat_read
{
print "Filter out lariat reads...\n";
my ($readlength,$anchorl)=@_;
open(IN,"./output/re_mapped.sam")or die "Can't open file \"./output/re_mapped.sam\" : $!"; #no gaps
open(OUT,">./output/lariat_read_backsplice")or die "$!";
open(OUT2,">./output/re_splice_sites_read_mismatch")or die "$!";
my $seqlength=$readlength-$anchorl;

my $line;
my $strand;
my $Aext=0;
my $mismatchnum=0;
my $backid="";
my $readid="";
my $leftmost=0;
my $readlength=0;#update
while($line=<IN>)
{
 my @a=split(/\t/,$line);
 $readlength=length($a[9]);
 if(($a[1]&0x0010) == 0){$strand="+";}
 else{$strand="-";}
 $backid=$a[2];
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
     	  if($mispos[$i]==($seqlength-$leftmost+1) or $mispos[$i]==($seqlength-$leftmost+2))
	  {$flag=1;#print "$line\n";!!!!wrong!!!!! cigar is read pos not refseq pos!!!!
	  }
   	 if($misca[$i+1] eq "N")#misca from 1
         {
           $Q+=2;#
         }else{
	   $Q+=ord($samphred[$mispos[$i]-1])-33;	
              } 
        }
   	if($flag==1)
   	{
          print OUT "$readid\t$backid\n";	
	}
  }#has mismatch 
 if($strand eq "+"){$Aext=$seqlength-$leftmost;}else{$Aext=$readlength-($seqlength-$leftmost+1)-1;}
           print OUT2 "$readid\t$backid\t$Aext\t";
	   print OUT2 "$strand:A::B::$score:$mismatchnum:$Q\t$readlength\n";
}#while

close IN;
close OUT;
close OUT2;
}

sub merge_hg19_trans
{
print "Merging mapped sequences...\n";
my ($readlength)=@_;	
open(IN1,"./output/mapped_sample_k20.sam") or die "Can't open file \"./output/mapped_sample_k20.sam\" : $!";#hg19
open(IN2,"./output/mapped_sample_trans_hg19") or die "Can't open file \"./output/mapped_sample_trans_hg19\" : $!"; #trans changed2 hg19 pos
open(OUT,">./output/merge_hg19_trans") or die "$!";
open(OUT2,">./output/unique_merge_hg19_trans") or die "$!";
$count=0;
if($line1=<IN1>)
{
 @a=split(/\t/,$line1);
 $a[0]=~/\S+\.(\d+)_(\d)/;
 $id1=$1.$2;#change2num;
 $readid1{$id1}=1;
 #print "hg19 read: $id1\n";
 if(($a[1]&0x0010)==0){$strand1="+";}else{$strand1="-";}
 $chr1=$a[2];
 $start1=$a[3]-1;
 $end1=$start1+$readlength;
 $line1=~/XM:i:(\S+)/;
 $mismatchnum1=$1;
 if($mismatchnum1!=0)#has mismatch
 {
 @phred=split(//,$a[10]);
 $line1=~/AS:i:(\S+)/;
 $score1=$1;
 $line1=~/MD:Z:(\S+)/;
 $Q=&get_Q($1,@phred);
 }else{$Q=0;$score1=0;}
 $merge{"$id1\t$chr1\t$strand1\t$start1"}="$a[0]\t$strand1\t$chr1\t$start1\t$end1\t$score1\t$mismatchnum1\t$Q\n";
 $idlist{"$id1\t$chr1\t$strand1\t$start1"}=$id1; 
 $score{"$id1\t$chr1\t$strand1\t$start1"}=$score1;
 
 $inpos=tell(IN2);
 while($line2=<IN2>)
 {
  @b=split(/\t/,$line2);
  $b[0]=~/\S+\.(\d+)_(\d)/;
  $id2=$1.$2;
    if($id2>$id1)
    {       seek(IN2,$inpos,0); 
	    last;}
    else{
	  $inpos=tell(IN2);#update inpos
      }
      #print "hg19:$id1-trans:$id2\n";
      #print "trand read: $id2\n";
  $strand2=$b[1];
  $chr2=$b[2];
  $start2=$b[3];
  $end2=$b[4];
  $score2=$b[5];
  if(!exists($merge{"$id2\t$chr2\t$strand2\t$start2"}))
  {	  
    $merge{"$id2\t$chr2\t$strand2\t$start2"}=$line2;
    $idlist{"$id2\t$chr2\t$strand2\t$start2"}=$id2;
    $score{"$id2\t$chr2\t$strand2\t$start2"}=$score2;
  }elsif($score2>$score{"$id2\t$chr2\t$strand2\t$start2"}){
    $merge{"$id2\t$chr2\t$strand2\t$start2"}=$line2;
    $idlist{"$id2\t$chr2\t$strand2\t$start2"}=$id2;
    $score{"$id2\t$chr2\t$strand2\t$start2"}=$score2;
  }
 }#while line2
}
#################################
while($line1=<IN1>)
{
 $count+=1;
 if($count%10000==0){print "$count\n";} 
 @a=split(/\t/,$line1);
 $a[0]=~/\S+\.(\d+)_(\d)/;
 $id1=$1.$2;#chang2num;
 #print "hg19 read: $id1\n";
  if(!exists($readid1{$id1}))
  {
	  #  print "hg19 new read :$id1\n";
    #print last ones#######################
    %numid=();
    foreach $m (sort{$idlist{$a}<=>$idlist{$b}}keys %idlist) #sort by id
    {
     print OUT $merge{$m};
     if(!exists($numid{$idlist{$m}}))
     {
      $numid{$idlist{$m}}=1;
     }else{$numid{$idlist{$m}}+=1;}
    }   
    foreach $m (sort{$idlist{$a}<=>$idlist{$b}}keys %idlist)
    {
      if($numid{$idlist{$m}}==1)
      {
      print OUT2 $merge{$m};
      }
    }
    %readid1=();
    $readid1{$id1}=1;
    %readid2=();
    $readid2{$id2}=1;
    %merge=();####wrong!!!!!!
    %score=();
    %idlist=(); 
  }#new one print over

   if(($a[1]&0x0010)==0){$strand1="+";}else{$strand1="-";}
   $chr1=$a[2];
   $start1=$a[3]-1;
   $end1=$start1+$readlength;
   $line1=~/XM:i:(\S+)/;
   $mismatchnum1=$1;
   if($mismatchnum1!=0)
   {
    @phred=split(//,$a[10]);    
    $line1=~/AS:i:(\S+)/;
    $score1=$1;
    $line1=~/MD:Z:(\S+)/;
    $Q=&get_Q($1,@phred);
   }else{$Q=0;$score1=0;}
   if(!exists($merge{"$id1\t$chr1\t$strand1\t$start1"}))
   {#print "no $chr1\t$strand1\t$start1 add one with score $score1\n";	 
   $merge{"$id1\t$chr1\t$strand1\t$start1"}="$a[0]\t$strand1\t$chr1\t$start1\t$end1\t$score1\t$mismatchnum1\t$Q\n";
   $score{"$id1\t$chr1\t$strand1\t$start1"}=$score1;
    $idlist{"$id1\t$chr1\t$strand1\t$start1"}=$id1;
   }elsif($score1>$score{"$id1\t$chr1\t$strand1\t$start1"}){
	  #print "has $chr1\t$strand\t$start1\t but score is big update score as $score1\n";
    $merge{"$id1\t$chr1\t$strand1\t$start1"}="$a[0]\t$strand1\t$chr1\t$start1\t$end1\t$score1\t$mismatchnum1\t$Q\n";
    $idlist{"$id1\t$chr1\t$strand1\t$start1"}=$id1;
    $score{"$id1\t$chr1\t$strand1\t$start1"}=$score1;}
   
 #$inpos=tell(IN2); only need first time
 while($line2=<IN2>)
 {
  @b=split(/\t/,$line2);
  $b[0]=~/\S+\.(\d+)_(\d)/;
  $id2=$1.$2;
    if($id2>$id1)
    {
      seek(IN2,$inpos,0);last;
    }
    else{
	  $inpos=tell(IN2);#update inpos
      }
   #print "hg19:$id1-trans:$id2\n";
  #print "trand read: $id2 $line2\n";
  $strand2=$b[1];
  $chr2=$b[2];
  $start2=$b[3];
  $end2=$b[4];
  $score2=$b[5];
  if(!exists($merge{"$id2\t$chr2\t$strand2\t$start2"}))
  {
    #print "no trans $chr2\t$strand2\t$start2 add one with score $score2\n";  
    $merge{"$id2\t$chr2\t$strand2\t$start2"}=$line2;
    $idlist{"$id2\t$chr2\t$strand2\t$start2"}=$id2;    
    $score{"$id2\t$chr2\t$strand2\t$start2"}=$score2;
  }elsif($score2>$score{"$id2\t$chr2\t$strand2\t$start2"}){
	  #print "has trans $chr2\t$strand2\t$start2\t but score is big update score as $score2\n";
    $merge{"$id2\t$chr2\t$strand2\t$start2"}=$line2;
    $score{"$id2\t$chr2\t$strand2\t$start2"}=$score2;
    $idlist{"$id2\t$chr2\t$strand2\t$start2"}=$id2;
  }
  }#while line2
}#while line1

print "last and left\n";
#print last ones and line2 left ##################################################
    %numid=();
    foreach $m (sort{$idlist{$a}<=>$idlist{$b}}keys %idlist) #sort by id
    {
     print OUT $merge{$m};
     if(!exists($numid{$idlist{$m}}))
     {
      $numid{$idlist{$m}}=1;
     }else{$numid{$idlist{$m}}+=1;}
    }   
    
    foreach $m (sort{$idlist{$a}<=>$idlist{$b}}keys %idlist)
    {
      if($numid{$idlist{$m}}==1)
      {
      print OUT2 $merge{$m};
      }
   }
   %merge=();
   %score=();
   $idlist=();

while($line2=<IN2>)
{
 if($count%10000==0){print "$count\n";}   	
  @b=split(/\t/,$line2);
  $b[0]=~/\S+\.(\d+)_(\d)/;
  $id2=$1.$2;
  $strand2=$b[1];
  $chr2=$b[2];
  $start2=$b[3];
  $end2=$b[4];
  $score2=$b[5];
  if(!exists($merge{"$id2\t$chr2\t$strand2\t$start2"}))
  {
	  #print "no $chr2 $strand2 $start2\n";
    $merge{"$id2\t$chr2\t$strand2\t$start2"}=$line2;
    $score{"$id2\t$chr2\t$strand2\t$start2"}=$score2;
    $idlist{"$id2\t$chr2\t$strand2\t$start2"}=$id2;
  }elsif($score2>$score{"$id2\t$chr2\t$strand2\t$start2"}){
	  #print "has trans $chr2\t$strand2\t$start2\t but score is big update score as $score2\n";
    $merge{"$id2\t$chr2\t$strand2\t$start2"}=$line2;
    $score{"$id2\t$chr2\t$strand2\t$start2"}=$score2;
    $idlist{"$id2\t$chr2\t$strand2\t$start2"}=$id2;    
    }
}
  
    %numid=();
    foreach $m (sort{$idlist{$a}<=>$idlist{$b}}keys %idlist) #sort by id
    {
     print OUT $merge{$m};
     if(!exists($numid{$idlist{$m}}))
     {
      $numid{$idlist{$m}}=1;
     }else{$numid{$idlist{$m}}+=1;}
    }   
    
    foreach $m (sort{$idlist{$a}<=>$idlist{$b}}keys %idlist)
    {
      if($numid{$idlist{$m}}==1)
      {
      print OUT2 $merge{$m};
      }
   }

close IN1;
close IN2;
close OUT;
close OUT2;
print "Merging mapped sequences over!\n";

}

1;
