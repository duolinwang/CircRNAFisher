#!/usr/bin/env perl

package extension_line;
use strict;
use List::Util qw/max/;
use List::Util qw/min/;

sub transtorev{
my ($utrbase)=@_;
my @base=split(//,$utrbase);
for(my $n=0;$n<=$#base;$n++){
if($base[$n] eq 'A'){$base[$n]='T';}
elsif($base[$n] eq 'T'){$base[$n]='A';}
elsif($base[$n] eq 'C'){$base[$n]='G';}
elsif($base[$n] eq 'G'){$base[$n]='C';}
elsif($base[$n] eq 'a'){$base[$n]='t';}
elsif($base[$n] eq 't'){$base[$n]='a';}
elsif($base[$n] eq 'c'){$base[$n]='g';}
elsif($base[$n] eq 'g'){$base[$n]='c';}
}
my $utrbase=join('',@base);
return scalar reverse($utrbase);
}
###########################################
our %hg19=();
sub readhg19{
my ($chrid,$refseq)=@_;
my $filename="$refseq".$chrid.".fa";
print "open reference $filename\n";
if(!open(IN2,$filename)){print "Don't have reference file $filename\n";return 0;};
my $l=<IN2>;###read first line#######
$hg19{$chrid}='';###global####3
while($l=<IN2>)
{
chomp($l);
$hg19{$chrid}.=$l;
}
close IN2;
return 1;
}

sub extension_line
{
our $testmultiple=0;
our $testunique=0;
print "candidate anchors extension...\n";
my ($refseq,$misextend,$mismatch_thres,$anchorl)=@_; 
open(IN1,"./output/linear_anchors.sam")or die "Can't open file \"./output/linear_anchors.sam\" : $!";
open(OUT,">./output/mapped_sample_linear") or die "$!";
#open(IN1,"./output/test_linear_anchors.sam")or die "$!";
#open(TEST1,">./output/can_extensionread") or die "$!";
#5'-3' intron motif GT-AG,GC-AG,AT-AC
my @splice_motif_3=("AG","AC","AG","GC","AC","AT");#3' intron motif B 
my @splice_motif_5=("GT","CT","GC","CT","AT","GT");#5' intron motif A
my @motif_score=(2*$misextend+2,2*$misextend+2,($misextend+1),($misextend+1),$misextend,$misextend);
my $mismatchscore=-1;
my $thres=6;
my (@a,@mismatchposA,@mismatchposB,@ss,@mismatchnum_record,@mismatchpos_record,@maxscore,@splice_site,@flankA,@flankB,@mispos_anchorA,@mispos_anchorB);
my ($anchor_misnum,$anchor_misnumA,$anchor_misnumB,$iA,$iB,$readid,$pos,$orread,$cur,$rcur,$blength_flag,$mismatchA,$mismatchB,$mismatchnum,$mismatchpos,$pos,$posA,$posB,$AB,$flag_motif,$poshg19_A,$poshg19_B,$A_cur,$B_cur,$splice_n,$bk,$bk_rev,$orread_rev,$mA_n,$mB_n,$Acur,$Bcur,$line,$readbk_s,$readbk_e);
while($line=<IN1>)
{
#print "$line";	
 @a=split(/\t/,$line);
if($a[2]!~/^(chr\d+|chrX|chrY|chrM)$/){next;}######filter out chr###############
 $pos=$a[3];###1-based need to change to offset 0
if($a[0]=~/(\S+_\d)_A_(\S+)/) ###caution later if add read to B!!!
{
$readid=$1;
$AB='A';####from A or B;
$orread=$2;
$line =~ /XM:i:(\S+)/;
$anchor_misnumA=$1;
@mispos_anchorA=();
if($anchor_misnumA!=0)
{
	$line=~/MD:Z:(\S+)/;
	my @miscode=split(/\D/,$1);
	$mispos_anchorA[0]=$miscode[0];#0 based from anchor A
	for(my $i=1;$i<=$#miscode-1;$i++)
	{
 	$mispos_anchorA[$i]=$mispos_anchorA[$i-1]+$miscode[$i]+1;
	}
}
}
else{
$AB='B';
$line =~ /XM:i:(\S+)/;
$anchor_misnumB=$1;
@mispos_anchorB=();
if($anchor_misnumB!=0)
{
	$line=~/MD:Z:(\S+)/;
	my @miscode=split(/\D/,$1);
	$mispos_anchorB[0]=$miscode[0];#0 based from anchor B
	for(my $i=1;$i<=$#miscode-1;$i++)
	{
 	$mispos_anchorB[$i]=$mispos_anchorB[$i-1]+$miscode[$i]+1;
	}
}
}
#print "$readid\n";
if(!exists($hg19{$a[2]}))
{
  if(&readhg19($a[2],$refseq)==0){next;}
}

if(($a[1]&0x0010)==0)####strand +
{
   my $strand="+";
   if($AB eq 'A')
   {   $mismatchA=0;
       $iA=-1;  #last match pos(start=0)
       @mismatchposA=();########first time
       my $blength_flag=0;
       $posA=$pos;
       my ($cur,$or);
      #$test=substr($hg19{$a[2]},$posA-1+20,90);
      #print "$test\n";
     while($mismatchA<=$misextend)###
     {
        $iA+=1;
        if($iA==(length($orread)-$anchorl)){$blength_flag=1;last;}
        $cur=substr($hg19{$a[2]},$posA-1+$anchorl+$iA,1);
        $or=substr($orread,$iA+$anchorl,1);
	#print "$cur----$or\n";
        if($or!~/$cur/i)
        {
         $mismatchposA[$mismatchA]=$iA+$anchorl;###read pos not hg19 from left to right! from 0
	 #print "A:$mismatchA:$mismatchposA[$mismatchA]\n"; 
         $mismatchA+=1;
        }
     }#while mismatchA
     if($blength_flag!=1)#break mimatch>=3
     {
        $iA-=1;
        $mismatchA-=1; 
        #while(($mismatchA-1)>=0 and $mismatchposA[$mismatchA-1]==$iA+20 and $iA!=-1){$iA-=1;$mismatchA-=1;}###for end-end last iA is match pos
        $mA_n=$iA+1;#######match number##########
     }
     else{$mA_n=$iA+1;}
     #print "mA=$mA_n\n";
   }#if $AB= 'A'
   else###($AB eq 'B')
   {
       $mismatchB=0;
       @mismatchposB=();
       $iB=-1;#last match pos(start=0)
       my $blength_flag=0;
       $posB=$pos;
       if(($posB-$posA)<(length($orread)-$anchorl)){next;}#gap
       my ($cur,$or);
      #$test=substr($hg19{$a[2]},($posB-2-90),91);
      #print "$test\n";
     while($mismatchB<=$misextend)
     {
        $iB+=1;
        if($iB==(length($orread)-$anchorl)){$blength_flag=1;last;}#
        $cur=substr($hg19{$a[2]},$posB-2-$iB,1);
        $or=substr($orread,-$iB-$anchorl-1,1);   #-$iB-21 because from orread substr!! not pos!!
	#print "$cur----$or\n";
        if($or!~/$cur/i)
        {
         $mismatchposB[$mismatchB]=length($orread)-$anchorl-1-$iB;#read pos from left to right!! not iB
	 #print "B:$mismatchB:$mismatchposB[$mismatchB]\n"; 
         $mismatchB+=1;
        }
     }#while mismatchB
     if($blength_flag!=1)#break mimatch>3
     {   
        $iB-=1;
        $mismatchB-=1; 
	while(($mismatchB-1)>=0 and $mismatchposB[$mismatchB-1]==length($orread)-$anchorl-1-$iB and $iB !=-1){$iB-=1;$mismatchB-=1;}##for end-end last iB is match pos
        $mB_n=$iB+1;##match number for anchor B
     }
     else{next;}
       #  print "mB=$mB_n\n";
   ######B over need to extend##########
   if(($mA_n+$mB_n)>=(length($orread)-2*$anchorl) and $mA_n<=(length($orread)-2*$anchorl-1+$thres) and $mB_n<=(length($orread)-2*$anchorl-1+$thres)) ###has possible splice site#######$mA_n mB——n不能过于长
   {
     $readbk_s=length($orread)-$iB-$anchorl-1;###orread pos (from 0)not hg19 from B end include B
     $readbk_e=$iA+$anchorl;  ###read pos not hg19!!! from 0 A end include A 
     ##############detect splice site############
      @splice_site=();
      @flankA=();
      @flankB=();
      $splice_n=0;
      @maxscore=();
      @mismatchnum_record=();
      @mismatchpos_record=();
      $anchor_misnum=$anchor_misnumA+$anchor_misnumB;
        	my $tempA=();
		for(my $mpi<=0;$mpi<=$#mispos_anchorA;$mpi++)
		{
		 $tempA.=$mispos_anchorA[$mpi].";";
		}
		my $tempB=();
	        for(my $mpi<=0;$mpi<=$#mispos_anchorB;$mpi++)
		{
		 $mispos_anchorB[$mpi]+=length($orread)-$anchorl;
		 $tempB.=$mispos_anchorB[$mpi].";";
		}
     #######!!!! $rcur is read cur from 0 offset 0 not hg19 can be trans to hg19 cur
     for($rcur=$readbk_s;$rcur<=$readbk_e+1;$rcur++)#rcur=$readbk_e+1 bk=$readbk_e no B extension
     {
         $bk=$rcur-1;
	 #print "rcur=$rcur\n";
	 #if(!($bk ~~ @mismatchposA or (($bk+1) ~~ @mismatchposB)))
         {
         ########################################
          $Acur=$rcur+$posA-1; ###A hg19 cur
          $Bcur=$posB-2-(length($orread)-$anchorl-1-$rcur); ##B hg19 cur 
         $flankA[$splice_n]=substr($hg19{$a[2]},$Acur,2);###strand +16
         $flankB[$splice_n]=substr($hg19{$a[2]},$Bcur-2,2); ####strand +
         #print "flanka:$flankA-------flankb:$flankB\n"; 
         #print "flanka:$flankA[$splice_n]-------flankb:$flankB[$splice_n]\n"; 
               my $mismatchnumA=$mismatchA;
               my $mismatchnumB=$mismatchB;
               $mismatchpos=();
              #$simple_flag=0;########if extension has aaaacccccggggtttt
	      #print "before:A:$mismatchnumA:B:$mismatchnumB\n";
              for(my $mpi=0;$mpi<=($mismatchA-1);$mpi++)
              {
                  if($mismatchposA[$mpi]>=$rcur){$mismatchnumA-=1;}
		  else{$mismatchpos.=$mismatchposA[$mpi].";";}
              }
              for(my $mpi=0;$mpi<=($mismatchB-1);$mpi++)
              { 
                if($mismatchposB[$mpi]<$rcur){$mismatchnumB-=1;}
		else{$mismatchpos.=$mismatchposB[$mpi].";";}
              }
               $mismatchnum=$mismatchnumA+$mismatchnumB+$anchor_misnum;
	       # print "mismatchnum=$mismatchnum\t$mismatchnumA\t$mismatchnumB\n";
         if($mismatchnum<=$mismatch_thres)#mismatch filter
         {
              $flag_motif=-1;
         	for(my $spmotif=0;$spmotif<=$#splice_motif_3;$spmotif++)
         	{  
           	 if((($flankA[$splice_n] =~/$splice_motif_5[$spmotif]/i) and ($flankB[$splice_n] =~/$splice_motif_3[$spmotif]/i)))
           	{
               	$flag_motif=$spmotif;
               	#print "find $flankA[$splice_n] $flankB[$splice_n]\n";
               	last; #find jump out
           	}#if flank
         	}#for motif
              	if($flag_motif !=-1)
              	{
              		$maxscore[$splice_n]=$mismatchnum*$mismatchscore+$motif_score[$flag_motif];
              	}else{$maxscore[$splice_n]=$mismatchnum*$mismatchscore;}
		#if($maxscore[$splice_n]>=0) #if maxscore>=0 pass for non motif mismatch==0
	      	{
              	$splice_site[$splice_n]=$rcur;# $bk is on the left of splice_site 
		$mismatchnum_record[$splice_n]=$mismatchnum; 
		$mismatchpos_record[$splice_n]=$tempA.$mismatchpos.$tempB;
               	$splice_n+=1;
              	}
          }#mismatchnum past filter
          }#if not mismatch at bk left and right not lariat 
          }#for $rcur  
              #########test##########################################################
              #return maxscore backsplice if multiple return all#record multiple same score to see number######
              if($splice_n>=1)
              {
              	my $max=max(@maxscore);
	      	my $printline="";
	      	my $minstart=0;
		my $maxi=0;
                #print OUT "range=$range==================================================\n";
              	for(my $si=0;$si<$splice_n;$si++)
              	{
                	if($maxscore[$si]==$max)
                	{
			 my $poshg19_A=$splice_site[$si]+$posA-1;# on the left of bk
			 my $poshg19_B=$posB-2-(length($orread)-$anchorl-1-$splice_site[$si]);# on the right of bk
			 my $bk=$splice_site[$si]-1;
			 my $hg19start=$poshg19_A; ###same as paper read start 0 based 
			 my $hg19end=$poshg19_B;   ###same as paper read end not include so not -1
		         if($maxi==0){$minstart=$hg19start;}
			 if($hg19start<$minstart or $maxi==0)
			 {
		 	  $printline="$readid\t$strand\t$a[2]\t$hg19start\t$hg19end\t*\t$mismatchnum_record[$si]\t$mismatchpos_record[$si]\n";#bk is read pos bk on the bk left!
		         }
			 $maxi+=1;
		        }#if == maxreturn all splice which same as max
              	}#for
		#if($maxi==1){$testunique+=1;}else{$testmultiple+=1;}
	      	print OUT $printline;
		#print $printline;
		#print "$testunique $testmultiple\n";

              }
    }#if(($mA_n+$mB_n)>=(length($orread)-40)) ###has possible splice site#######
   }#elsif B  AB 
}#if strand +
###############################################################################################################
###############################################################################################################
###############################################################################################################
else####strand - $orread_rev=&transtorev($orread); change A-B every B compute as A
{
   $orread_rev=&transtorev($orread);
   my $strand="-";
   if($AB eq 'A')###($AB eq 'A')
   {
       $mismatchA=0;
       @mismatchposA=();########first time
       $iA=-1;#last match pos(start=0)
       my $blength_flag=0;
       $posA=$pos;
       my ($cur,$or);
      #$test=substr($hg19{$a[2]},($posA-2-90),91);
      #print "hg19A:$test\n";
     while($mismatchA<=$misextend)
     {
        $iA+=1;
        if($iA==(length($orread_rev)-$anchorl)){$blength_flag=1;last;}
        $cur=substr($hg19{$a[2]},$posA-2-$iA,1);
        $or=substr($orread_rev,-$iA-$anchorl-1,1);
	#print "$cur----$or\n";
        if($or!~/$cur/i)
        {
         $mismatchposA[$mismatchA]=length($orread_rev)-$iA-$anchorl-1;###read pos not hg19 from B to A not(A to B left to right)
	 #print "mismatchposA:$mismatchA:$mismatchposA[$mismatchA]\n";
         $mismatchA+=1;
        }      
     }#while mismatchA
     if($blength_flag!=1)#break mimatch>3
     {
        $iA-=1;
        $mismatchA-=1;
        #while(($mismatchA-1)>=0 and $mismatchposA[$mismatchA-1]==length($orread_rev)-$iA-$anchorl-1 and $iA!=-1){$iA-=1;$mismatchA-=1;}###for end-end last iA is last match pos
        $mA_n=$iA+1;##match number for anchor B
	#print "entendA=$mA_n\n";  
     }
     else{$mA_n=$iA+1;}
          #print "mA=$mA_n\n";
    }#if A
   else ############if $AB=B
   {   $mismatchB=0;
       @mismatchposB=();########first time
       $iB=-1;  #last match pos(start=0)
       $blength_flag=0;
       $posB=$pos;
       if(($posA-$posB)<(length($orread)-$anchorl)){next;}#gap
       my ($cur,$or);
      #$test=substr($hg19{$a[2]},$posB-1+20,90);
      #print "hg19B:$test\n";
     while($mismatchB<=$misextend)
     {
        $iB+=1; 
        if($iB==(length($orread_rev)-$anchorl)){$blength_flag=1;last;}
        $cur=substr($hg19{$a[2]},$posB-1+$anchorl+$iB,1);
        $or=substr($orread_rev,$iB+$anchorl,1);
	#print "$cur----$or\n";
        if($or!~/$cur/i)
        {
         $mismatchposB[$mismatchB]=$iB+$anchorl;#rev_read pos if A of orread is on left then this pos is from right to left
	 #print "mismatchposB:$mismatchB:$mismatchposB[$mismatchB]\n";
         $mismatchB+=1;
        }     
     }#while mismatch
     if($blength_flag!=1)#break mimatch>3
     {
        $iB-=1;### last match pos(from 0)
        $mismatchB-=1;
	while(($mismatchB-1)>=0 and $mismatchposB[$mismatchB-1]==($iB+$anchorl) and $iB!=-1){$iB-=1;$mismatchB-=1;}###for end-end last iB is match pos
        $mB_n=$iB+1;#######match number##########
     }
     else{next;}
   ######B over need to extend##########
   #print "extension:$mA_n-$mB_n\n";
   if(($mA_n+$mB_n)>=(length($orread_rev)-$anchorl*2) and $mA_n<=(length($orread_rev)-$anchorl*2-1+$thres) and $mB_n<=(length($orread_rev)-2*$anchorl-1+$thres)) ###has possible splice site#######
   {
      $readbk_s=length($orread_rev)-$iA-$anchorl-1;###orread_rev pos (from 0)
      $readbk_e=$iB+$anchorl;
     #$test=substr($orread_rev,$readbk_s,$readbk_e-$readbk_s+1);
     ##############detect splice site############
      @splice_site=();
      @flankA=();
      @flankB=();
      $splice_n=0;
      @maxscore=();
      @mismatchnum_record=();
      @mismatchpos_record=();
      $anchor_misnum=$anchor_misnumA+$anchor_misnumB; 
                        my $tempB=();
			for(my $mpi<=0;$mpi<=$#mispos_anchorB;$mpi++)
			{
		 	$tempB.=$mispos_anchorB[$mpi].";";
			}
			my $tempA=();
	        	for(my $mpi<=0;$mpi<=$#mispos_anchorA;$mpi++)
			{
		 	$mispos_anchorA[$mpi]+=length($orread)-$anchorl;
		 	$tempA.=$mispos_anchorA[$mpi].";";
			}
		
     for($rcur=$readbk_s;$rcur<=$readbk_e+1;$rcur++)#rcur=$readbk_e+1 bk=$readbk_e no B extension
     {
       $bk_rev=$rcur-1;
       #print "rcur=$rcur\n";
       #if(!($bk_rev ~~ @mismatchposB or (($bk_rev+1)~~ @mismatchposA)))#bk read pos on the backsplice left,0-based	     
       { #print "bk no mismatch\n";
        $Bcur=$rcur+$posB-1;#B hg19 cur
        $Acur=$posA-2-(length($orread_rev)-$anchorl-1-$rcur);#A hg19 cur
       $flankB[$splice_n]=substr($hg19{$a[2]},$Bcur,2);###strand +
       $flankA[$splice_n]=substr($hg19{$a[2]},$Acur-2,2); ####strand +
       #print "not rev:$flankA$splice_n---$flankB$splice_n\n";
       $flankB[$splice_n]=&transtorev($flankB[$splice_n]);###turn to - revcomplement seq
       $flankA[$splice_n]=&transtorev($flankA[$splice_n]);
       #print "flankA=$flankA[$splice_n],flankB=$flankB[$splice_n]\n";
        my $mismatchnumA=$mismatchA;
        my $mismatchnumB=$mismatchB;
        $mismatchpos=();
	       for(my $mpi=0;$mpi<=($mismatchB-1);$mpi++)###
              {  
                 if($mismatchposB[$mpi]>=$rcur){$mismatchnumB-=1;}
		 else{$mismatchpos.=$mismatchposB[$mpi].";";}#mismatchpos from B to A
              }
              for(my $mpi=0;$mpi<=($mismatchA-1);$mpi++)
              {
                  if($mismatchposA[$mpi]<$rcur){$mismatchnumA-=1;}#rcur and mismatchpos is read pos
		  else{$mismatchpos.=$mismatchposA[$mpi].";";}
              }
             
               $mismatchnum=$mismatchnumA+$mismatchnumB+$anchor_misnum;
	       #print "mismatchpos=$mismatchpos\n";
	       # print "mismatchnum=$mismatchnum\t$mismatchnumA\t$mismatchnumB\n";
              if($mismatchnum<=$mismatch_thres)#mismatch filter
             {
              	        $flag_motif=-1;
              		for(my $spmotif=0;$spmotif<=$#splice_motif_3;$spmotif++)
             		{  
               		 if((($flankA[$splice_n] =~/$splice_motif_5[$spmotif]/i) and ($flankB[$splice_n] =~/$splice_motif_3[$spmotif]/i)))
              		 {
               			$flag_motif=$spmotif;
               			#print "find motif and mismatchnum=$mismatchnum\n";
               			last; #find jump out
              		 }#if flank
            		}#for motif

              		if($flag_motif !=-1)
              		{
              		 $maxscore[$splice_n]=$mismatchnum*$mismatchscore+$motif_score[$flag_motif];
             		}else{$maxscore[$splice_n]=$mismatchnum*$mismatchscore;}
			#	if($maxscore[$splice_n]>=0) #if maxscore>=0 pass for non motif mismatch==0
			{
              		$splice_site[$splice_n]=$rcur;# $bk is on the left of splice_site 
		        $mismatchnum_record[$splice_n]=$mismatchnum; 
			$mismatchpos_record[$splice_n]=$tempB.$mismatchpos.$tempA;	
			$splice_n+=1;
                	}
      		} #filter by mismatchnum
               }#if not mismatch at bk left and right not lariat 
              }#for $rcur
             if($splice_n>=1)
             {
                my $max=max(@maxscore);
         	my $printline="";
	      	my $minstart=0;
		my $maxi=0;
                for(my $si=0;$si<$splice_n;$si++)
                {
			if($maxscore[$si]==$max)
              		{ 
				my $poshg19_B=$splice_site[$si]+$posB-1;# on the right of bk
				my $poshg19_A=$posA-2-(length($orread)-$anchorl-1-$splice_site[$si]);# on the right of bk
				#my $bk_rev=$splice_site[$si]-1;
				#my $bk=length($orread_rev)-1-$bk_rev-1;##-1 again bs bk on the left of orread !!! from A to B 5'---3'
                 		my $hg19start=$poshg19_B; ###same as paper read start 0 based 
				my $hg19end=$poshg19_A;   ###same as paper read end not include so not -1
				if($maxi==0){$minstart=$hg19start;}#select left most
				if($hg19start<$minstart or $maxi==0)
				 {   
				 $printline="$readid\t$strand\t$a[2]\t$hg19start\t$hg19end\t*\t$mismatchnum_record[$si]\t$mismatchpos_record[$si]\n";#bk is read pos bk on the bk left!
			         }
				$maxi+=1;
              		}#if max
              	}#for
		#if($maxi==1){$testunique+=1;}else{$testmultiple+=1;}
		print OUT $printline;
		#print $printline;
		#print "$testunique $testmultiple\n";
		
              }#if >=1;
    }#if(($mA_n+$mB_n)>=(length($orread_rev)-40)) ###has possible splice site#######
   }#elsif B  AB 
}#if strand -
}###while read line
%hg19=();
close IN1;
close OUT;
}

#&extension_line("/home/wangdu/data/circRNA/hg19/",3,3,20);
1;
