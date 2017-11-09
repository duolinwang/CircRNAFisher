#!/usr/bin/env perl

package find_circ;
use strict;

sub ifcirc{
  my($flaga,$flagb,$chra,$chrb,$posa,$posb,$find_range,$anchorl)=@_;
  if($chra ne $chrb)
  {
   return 0;
  }
  if(($posa-$posb)>$find_range or ($posa-$posb)<-$find_range){return 0;}
  if(($flaga&0x0010)==0 and ($flagb&0x0010)==0) #this strand +
  {
    if($posa-$posb>=$anchorl){return 1;}
    else{return 0;}         
  }
  elsif(($flaga&0x0010)!=0 and ($flagb&0x0010)!=0)
  {
    if($posb-$posa>=$anchorl){return 1;}# because use + pos 
    else{return 0;}
  }
  else{return 0;}# one strand+ the other strand-
}

sub ifline{
  my($flaga,$flagb,$chra,$chrb,$posa,$posb,$find_range,$anchorl)=@_;
  if($chra ne $chrb)
  {
   return 0;
  }
  if(($posa-$posb)>$find_range or ($posa-$posb)<-$find_range){return 0;}
  if(($flaga&0x0010)==0 and ($flagb&0x0010)==0) #this strand +
  {
    if($posb-$posa>=$anchorl){return 1;}
    else{return 0;}         
  }
  elsif(($flaga&0x0010)!=0 and ($flagb&0x0010)!=0)
  {
    if($posa-$posb>=$anchorl){return 1;}# because use + pos 
    else{return 0;}
  }
  else{return 0;}# one strand+ the other strand-
}


sub find_circ
{
print "Detecting candidate anchors which in circular mode...\n";
my ($find_range,$anchorl)=@_;
open(IN,"./output/mapped_realign_anchors_k20.sam")or die "can't open file \"./output/mapped_realign_anchors_k20.sam\" : $!";
open(OUT,">./output/candidate_anchors.sam")or die "$!";
open(OUT2,">./output/linear_anchors.sam")or die "$!";
#read the first id############
my $can_pairs=0;
my $idnum_a=-1;
my $idnum_b=-1;
my @a=();
my @id_A=();
my @id_B=();
my %id_cur=();
my $id="";
my $na="";
my @A=();
my @B=();
my $isflag=0;
my $line=<IN>;
my ($i,$j);
chomp($line);
@a=split(/\t/,$line);
$a[0]=~/(\S+_\d)_([AB])/;
$id_cur{$1}=1;
if($2 eq "A")
{
$idnum_a+=1;
$id_A[$idnum_a]=$line;
}
else
{
 $idnum_b+=1;
 $id_B[$idnum_b]=$line;
}
######begin#####################################
while($line=<IN>)
{
 chomp($line);
 @a=split(/\t/,$line);
 $a[0]=~/(\S+_\d)_([AB])/;
 $id=$1;$na=$2;
 #print "$id"."_"."$na\n";
 if(exists($id_cur{$id}))
 {
   if($na eq "A") 
   {$idnum_a+=1;
    $id_A[$idnum_a]=$line;
   }
   elsif($na eq "B")
   {$idnum_b+=1;
    $id_B[$idnum_b]=$line;
   }
 }
 else#### not exists $id_cur{$id}
 { 
   if($idnum_a!=-1 and $idnum_b!=-1)##### not single mapped anchor
   { ####$id_A $id_B##########
     for($i=0;$i<=$idnum_a;$i++)
     {
       @A=split(/\t/,$id_A[$i]);
       for($j=0;$j<=$idnum_b;$j++)
        {
           @B=split(/\t/,$id_B[$j]);
           if($isflag=&ifcirc($A[1],$B[1],$A[2],$B[2],$A[3],$B[3],$find_range,$anchorl))
            {
              $can_pairs+=1; 
              print OUT "$id_A[$i]\n";
              print OUT "$id_B[$j]\n";
            }elsif(&ifline($A[1],$B[1],$A[2],$B[2],$A[3],$B[3],$find_range,$anchorl))
	     {
	      print OUT2 "$id_A[$i]\n";
              print OUT2 "$id_B[$j]\n";
             }
        }#for j
     }#for i   
    }#if 
  #######read next after pre has processed###########
  $idnum_a=-1;$idnum_b=-1;@id_A=();@id_B=();%id_cur=();
  $id_cur{$id}=1;
  if($na eq "A")
  { 
   $idnum_a+=1;
   $id_A[$idnum_a]=$line;
  }
  else
  {
   $idnum_b+=1;
   $id_B[$idnum_b]=$line;
  }#else B
  ###########read next############# 
 }#else #### not exists $id_cur{$id}
}
#####when to the end the last has to do one time
if($idnum_a!=-1 and $idnum_b!=-1)
   { ####$id_A $id_B##########
     for($i=0;$i<=$idnum_a;$i++)
     {
       @A=split(/\t/,$id_A[$i]);
       for($j=0;$j<=$idnum_b;$j++)
        {  @B=split(/\t/,$id_B[$j]);
           if($isflag=&ifcirc($A[1],$B[1],$A[2],$B[2],$A[3],$B[3],$find_range,$anchorl))
            { $can_pairs+=1; 
              print OUT "$id_A[$i]\n";
              print OUT "$id_B[$j]\n";
            }elsif(&ifline($A[1],$B[1],$A[2],$B[2],$A[3],$B[3],$find_range,$anchorl))
	        {
	          print OUT2 "$id_A[$i]\n";
              print OUT2 "$id_B[$j]\n";
            }
        }#for j
     }#for i   
  }#if 
###########the last has to do one time###################################
close IN;
close OUT;
close OUT2;
print "Find out $can_pairs pairs of anchors in circular mode!\n";
}

1;
