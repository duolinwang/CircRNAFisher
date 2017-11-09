#!/usr/bin/env perl

sub ifdiscordant{
  my($flaga,$flagb,$chra,$chrb,$posa,$posb,$anchorl)=@_;
  if($chra ne $chrb)
  {
   return 0;
  }
  if($flaga eq "+" and $flagb eq "-") #A + B -
  {
    if($posa-$posb>=$anchorl){return 1;}
    else{return 0;}         
  }
  elsif($flaga eq "-" and $flagb eq "+")#A - B +
  {
    if($posb-$posa>=$anchorl){return 1;}# because use + pos 
    else{return 0;}
  }
  else{return 0;}# one strand+ the other strand-
}

sub find_discordant_readpair
{
print "Detecting discordant readpairs...\n";
my ($unique,$readlength)=@_;
my $infile="";
if($unique ==1 ){$infile="./output/unique_merge_hg19_trans";}
else{$infile="./output/merge_hg19_trans";}
open(IN,$infile)or die "can't open file $infile : $!";
open(OUT,">./output/discordant.bed")or die "$!";
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
my $countline=1;
my ($i,$j);
chomp($line);
@a=split(/\t/,$line);
$a[0]=~/(\S+)_(\d)/;
$id_cur{$1}=1;
if($2 eq "1")
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
 $countline+=1;
 if($countline%10000==0){print "$countline\n";}
 chomp($line);
 @a=split(/\t/,$line);
 $a[0]=~/(\S+)_(\d)/;
 $id=$1;$na=$2;
 #print "$id"."_"."$na\n";
 if(exists($id_cur{$id}))
 {
   if($na eq "1") 
   {$idnum_a+=1;
    $id_A[$idnum_a]=$line;
   }
   elsif($na eq "2")
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
           if($isflag=&ifdiscordant($A[1],$B[1],$A[2],$B[2],$A[3],$B[3],$readlength))
            {
              $can_pairs+=1; 
              print OUT "$A[2]\t$A[3]\t$A[4]\t$A[0]\t255\t$A[1]\n";
              print OUT "$B[2]\t$B[3]\t$B[4]\t$B[0]\t255\t$B[1]\n";
            }
        }#for j
     }#for i   
    }#if 
  #######read next after pre has processed###########
  $idnum_a=-1;$idnum_b=-1;@id_A=();@id_B=();%id_cur=();
  $id_cur{$id}=1;
  if($na eq "1")
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
           if($isflag=&ifdiscordant($A[1],$B[1],$A[2],$B[2],$A[3],$B[3],$readlength))
            { $can_pairs+=1; 
              print OUT "$A[2]\t$A[3]\t$A[4]\t$A[0]\t255\t$A[1]\n";
              print OUT "$B[2]\t$B[3]\t$B[4]\t$B[0]\t255\t$B[1]\n";
            }
        }#for j
     }#for i   
  }#if 
###########the last has to do one time###################################
close IN;
close OUT;
print "Find out $can_pairs discordant read pairs!\n";
}

1;

