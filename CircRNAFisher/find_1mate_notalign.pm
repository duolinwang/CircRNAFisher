#!   /usr/bin/perl   
package find_1mate_notalign;
#use strict;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(find_1mate_notalign);
our @EXPORT_OK = qw(find_1mate_notalign);


sub find_1mate_notalign
{
my ($reflag)=@_;
open(IN,"./output/".$reflag."mate_not_align")or die "can't open file \"./output/".$reflag."mate_not_align\" : $!";
open(OUT1,">./output/".$reflag."2mate_hassame_pos")or die "$!";
open(OUT2,">./output/".$reflag."2mate_nosame_pos")or die "$!";
open(RECORD,">./output/".$reflag."2mate_realscore_nullscore")or die "$!";
$num_2mate=0;
$line=<IN>;
chomp($line);
$line=~/(\S+\.\d+)_(\d)/;
$rid=$1;
$rAB=$2;
@a=split(/\t/,$line);
$pos{$rid."_".$rAB}=$a[2];
$posnumr=$a[1];
while($line=<IN>)
{
 $line=~/(\S+\.\d+)_(\d)/;
 chomp($line);
 $id=$1;
 $AB=$2;
 @a=split(/\t/,$line);
 $pos{$id."_".$AB}=$a[2];
 $posnum=$a[1];
 if($id eq $rid)
 {
    $num_2mate+=1;
    @b=split(/:/,$pos{$id."_".$AB});
    @br=split(/:/,$pos{$rid."_".$rAB});
    $samepos_num=0;
    @samepos1=();
    #@samepos2=();
    @read1score=();
    @read2score=();
    @null_model=();
    $nulli=0;
    for($bi=0;$bi<=$#b;$bi++)
    {
      @sa=split(/,/,$b[$bi]);       
      $read2score[$bi]=$sa[4];
      for($bri=0;$bri<=$#br;$bri++)
      {
        @sb=split(/,/,$br[$bri]);
        if(($sa[0] ne $sb[0]) and ($sa[1] eq $sb[1]) and ($sa[2] eq $sb[2]) and ($sa[3] eq $sb[3]))###for strand fr!!!!
       { $samepos1[$samepos_num]=$br[$bri].",".$read2score[$bi]; 
	 #$samepos2[$samepos_num]=$b[$bi];
         $samepos_num+=1;
       }#if
        $br[$bri]=~/\S+,\S+,\S+,\S+,(\S+)/;
	$null_model[$nulli]="$1,$read2score[$bi]";
	$nulli+=1;
      }#for bri
    }#forbi
   if($samepos_num==0)##no same pos
   {
     $pposr=$pos{$rid."_".$rAB};
     $ppos=$pos{$id."_".$AB};
     print OUT2 "$rid"."_$rAB\t$posnumr\t$pposr\n";
     #print OUT2 "$id"."_$AB\t$posnum\t$ppos\n";
   }#if
   else
   { $ppos1=join(":",@samepos1);
     #$ppos2=join(":",@samepos2);
     $t=join(",",@read2score);
     print OUT1 "$rid"."_$rAB\t$t\t$ppos1\n";
     print RECORD "$rid"."_$rAB\t$ppos1\n";
     # print OUT1 "$id"."_$AB\t$t\t$ppos2\n";
     $hassame{$rid."_".$rAB}=1;
     $hassame{$id."_".$AB}=1;
     for($nulli=0;$nulli<$#null_model;$nulli++)
     {
     print RECORD "$null_model[$nulli]\t";    
     }
     print RECORD "$null_model[$nulli]\n";
   }#else
  }#if $rid=$id
 else
 {
  $rid=$id;
  $rAB=$AB;
  $posnumr=$posnum;
 }
}
close(IN);
close(OUT1);
close(OUT2);
}


1;
