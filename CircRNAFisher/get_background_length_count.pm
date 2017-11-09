#!/usr/bin/env perl

package get_background_length_count;
#use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(get_background_length_count);
our @EXPORT_OK = qw(get_background_length_count);

sub merge #only 3 way at last @merge since sorted mappable!! 1)s in e in 2)s in e out 3)s out e out
{
 my($s,$e,@merged)=@_;
 my @returnmerged=();
 my @left=();
 #3)
 if($s>=($merged[$#merged]+1))
 {
  @returnmerged=(@merged,$s,$e);
  return @returnmerged;
 }
 elsif($e<=$merged[$#merged])# 1)
 {
  return @merged;
 }
 else
 {
   if($s<$merged[$#merged-1]){die "s=$s @merged\n";}
   @left=(@merged[0..($#merged-1)],$e);
   return @left;
 }
}

#compute length and count for each chr
sub print_result
{
 my($chr,$readcount,@chr_length,$readlength)=@_;
 my $mappable_base=0;# whole
 my $temp_mappable=0;# record bin
 my $bincount=scalar(@chr_length);
 for(my $j=0;$j<$bincount;$j+=2)
 {
   my $pos1=$chr_length[$j];
   my $pos2=$chr_length[$j+1];
   my $temp_mappable=$pos2-$pos1;
   if($temp_mappable>=$readlength)
   {
    $mappable_base+=$temp_mappable-$readlength+1;
   }
   else
   {
   $mappable_base+=1; #<100 but backsplice sites so can +1;
   }
 }#for
print OUT "$chr\t$mappable_base\t$readcount\n";
print OUT2 "$chr\t";
for(my $c=0;$c<$#chr_length;$c++)
{
print OUT2 "$chr_length[$c]:";
#print "$chr_length[$c]:";
}
print OUT2 "$chr_length[$#chr_length]\n";
#print "$chr_length[$#chr_length]\n";
#print "$chr\t$mappable_base\t$readcount\n";
}

sub get_background_length_count
{
print "Generating background for the call-peak step...\n";
my $readlength=$_[0];
open(IN1,"./output/backsplice_bed.bed")or die "Can't open file \"./output/backsplice_bed.bed\":$!";
open(IN2,"./output/valid_disorder.bed")or die "Can't open file \"./output/valid_disorder.bed\":$!";
open(IN,"./output/mappable_sorted.bed")or die "Can't open file \"./output//mappable_sorted.bed\":$!";
open(OUT,">./output/background_length_count.txt")or die "$!";
open(OUT2,">./output/mappable_pos.txt")or die "$!";

while($line=<IN1>)
{
@a=split(/\t/,$line);
$readcount{$a[0]}+=0.5;
}

while($line=<IN2>)
{
@a=split(/\t/,$line);
$readcount{$a[0]}+=1;#read1 +1 and next line read2 +1
}


#########################################################
$time=0;
%readid=();
$readcount=0;
while($line=<IN>)
{
 @a=split(/,/,$line);
 $chr=$a[0];
 $pos1=$a[1];
 $pos2=$a[2];
#print "processing $chr,$pos1,$pos2\n";
   if(!exists($chr_count{$chr}))
   {
    if($time!=0)
    {   
        close OUT;
        close OUT2;
        open(OUT,">>./output/background_length_count.txt")or die "$!";
        open(OUT2,">>./output/mappable_pos.txt")or die "$!";
	#print "$readcountin";
        &print_result($chr_pre,$readcount{$chr_pre},@chr_length,$readlength);# print chr_pre
        print "$chr_pre done!\n";
 	%readid=();
    }
    $time=1;
    $chr_count{$chr}=1;#update
    @chr_length=();
    #print "1 pos1=$pos1 pos2=$pos2\n";
    @chr_length=($pos1,$pos2);#update
    #print "1 @chr_length\n";
    $chr_pre=$chr;
    $readcount=1;
   }
   else
   { 
    $readcount+=1;
    if($pos1>$chr_length[$#chr_length])
    {
      push(@chr_length,$pos1);
      push(@chr_length,$pos2);
    }
    elsif($pos2<=$chr_length[$#chr_length])# 1)
    {
       ;
    }
    else
    {
      if($pos1<$chr_length[$#chr_length-1]){die "s=$pos1 @chr_length\n";}
       pop(@chr_length);
       push(@chr_length,$pos2);
    }
   }
}
&print_result($chr_pre,$readcount{$chr_pre},@chr_length,$readlength);# print chr_pre
print "$chr_pre done!\n";

close IN1;
close IN2;
close IN;
close OUT;
close OUT2;
}

1;


