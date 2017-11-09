#!/usr/bin/env perl

package get_splice_site_plot_strand;

#use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(get_splice_site_plot_strand);
our @EXPORT_OK = qw(get_splice_site_plot_strand);

sub get_splice_site_plot_strand
{
my ($unique,$readlength)=@_;
my $in1;
if($unique){$in1="./output/re_pass_paired_end_read_unique.txt";}else{$in1="./output/re_pass_paired_end_read.txt";}
print "IN1=$in1\n";
open(IN1,$in1)or die "Can't open file $in1 : $!";
open(IN2,"./output/re_splice_sites_read_mismatch")or die "Can't open file \"./output/re_splice_sites_read_mismatch\" : $!";
open(OUT,">./output/re_spliceplot_withstand.txt")or die "$!";

while($line1=<IN1>)
{
chomp($line1);
@a=split(/\t/,$line1);
@b=split(/:/,$a[2]);
for($i=0;$i<=$#b;$i++)
{
 @c=split(/,/,$b[$i]);
 $key=join(",",@c[0..3]); 
 $passsplit{$key}=$c[5];
 $passsread_mq{$a[0]}{$key}=$c[4];
 $passread{$a[0]}[$i]=$key;
}
}
while($line2=<IN2>)
{
  if($line2 !~ /range/)
  {
  chomp($line2);
  @l1=split(/\t/,$line2);
  @ll1=split(/,/,$l1[1]);
  if($ll1[0]!~/^(chr\d+|chrX|chrY|chrM)$/){<IN2>;<IN2>;<IN2>;next;}###过滤掉chrUn
  @l2=split(/:/,$l1[3]);###read flank!!########
  $strand=$l2[0];
  $score=$l2[5];
  $mismatchnum=$l2[6];
  $readid=$l1[0];
  $key=$l1[1];
  $key=$strand.",".$key;
  $Aextend=$l1[2]+1;
  $Bextend=$readlength-$Aextend; 
  if($strand eq '+'){$left=$Bextend;$right=$Aextend;}else{$left=$Aextend;$right=$Bextend;}

  if(exists($passsplit{$key}) and exists($passread{$readid}) and ($key ~~@{$passread{$readid}}))####加这个就是所有通过了paired end的read和splice才算
  {
    @strdlist=split(/,/,$key);
    $strand=$strdlist[0];
    $mq=$passsplit{$key};
    $eachmq=$passsread_mq{$readid}{$key};
    $key=~s/[\+|-],//; #delete strand;
    if(!exists($myread{$key}))
    { 
      $myreadnum{$key}=1;
      $myread{$key}{0}=$readid.":".$strand.":".$score.":".$mismatchnum.":".$eachmq.":".$left.",".$right; #this way for sort later 
      $myreadleft{$key}{0}=$left;
      $backmq{$key}=$mq;
      $backscore{$key}=$score;
      $backmismatchnum{$key}=$mismatchnum;
    }
    else
    {
      $myreadnum{$key}+=1;      
      $myread{$key}{$myreadnum{$key}}=$readid.":".$strand.":".$score.":".$mismatchnum.":".$eachmq.":".$left.",".$right;
      $myreadleft{$key}{$myreadnum{$key}}=$left;
      $backmq{$key}=$mq;
      $backscore{$key}+=$score;
      $backmismatchnum{$key}+=$mismatchnum;
    }
  }
}#if range
}#while IN2

foreach $t (sort{$myreadnum{$b}<=>$myreadnum{$a}}keys (%myread))
{
 print OUT "$t\t";
 print OUT "$myreadnum{$t}\t";
 $fn=0;
 foreach $tt (sort{$myreadleft{$t}{$a}<=>$myreadleft{$t}{$b}} keys %{$myread{$t}})
 {
  $fn+=1;
  if($fn<$myreadnum{$t})
  {
   print OUT "$myread{$t}{$tt};";
  }
  else
  {
  print OUT "$myread{$t}{$tt}\t";
  }
 }
 $score=$backscore{$t}/$myreadnum{$t};
 $s=sprintf "%0.2f",$score; 
 print OUT "$s\t";
 $mismatchnum=$backmismatchnum{$t}/$myreadnum{$t};
 $s=sprintf "%0.2f",$mismatchnum;
 print OUT "$s\t";
 $s=sprintf "%0.2f",$backmq{$t};
 print OUT "$s\n";
}
close IN1;
close IN2;
close OUT;
}
#&get_splice_site_plot_strand("re_",1,76);

1;
