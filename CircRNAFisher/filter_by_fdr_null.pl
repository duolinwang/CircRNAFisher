use List::Util qw/sum min/; 
open(IN,"./output/record_nullmodel")or die "!!!";
open(IN2,"./output/re_record_mate_nullscore")or die "!!!";
open(OUT1,">./output/nullmodel_score");
open(OUT2,">./output/realmodel_score");

$count=0;
while($line=<IN>)
{
$count+=1;
if($count%10000==0){print "$key\t$read1score\t$read2score\n";}
#$line=~/[\+|-],(\S+,\d+,\d+),(\d+)\t(\d+)/;
chomp($line);
@a=split(/\t/,$line);
@b=split(/,/,$a[0]);
$key="$b[1],$b[2],$b[3]";
$read1score=$b[4];
$read2score=$a[1];
$score=$read1score+$read2score;
#print "$key\t$read1score\t$read2score\n";

if(!exists($bk{$key}))
{
 $bk{$key}=$score;	
 $bk1{$key}=$read1score;
 $bk2{$key}=$read2score;
 $numbk{$key}=1;
}else{$bk{$key}+=$score;$bk1{$key}+=$read1score;$bk2{$key}+=$read2score;$numbk{$key}+=1;}
}
foreach $k (keys %bk)
{
 $aver=$bk{$k}/$numbk{$k};
 $aver1=$bk1{$k}/$numbk{$k};
 $aver2=$bk2{$k}/$numbk{$k};
 print OUT1 "$k\t$aver\t$aver1\t$aver2\n";
}



############real model####################3
$count=0;
while($line=<IN2>)
{
$count+=1;
if($count%10000==0){print "$key\t$real1score\t$real2score\n";}
chomp($line);
@a=split(/\t/,$line);
@b=split(/:/,$a[2]);
for($i=0;$i<=$#b;$i++)
{
 @c=split(/,/,$b[$i]);
 $key=join(",",@c[1..3]);
 $real1score=$c[4];
 $real2num=scalar(@c)-5;
 if($real2num==1){$real2score=$c[5];}
 else{$real2score=sum(@c[5..$#c])/$real2num;}
 $score=$real1score+$real2score;
 if(!exists($bk_real{$key}))
 {
  $bk_real{$key}=$score;
  $bk_real1{$key}=$real1score;
  $bk_real2{$key}=$real2score;
  $numbk_real{$key}=1;
 }else{$bk_real{$key}+=$score; $bk_real1{$key}+=$real1score;
  $bk_real2{$key}+=$real2score;$numbk_real{$key}+=1;}
}
}

 
foreach $k (keys %bk_real)
{
 $aver=$bk_real{$k}/$numbk_real{$k};
 $aver1=$bk_real1{$k}/$numbk_real{$k};
 $aver2=$bk_real2{$k}/$numbk_real{$k};
 print OUT2 "$k\t$aver\t$aver1\t$aver2\n";
}

