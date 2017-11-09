#!/usr/bin/env perl

package call_peak;

use threads;
use Math::CDF;
use Statistics::R;
#use strict;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(peak_count compute_pvalue  padjust_and_filter);
our @EXPORT_OK = qw(peak_count compute_pvalue padjust_and_filter);
 
sub peak_count
{	
     print "Generating peaks...\n"; 
     my ($flariat)=@_;	
     if($flariat==1){$in1="./output/re_spliceplot_withstrand_detetelariat.txt"}else{$in1="./output/re_spliceplot_withstand.txt";}
     open(IN1,$in1)or die "Can't open file \"$in1\" : $!";
     print "open $in1\n";
     #open(IN1,"./output/merged_backsplice.txt")or die "Can't open file \"./output/merged_backsplice.txt\": $!";
     open(IN2,"./output/pairedreads_match_splicesite_within_fragment")or die "Can't open file \"./output/pairedreads_match_splicesite_within_fragment\": $!";
     open(OUT,">./output/peak_count") or die "$!";
     my $line;
     while($line=<IN1>)
     {
     chomp($line);
     @a=split(/\t/,$line);
     $back_count{$a[0]}=$a[1];
     }
     
     while($line=<IN2>)
     {
      chomp($line);
      @a=split(/\t/,$line);
      @b=split(/:/,$a[2]);
      $sum=0;
      for($i=0;$i<=$#b;$i++)
      {	 
       $sum+=$back_count{$b[$i]};
      }
      for($i=0;$i<=$#b;$i++)
      {
        if(!exists($dis_count{$b[$i]}))
        {
        $dis_count{$b[$i]}=$back_count{$b[$i]}/$sum;
        $dis_count_or{$b[$i]}=1;
        }
        else
        {    
        $dis_count{$b[$i]}+=$back_count{$b[$i]}/$sum;
        $dis_count_or{$b[$i]}+=1;
        }
      }
     }
     
     foreach $k (sort{$back_count{$b}<=>$back_count{$a}}keys %back_count)
     {
      if(exists($dis_count{$k}))
      {
        $count=$back_count{$k}+$dis_count{$k};
      }
      else{$count=$back_count{$k};$dis_count{$k}=0;$dis_count_or{$k}=0;}
      print OUT "$k\t$count\t$back_count{$k}\t$dis_count{$k}\t$dis_count_or{$k}\n";
     }
     
     close IN1;
     close IN2;
     close OUT;
}



sub compute_pvalue
{
        my ($fragmentsize,$readlength)=@_;
        open(IN1,"./output/peak_count")or die "Can't open file \"./output/peak_count\": $!";	
        open(IN2,"./output/background_length_count.txt")or die "Can't open file \"./output//background_length_count.txt\": $!";
	open(OUT,">./output/peak_pvalue.txt")or die "$!";

        my ($line2,@a,%bigsize,%wholenum,$line1,$chr,$peaknum,$n,$N,$p,$mean,$var,$key);
	my $window=($fragmentsize-$readlength)*2;
	while($line2=<IN2>)
	{
		chomp($line2);
		@a=split(/\t/,$line2);
		$bigsize{$a[0]}=$a[1];
		$wholenum{$a[0]}=$a[2];
	}
	my $linenum=0;
	while($line1=<IN1>)
	{
		chomp($line1);
		$linenum+=1;
		@a=split(/\t/,$line1);
		@b=split(/,/,$a[0]);
		$chr=$b[0];
		$peaknum=$a[1];#junction read and discordant read
		$backnum=$a[2];
        ###z-score
	$n=$peaknum;
	$N=$wholenum{$chr};
	$p=($window-$readlength+1)/$bigsize{$chr};
	$var=sqrt($N*$p*(1-$p));
	$mean=$N*$p;
	#$z=($n-$mean)/$var;
	$key=$a[0];
        $pvalue=1-&Math::CDF::pbinom($n, $N, $p);
        print OUT "$key\t$n\t$a[2]\t$a[3]\t$a[4]\t$mean\t$pvalue\n";
        }
close IN1;
close IN2;
close OUT;
}



sub padjust_and_filter
{
	my ($cutoff,$method)=@_;
        open(IN,"./output/peak_pvalue.txt")or die "Can't open file ./output/peak_pvalue : $! please run this pipeline step by step!";
        open(IN0,"./output/backsplice_transcripts.txt")or die "Can't open file ./output/backsplice_transcripts : $! please run this pipeline step by step!";
        $in1="./output/re_spliceplot_withstand.txt";
        open(IN1,$in1)or die "Can't open file \"$in1\" : $!";
        print "open $in1\n";
	open(OUT,">./output/backsplice_passcutoff.txt")or die "$!";
	while($line=<IN0>)
	{
	 chomp($line);
	 @a=split(/\t/,$line);
	 $back_trans{$a[0]}=join("\t",@a[1..$#a]);
	}
	close IN0;
	while($line=<IN1>)
	{
	 chomp($line);
	 @a=split(/\t/,$line);
         $a[3]=sprintf "%0.2f",$a[3];
         $a[4]=sprintf "%0.2f",$a[4]; 
	 $back_mq{$a[0]}=join(",",@a[3..5]);
	}
	close IN1;
	while($line=<IN>)
	{
	  chomp($line);
	  my @a=split(/\t/,$line);
	  $pvalue=$a[6];
	  @b=split(/,/,$a[0]);
	  $chr=$b[0];
	  if(!exists($pvaluelist{$chr}))
	  {
	   $pvaluelist{$chr}[0]=$pvalue;
	   $backsplice{$chr}[0]=$line;
	   $backkey{$chr}[0]=$a[0];
	  }
	  else
	  {
 	   $index=scalar(@{$pvaluelist{$chr}});
 	   $pvaluelist{$chr}[$index]=$pvalue;
 	   $backsplice{$chr}[$index]=$line;
	   $backkey{$chr}[$index]=$a[0];
          }
	}
        close IN;
        print "Adjusting p-values for multiple comparisons by method \"$method\"\n";
	my %filter_line=();	
	foreach $chr (keys %pvaluelist)
	{
	my $R = Statistics::R->new();
	$R->run('rm(p); rm(qvalue); p<-1');
        $index=scalar(@{$pvaluelist{$chr}});
        for($i=0;$i<$index;$i++)
        {      $in=$i+1;
	       $p='p['.$in.']';
               $R->set($p,$pvaluelist{$chr}[$i]);
        }
	$R->set('m',"$method");
	$R->run('qvalue<-p.adjust(p, method = m)');
	$qvalue=$R->get('qvalue');
        $R->stop();
	delete($pvaluelist{$chr});
     for($i=0;$i<$index;$i++)
     {
      if($$qvalue[$i] eq ''){$$qvalue[$i]=0;}
            
      if($$qvalue[$i]<=$cutoff and exists($back_trans{$backkey{$chr}[$i]}))
      {
  	$filter_line{$backsplice{$chr}[$i]}="$$qvalue[$i]\t$back_mq{$backkey{$chr}[$i]}\t$back_trans{$backkey{$chr}[$i]}";
      }
     }
     delete($back_trans{$chr});
     delete($backkey{$chr});
     }#foreach chr
   	foreach $key (sort{$filter_line{$a}<=>$filter_line{$b}} keys %filter_line)
	{
	  print OUT "$key\t$filter_line{$key}\n";
	}
        close OUT;
}

1;



