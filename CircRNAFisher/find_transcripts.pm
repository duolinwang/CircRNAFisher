#!/usr/bin/env perl

package find_transcripts;

#use strict;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(find_transcripts);
our @EXPORT_OK = qw(find_transcripts);
 

sub sort_annotated_exons
{
my ($known_exons)=@_;	
open(IN,"$known_exons")or die "can't open file $known_exons ";
open(OUT,">./output/known_exons_sort")or die "can't open file ./outpue/known_exons_sort";
$line=<IN>;#head
while($line=<IN>)
{
@a=split(/\t/,$line);
$trans{$a[1]}{$a[3]}{$a[4]}{$a[0]}=$line;
}

foreach $chr (sort keys %trans)#chr        
{
 #print TEST "$chr\n";
 foreach $start (sort{$a<=>$b} keys %{$trans{$chr}})
 {
  #print TEST2 "$chr\t$start\n";
  foreach $end (sort{$a<=>$b} keys %{$trans{$chr}{$start}})
  {
    foreach $key (keys %{$trans{$chr}{$start}{$end}}) 
    {
     print OUT $trans{$chr}{$start}{$end}{$key};
    }
  }
 }
}
close IN;
close OUT;
}

sub sort_my_circRNAs
{
$in="./output/re_spliceplot_withstand.txt";
open(IN,$in)or die "can't";
open(OUT,">./output/sort_my_backsplice");
while($line=<IN>)
{
@a=split(/\t/,$line);
@b=split(/,/,$a[0]);
$pos{$b[0]}{$b[1]}{$b[2]}=$a[0];
}

foreach $chr (sort keys %pos)
{
 foreach $start (sort{$a<=>$b} keys %{$pos{$chr}})
 {
   foreach $end (sort{$a<=>$b} keys %{$pos{$chr}{$start}})
    {
     {print OUT "$pos{$chr}{$start}{$end}\n";}
    }
 }

}
close IN;
close OUT;
}

sub find_transcripts
{
my ($known_exons)=@_;
print "adding annotations to circRNAs....\n";
&sort_annotated_exons($known_exons);
&sort_my_circRNAs();
open(IN1,"./output/known_exons_sort")or die "$!";
open(IN2,"./output/sort_my_backsplice")or die "$!";
open(OUT,">./output/backsplice_transcripts.txt")or die "$!";
open(TEST,">./output/lariat2.txt")or die "$!";
$pos=tell(IN1);
while($line1=<IN1>)
{
 chomp($line1);
 @a=split(/\t/,$line1);
 if(!exists($knownexons{$a[1]}))
 {
 $knownexons{$a[1]}[0]=$line1;
 }else{
  $num=scalar(@{$knownexons{$a[1]}});
  $knownexons{$a[1]}[$num]=$line1;}

}
%print_chr=();

while($line2=<IN2>)
{
chomp($line2);
@b=split(/,/,$line2);
$chr=$b[0];
  if(!exists($print_chr{$chr}))
  {
   print "processing $chr...\n";
   $print_chr{$chr}=1;
   $num=scalar(@{$knownexons{$chr}});
  }
$re_start=$b[1];
$re_end=$b[2];
$can_num=0;
@can_trans_conbk=();
@can_trans_exonlength=();
@can_trans_CDS=();
@can_trans_id=();
@exonstart_list=();
@exonend_list=();
%geneid=();
@lariat=();
for($ki=0;$ki<$num;$ki++)
{
    @uc=split(/\t/,$knownexons{$chr}[$ki]);
    $geneid{$uc[0]}=$uc[10];
    @exonstarts=split(/,/,$uc[8]);
    @exonends=split(/,/,$uc[9]);
    $can_trans_conbk[$can_num]=0;
    $can_tbrans_exonlength[$can_num]=0;
    $can_trans_CDS[$can_num]=0;
    $exonstart_list[$can_num]=-1;
    $exonend_list[$can_num]=-1;
    $lariat[$can_num]=0;
    $ins=-1;$os=-1;$ine=-1;$oe=-1;
    $sdone=0;
    $edone=0;
    $findoverlap=0;
    $start=$re_start;
    $end=$re_end; 
    if($uc[3]<=$re_start)
    {
      if($uc[4]>=$re_end)
      {
       $findoverlap=1;
      }#if($uc[4]>=$end)
      else#uc4<end
      {  if($uc[4]<=$re_start){next;}
         else{$findoverlap=2;   
              $end=$uc[4];
             }
      }#elsif $uc[4]<$end
   }# if($uc[3]<=$re_start)
   else
   {
    if($uc[3]>=$re_end)
    {  
       last;    
    }elsif($uc[4]>=$re_end)
    {$findoverlap=3;
    $start=$uc[3];
   }else#splice site 
    {
    $findoverlap=4;
    $start=$uc[3];
    $end=$uc[4];
    }
   }
    if($findoverlap!=0)
    {
      $can_trans_id[$can_num]=$uc[0]; 
      for($es=0;$es<=$#exonstarts;$es++)
      {
	     if($re_start==$exonends[$es] or $re_end==$exonstarts[$es]){$lariat[$can_num]=1;} 	      
             if($start==$exonstarts[$es]){
                                  if($findoverlap!=3 and $findoverlap!=4){$can_trans_conbk[$can_num]+=1;}
                                  $ins=$es;}
             if($end==$exonends[$es]){
                                  if($findoverlap!=2 and $findoverlap!=4){$can_trans_conbk[$can_num]+=1;}
                                  $ine=$es;}
              ############################################################ 
             if($start<$exonstarts[$es]) 
             {
        	if($sdone==0)
       		{
        	if($es==0){$os=0;}elsif($start>=$exonends[$es-1]){$os=$es;}else{$ins=$es-1;}
        	$sdone=1;     
       		}
     	     }elsif($es==$#exonstarts and $sdone==0)#
      	     {
          	$ins=$es;$sdone=1; 
             }
            if($end<$exonends[$es])#
      	    {
      		if($edone==0)
      		{
		if($end<=$exonstarts[$es]){$oe=$es-1;}else{$ine=$es;}
                $edone=1; 
     		}
     	    }
      }#for
      if($ins!=-1)
      {
       if($ine!=-1)
       {
          if($ins==$ine)
          {
           $can_trans_exonlength[$can_num]=$end-$start;
           $exonstart_list[$can_num]=$re_start;
           $exonend_list[$can_num]=$re_end;
          }
          else 
          {
           $can_trans_exonlength[$can_num]+=$exonends[$ins]-$start;
           $exonstart_list[$can_num]=$re_start;
           $exonend_list[$can_num]=$exonends[$ins]; 
           for($i=$ins+1;$i<=$ine-1;$i++)
           {
            $can_trans_exonlength[$can_num]+=$exonends[$i]-$exonstarts[$i];
            $exonstart_list[$can_num].=",".$exonstarts[$i];
            $exonend_list[$can_num].=",".$exonends[$i];
           }
           $can_trans_exonlength[$can_num]+=$end-$exonstarts[$ine];
           $exonstart_list[$can_num].=",".$exonstarts[$ine];
           $exonend_list[$can_num].=",".$re_end;
          }
       }#if($ine!=-1)
       elsif($ine==-1)
       {
           $can_trans_exonlength[$can_num]+=$exonends[$ins]-$start;
           $exonstart_list[$can_num]=$re_start;
           if($oe==$ins){$exonend_list[$can_num]=$re_end;}else{$exonend_list[$can_num]=$exonends[$ins];}
           for($i=$ins+1;$i<=$oe;$i++)
           {
            $can_trans_exonlength[$can_num]+=$exonends[$i]-$exonstarts[$i];
            $exonstart_list[$can_num].=",".$exonstarts[$i];
            if($i<$oe){$exonend_list[$can_num].=",".$exonends[$i];}else{$exonend_list[$can_num].=",".$re_end;} 
           }
       } 
     }#if($ins!=-1)
     else#start out
     {
      if($ine!=-1)#3
      {
       if($os==$ine)
       {
        $can_trans_exonlength[$can_num]=$end-$exonstarts[$os];
        $exonstart_list[$can_num]=$re_start;
        $exonend_list[$can_num]=$re_end;
       }
       else
       {
         for($i=$os;$i<=$ine-1;$i++)
         {
         $can_trans_exonlength[$can_num]+=$exonends[$i]-$exonstarts[$i];
         if($i==$os){$exonstart_list[$can_num]=$re_start;
                     $exonend_list[$can_num]=$exonends[$i];
                    }
         else{$exonstart_list[$can_num].=",".$exonstarts[$i];
              $exonend_list[$can_num].=",".$exonends[$i];
             }
         }
         $can_trans_exonlength[$can_num]+=$end-$exonstarts[$ine];
         $exonstart_list[$can_num].=",".$exonstarts[$ine];
         $exonend_list[$can_num].=",".$re_end;
       }#esle
      }else#end out3 start
      {
         if($os>$oe){$can_trans_exonlength[$can_num]+=0;
                     $exonstart_list[$can_num]=$re_start;
                     $exonend_list[$can_num]=$re_end;
                    }
         if($os==$oe){$can_trans_exonlength[$can_num]+=$exonends[$os]-$exonstarts[$os];
                     $exonstart_list[$can_num]=$re_start;
                     $exonend_list[$can_num]=$re_end;
                    }
         else{
               for($i=$os;$i<=$oe;$i++)
               {
                 $can_trans_exonlength[$can_num]+=$exonends[$i]-$exonstarts[$i];
                 if($i==$os){$exonstart_list[$can_num]=$re_start;
                             $exonend_list[$can_num]=$exonends[$i]; 
                    }
                 elsif($i==$oe){$exonstart_list[$can_num].=",".$exonstarts[$i];
                                $exonend_list[$can_num].=",".$re_end; 
                    } 
                 else{$exonstart_list[$can_num].=",".$exonstarts[$i];
                      $exonend_list[$can_num].=",".$exonends[$i];
                     }
            }#for
         }#else
      }
     }#else#s out
     if($start<=$uc[5]){$mins=$uc[5];}else{$mins=$start;}
     if($end<=$uc[6]){$mine=$end;}else{$mine=$uc[6];}
     if($end<=$uc[5] or $start>$uc[6]){ $can_trans_CDS[$can_num]=0;}
     else{$can_trans_CDS[$can_num]=$mine-$mins;}
     $can_num+=1;
    }#if(findoverlap) 
}#for ki
     if($can_num==0)
     {
       print OUT "$line2\tNone\n";
     }
     else
     {
           if($can_num==1)
           {
            if($can_trans_exonlength[0]==0 and $lariat[0]==1)
	    {
	    print TEST2 "$line2\t$can_trans_id[0]\t$can_trans_conbk[0]\t$can_trans_exonlength[0]\t$can_trans_CDS[0]\t$exonstart_list[0]\t$exonend_list[0]\t$geneid{$can_trans_id[0]}\n";	  
	    }else{
	    print OUT "$line2\t$can_trans_id[0]\t$can_trans_conbk[0]\t$can_trans_exonlength[0]\t$can_trans_CDS[0]\t$exonstart_list[0]\t$exonend_list[0]\t$geneid{$can_trans_id[0]}\n";
           }
           }
           else
           {
             $first=0;
             %forsort=();
             $has=0;
             for($i=0;$i<=$can_num-1;$i++)
             {
             $forsort{$can_trans_conbk[$i]}{$can_trans_exonlength[$i]}{$can_trans_CDS[$i]}=$i;
             if($can_trans_exonlength[$i] == 0 and $lariat[$i]==1 and $has==0)
             {$has=1;
              print TEST "$line2\t$can_trans_id[$i]\t$can_trans_conbk[$i]\t$can_trans_exonlength[$i]\t$can_trans_CDS[$i]\t$exonstart_list[$i]\t$exonend_list[$i]\t$geneid{$can_trans_id[$i]}\n";
	     }
             }
             foreach $conbk (sort{$b<=>$a} keys %forsort)
             {
                foreach $exonl (sort{$b<=>$a} keys %{$forsort{$conbk}})
                {  
                  foreach $cds (sort{$b<=>$a} keys %{$forsort{$conbk}{$exonl}})
                  {
                    if($first==0)
                    {
		       if($has==0)
	               {	       
                       $index=$forsort{$conbk}{$exonl}{$cds};
                       print OUT "$line2\t$can_trans_id[$index]\t$can_trans_conbk[$index]\t$can_trans_exonlength[$index]\t$can_trans_CDS[$index]\t$exonstart_list[$index]\t$exonend_list[$index]\t$geneid{$can_trans_id[$index]}\n";
	               }
		       $first=1;
                    }#if
                  }#foreach
                }#foreach
             }#foreach
            }#if($can_num==1)
       }#if($can_num!=0)
}#while line2
close IN1;
close IN2;
close OUT;
close TEST;
}

1;
