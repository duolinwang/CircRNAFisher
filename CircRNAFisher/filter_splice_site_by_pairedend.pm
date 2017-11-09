#!   /usr/bin/perl   
package filter_splice_site_by_pairedend;
#use strict;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(filter_splice_site_by_pairedend);
our @EXPORT_OK = qw(filter_splice_site_by_pairedend);
use List::Util qw/min max/; 

our $Inmatepos=0;
our $Inmagepos2=0;

sub read_mate
{
 my($find_id)=@_;
 my $da=0;
 my $db=0;
 my $matestrand;
 my $matetype;
 my $matepos;
 my $matechr;
 my $mate_alignscore;
 my $findk=20; #ps -k 15 so 20 is enough
 my $findflag=0;
 my $mateid=0;
 my %mate=(); ##to record mateid
 my %matehg19=();
 my @a;
 my $line;
 if($find_id =~/(\S+\.\d+_)1/){$find_id=$1."2";}
 elsif($find_id =~/(\S+\.\d+_)2/){$find_id=$1."1";}
 ($da=$find_id)=~s/_//;$da =~/\S+\.(\d+)/;$da=$1; 
  ##############################################################################
  $findk=20;
  while($line=<Inmate>)
  { chomp($line);
    @a=split(/\t/,$line);
    $mateid=$a[0];
    ($db=$mateid)=~s/_//;$db=~/\S+\.(\d+)/;$db=$1;
    if($da<$db){last;} ###can't find findid
    if($da>$db){$Inmatepos=tell(Inmate);}
    if(($mateid eq $find_id) and $findk>=1)
    { 
      $findflag=1;#both find
      $findk-=1;
      $Inmatepos=tell(Inmate);
      if(!exists($mate{$mateid}))
      {$mate{$mateid}=1;
       $matepos1=$a[3];
       $matechr=$a[2];
       $mate_alignscore=$a[7];
       $matestrand=$a[1];
       $matetype="1";
       $matepos2=$a[4];# right open
       ####strand + 
       $matehg19{"$mateid,$matestrand,$matechr,$matepos1,$matepos2"}=1;
      }
      else#######exists matepos mateid,mate multiple align
      {
       if(!exists($matehg19{"$mateid,$a[1],$a[2],$a[3],$a[4]"}))
       {
        $matepos1.="\t".$a[3];
        $matepos2.="\t".$a[4];
        $matechr.="\t".$a[2];
        $mate_alignscore.="\t".$a[7];
	$matestrand.="\t".$a[1];
        $matetype.="\t1";
        $matehg19{"$mateid,$a[1],$a[2],$a[3],$a[4]"}=1;
       }
      }###else exists matepos mateid
    }#find id 
  }#while line 
######################################################### 
 $findk=20;
 while($line=<Inmate2>)
 { chomp($line);
  @a=split(/\t/,$line);
  $mateid=$a[0];
  ($db=$mateid)=~s/_//;$db=~/\S+\.(\d+)/;$db=$1;
  if($da<$db){last;} ###can't find findid
  if($da>$db){$Inmatepos2=tell(Inmate2);}
  if(($mateid eq $find_id) and $findk>=1)
  { 
      $findflag=1;#both find
      $findk-=1;
      $Inmatepos2=tell(Inmate2);
      if(!exists($mate{$mateid}))
      {$mate{$mateid}=1;
       $matepos1=$a[3];
       $matechr=$a[2];
       $mate_alignscore=$a[7];
       $matestrand=$a[1];
       $matetype="1";
       $matepos2=$a[4];# right open
       ####strand + 
       $matehg19{"$mateid,$matestrand,$matechr,$matepos1,$matepos2"}=1;
      }
      else#######exists matepos mateid,mate multiple align
      {
       if(!exists($matehg19{"$mateid,$a[1],$a[2],$a[3],$a[4]"}))
       {
        $matepos1.="\t".$a[3];
        $matepos2.="\t".$a[4];
        $matechr.="\t".$a[2];
        $mate_alignscore.="\t".$a[7];
	$matestrand.="\t".$a[1];
        $matetype.="\t1";
        $matehg19{"$mateid,$a[1],$a[2],$a[3],$a[4]"}=1;
       }
      }###else exists matepos mateid
   }#find id 
  }#while line 
  #print "want to find mate $find_id\n";
  #print "mate_bk=$mate_bk{$find_id}\n"; 
 if(exists($mate_bk{$find_id}))
 {
     $findflag=1;
     #print "find !! $find_id\n";
     @mateline=split(/\t/,$mate_bk{$find_id});
     @mateposarray=split(/:/,$mateline[2]);
     for($matein=0;$matein<=$#mateposarray;$matein++)
     {
      $mateid=$find_id;
      @mateinfor=split(/,/,$mateposarray[$matein]);
      if(!exists($mate{$mateid}))
      {$mate{$mateid}=1;
       $matepos1=$mateinfor[2];
       $matechr=$mateinfor[1];
       $mate_alignscore=$mateinfor[4];
       $matestrand=$mateinfor[0];
       $matetype="2";
       $matepos2=$mateinfor[3];# right open
       ####strand + 
      }
      else#######exists matepos mateid,mate multiple align
      {
       
        $matepos1.="\t".$mateinfor[2];
        $matepos2.="\t".$mateinfor[3];
        $matechr.="\t".$$mateinfor[1];
        $mate_alignscore.="\t".$mateinfor[4];
	$matestrand.="\t".$mateinfor[0]; 
        $matetype.="\t2";	
      }###else exist
     }
 }

 seek(Inmate,$Inmatepos,0); ###update mate pos 
 seek(Inmate2,$Inmatepos2,0); ##update mate pos
 if($findflag==1)#find at least one
 { 
 return ($matestrand,$matechr,$matepos1,$matepos2,$mate_alignscore,$matetype);
 }
 else
 { 
 return (-1,-1,-1,-1,-1);#### mate not align
 }
}


sub filter_splice_site_by_pairedend
{
my ($reflag)=@_;
open(IN1,"./output/".$reflag."my_readcount_id_mismatch") or die "Can't open file \"./output/".$reflag."my_readcount_id_mismatch\" : $!";
open(Inmate,"./output/merge_hg19_trans") or die "Can't open file \"./output/merge_hg19_trans\" : $!";
open(Inmate2,"./output/mapped_sample_linear_score") or die "Can't open file \"./output/mapped_sample_linear_score\" : $!";
open(IN2,"./output/my_readcount_linear")or die "Can't open file ./output/my_readcount_linear :$!";
open(OUT1,">./output/".$reflag."mate_not_align") or die "$!";
open(OUT2,">./output/".$reflag."mate_align_right") or die "$!";
open(OUT3,">./output/".$reflag."mate_align_wrong") or die "$!";
open(RECORD3,">./output/record_realscore_nullmodel")or die "$!";
open(TEST,">./output/testlinearnull")or die "$!";
open(TEST2,">./output/test2mate")or die "$!";
open(TEST3,">./output/test_mate_infor")or die "$!";

while($line0=<IN2>)
{
  chomp($line0);
  my @b=split(/\t/,$line0);
  $linear_read{$b[0]}=$line0;
}

%mate_bk=();
my %splice_read_ini=();
while($line1=<IN1>)
{
 	
  chomp($line1);
  my @b=split(/\t/,$line1);
  $splice_read_ini{$b[0]}=$line1;
  #print "$b[0]=$line1\n";
}

%splice_read=();
foreach $readid (keys %splice_read_ini)
{
  $readid=~/(\S+)_(\d)/;
  #print "readid=$readid-$1-$2\n";
  if($2 eq "1"){$mateid="$1_2";}else{$mateid="$1_1";} 
  #print "$readid-$mateid\n";
  if(exists($splice_read_ini{$mateid}))
  {
	  #print "exists $mateid\n";
    if($2 eq "1")
    {    
      $splice_read{$readid}=$splice_read_ini{$readid}; 
      #print "record mate_bk $mateid==$splice_read_ini{$mateid}\n";
      $mate_bk{$mateid}=$splice_read_ini{$mateid};
      print TEST2 "$splice_read_ini{$readid}\n";
      print TEST2 "$splice_read_ini{$mateid}\n";            
    }
  }else
  {
   $splice_read{$readid}=$splice_read_ini{$readid};
  }
}
%splice_read_ini=();
#####keep read order by id######################mapped_sample is in order##########
foreach $readid (sort{ ($da=$a) =~s/_//;$da =~/\S+\.(\d+)/;$da=$1;($db=$b)=~s/_//;$db=~/\S+\.(\d+)/;$db=$1; $da<=>$db}keys %splice_read)
{
	#print "processing read $readid\n";
   	@b=split(/\t/,$splice_read{$readid});
  	($matestrand,$matechr,$matepos1,$matepos2,$mate_alignscore,$mate_type)=&read_mate($readid);#mate align pos and score
	print TEST3 "$readid#$matestrand#$matechr#$matepos1#$matepos2#$mate_alignscore#$mate_type\n";
  	###mate can not align####################################
  	if($matechr == -1)   
  	{
   	  print OUT1 "$splice_read{$readid}\n";
   	  next; ###read next read          
  	}
  ##################mate align#############################################
  	$fhasright=-1;
  	$fhasright2=-1;
  	$fhaswrong=-1;
  	%rightposnum=();
	@null_model=();
	$nulli=0;
  	@splice_pos=split(/:/,$b[2]);
	@mstrand=split(/\t/,$matestrand);
     	@mchr=split(/\t/,$matechr);
     	@mpos1=split(/\t/,$matepos1);
     	@mpos2=split(/\t/,$matepos2);
     	@mscore=split(/\t/,$mate_alignscore);
	@mtype=split(/\t/,$mate_type);
  	for($i=0;$i<=$#splice_pos;$i++)#
  	{
     		@c=split(/,/,$splice_pos[$i]);#contain junction read mismatch
     		$readscore=$c[4];#alignment score
		#mate map in line
     		for($mi=0;$mi<=$#mpos1;$mi++)#for each mate pos
     		{
      			if(($c[1] eq $mchr[$mi])and($c[0] ne $mstrand[$mi]))#same chr and reverse strand
      			{
       
        		   if(($c[2]<=$mpos1[$mi] and $c[3]>=$mpos2[$mi] and $mtype[$mi] eq "1")or($c[2]==$mpos1[$mi] and $c[3]==$mpos2[$mi] and $mtype[$mi] eq "2"))
        		   {
             			$fhasright+=1; 
             			if($#mpos1==0){$type="unique";}
             			else{$type="multiple";} 
        			if(!exists($rightposnum{$splice_pos[$i]}))
        			{
          				$fhasright2+=1;
          				$rightpos[$fhasright2]=$splice_pos[$i];
          				$rightposnum{$splice_pos[$i]}=1;#
	  				$rightscore[$fhasright2]=$splice_pos[$i].",".$mscore[$mi];#splice_pos contain readscore
        			}
        			else
        			{
	  			$rightscore[$fhasright2].=",".$mscore[$mi];#multiple mate pos support this one splice site
        			}        
        		   }
                                $splice_pos[$i]=~/\S+,\S+,\S+,\S+,(\S+)/;
				$null_model[$nulli]="$1,$mscore[$mi]";
				$nulli+=1;
			     
      			}#if
     		}#foreach mi
      	if($fhasright==-1)#
      	{
       		$fhaswrong+=1;
       		$wrongpos[$fhaswrong]=$splice_pos[$i];#
      	}	   
   	}#for splicepos
##################add linear null model##############################
if($fhasright>-1)
{
   if(exists($linear_read{$readid}))
   {
    my  @lineara=split(/\t/,$linear_read{$readid});
    my  @linear_pos=split(/:/,$lineara[2]); 
    for($i=0;$i<=$#linear_pos;$i++)#
    {
     	@c=split(/,/,$linear_pos[$i]);
     	$readscore=$c[4];#alignment score
     	
     	for($mi=0;$mi<=$#mpos1;$mi++)#for each mate pos
     	{
      	 if(($c[1] eq $mchr[$mi])and($c[0] ne $mstrand[$mi]))#same chr and reverse strand
      	 {       
          $linear_pos[$i]=~/\S+,\S+,\S+,\S+,(\S+)/;
	  $null_model[$nulli]="$1,$mscore[$mi]";
	  $nulli+=1;			     
      	 }#if
     	 }#foreach mi
    }
    print TEST "$linear_read{$readid}\n";
    }
}
#####################################################################
   if($fhasright>-1)#mate pos
   {
    $num=$fhasright2+1;# $fhasright+1 for test multiple mate support;
    print OUT2 "$readid\t$num\t";
    print RECORD3 "$readid\t";
    for($ri=0;$ri<$fhasright2;$ri++)
    {
      print OUT2 "$rightpos[$ri]:";
      print RECORD3 "$rightscore[$ri]:";
    }
    print OUT2 "$rightpos[$fhasright2]\n";
    print RECORD3 "$rightscore[$fhasright2]\n"; 
    for($nulli=0;$nulli<$#null_model;$nulli++)
    {
    print RECORD3 "$null_model[$nulli]\t";    
    }
    print RECORD3 "$null_model[$nulli]\n";
   }###has unless one right! in multiple pos and multiple mate pos
   else ######all pos is wrong is in OUT3
   { $num=$fhaswrong+1;
     print OUT3 "$readid\t$num\t";
     for($wi=0;$wi<$fhaswrong;$wi++)
     {
      print OUT3 "$wrongpos[$wi]:";
     }
     print OUT3 "$wrongpos[$fhaswrong]\n";
   }#else
 }#foreach

close IN1;
close Inmate;
close Inmate2;
close OUT1;
close OUT2;
close OUT3;
close RECORD3;
}

1;
