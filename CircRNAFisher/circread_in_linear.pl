open(IN1,"./output/mapped_sample_linear")or die "!!!";
open(IN2,"./output/splice_sites_read_mismatch")or die "!!!";
#open(IN2,"./output/test_splice_sites_read_mismatch")or die "!!!";
open(OUT,">./output/test_circread_inlinear");
open(OUT2,">./output/delete_line_splice_sites_read_mismatch");
our $Inpos=0;
sub ifline
{ 
my($find_id)=@_;
my $line;
my @a;
my $da,$db;
 ($da=$find_id)=~s/_//;$da =~/\S+\.(\d+)/;$da=$1; 
my $findk=20; #ps -k 15 so 20 is enough
my $findflag=0;
my %linear=();
while($line=<IN1>)
{ 
  chomp($line);
  @a=split(/\t/,$line);
  $lineid=$a[0];
  ($db=$lineid)=~s/_//;$db=~/\S+\.(\d+)/;$db=$1;
  if($da>$db){$Inpos=tell(IN1);}
  if($da<$db){last;} ###can't find findid
  if(($lineid eq $find_id) and $findk>=1)
  { 
      $findflag=1;
      $findk-=1;
      #$Inpos=tell(IN1);# if same don't update inpos! only da>db findid big update.
      if(!exists($linear{$find_id}))
      {
      $linear{$find_id}=$line;
      }else{$linear{$find_id}.=";".$line;}
  }
}
 seek(IN1,$Inpos,0); ###update mate pos 
 if($findflag==1)#find at least one
 { 
 return $linear{$find_id};
 }
 else
 { 
 return -1;
 }
}

while($line=<IN2>)
{
chomp($line);
@a=split(/\t/,$line);
$readid=$a[0];
$temp=&ifline($readid);
if($temp!=-1){print OUT "$line\t$temp\n";}else{print OUT2 "$line\n";}
}


