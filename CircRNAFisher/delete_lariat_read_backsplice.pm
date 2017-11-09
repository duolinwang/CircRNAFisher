
#!/usr/bin/env perl

package delete_lariat_read_backsplice; 


require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(delete_lariat_read_backsplice);
our @EXPORT_OK = qw(delete_lariat_read_backsplice);

sub delete_lariat_read_backsplice
{
open(IN1,"./output/lariat_read_backsplice")or die "can not open file:./output/lariat_read_backsplice\n";
$in2="./output/re_spliceplot_withstand.txt";
open(IN2,$in2)or die "can not open file $in2 \n";
open(IN3,"./output/lariat2.txt")or die "can not open file :./output/lariat2.txt\n";
open(OUT,">./output/re_spliceplot_withstrand_detetelariat.txt")or die "$!";

while($line=<IN1>)
{
chomp($line);
@a=split(/\t/,$line);
if(!exists($lariat_back{$a[1]}))
{
$lariat_back{$a[1]}=1;
}else{$lariat_back{$a[1]}+=1;}
}
close IN1;

while($line=<IN3>)
{
@a=split(/\t/,$line);
$back_lariat{$a[0]}=1;
}
close IN3;

while($line=<IN2>)
{
chomp($line);
@a=split(/\t/,$line);
if(!exists($back_lariat{$a[0]}))
{
  if(!exists($lariat_back{$a[0]}))
  {
  print OUT "$line\n";
  }
  else
  {
    $ri=$lariat_back{$a[0]}/$a[1];
    if($ri<0.8)
    {
     print OUT "$line\n";
    }
  }
}
}

}


#&delete_lariat_read_backsplice();
1;

