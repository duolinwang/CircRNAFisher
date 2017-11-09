#!/usr/bin/env perl

package print_parameter;
use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(print_parameter);
our @EXPORT_OK = qw(print_parameter);

sub print_parameter{
my ($parameter,$value)=@_;
my $find=0;
print "update parameter $parameter=$value\n";
open(IN,"+<./output/parameters.txt")or die "can't open file \"./output/parameters.txt\" : $!";
my $line;
my $inpos=tell(IN);
while($line=<IN>)
{
  if($line =~/$parameter/)
  {
    $find=1;    
    seek(IN,$inpos,0); 
    print IN "$parameter\t$value\n";
    <IN>;
  }
  $inpos=tell(IN);
}
if($find==0){print IN "$parameter\t$value\n";}
close IN;
}
