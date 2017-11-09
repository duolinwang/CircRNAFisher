#!/usr/bin/env perl
use strict;
use Getopt::Long;
use Pod::Usage;
use get_unique_mapped_reads;
use find_discordant_readpair;
use change2bed qw(get_backsplice_bed changemapped2bed);
use find_match_splice_site;
use get_background_length_count;
use call_peak qw(peak_count compute_pvalue);

#const
my $fragment=300; 
#my $max_fragmentsize=400;
my $unique=1;
my $help=0;
GetOptions("fragment=i"		=> \$fragment,
           "unique=i"           => \$unique,
           "h|help" 		=> \$help) or pod2usage(-exitval => 2, -verbose => 2);

print "---------------------------------------\n";
print "Running in step select_peak_circRNA ...\n";

open(IN,"./output/parameters.txt")or die "ERROR: haven't product parameters.txt, please goto the first step!";
my $line;
my $f1=0;
my $f2=0;
my $readlength=0;
my $flariat=0;
while($line=<IN>)
{
 chomp($line);	
 my @paras=split(/\t/,$line);
 if($paras[0] eq "readlength"){$readlength=$paras[1];$f1=1;}
 if($paras[0] eq "lariat"){$flariat=$paras[1];$f2=1;}
}
close IN;
if($f1==0){die "ERROR: no readlength parameter, something wrong in the find_unmapped_reads step !";}
if($f2==0){die "ERROR: no flariat parameter, run circRNA_junction_read.pl step first !";}


&find_discordant_readpair($unique,$readlength);
&get_backsplice_bed($flariat);
&changemapped2bed($unique);

my $command="";
   $command="cat ./output/mapped_hg19.bed ./output/backsplice_bed.bed > ./output/mappable.bed";
   system($command);
   $command="awk -F \'\\t\' \'{print \$1\",\"\$2\",\"\$3\",\"\$4\",\"\$5\",\"\$6}\' ./output/mappable.bed >./output/mappable_sorted.bed";   
   system($command); 
   $command="sort -t ',' -k1.3,1 -k2n -k3n ./output/mappable_sorted.bed -o ./output/mappable_sorted.bed";
   system($command); 

&find_match_splice_site($fragment,$readlength,$flariat);
&get_background_length_count($readlength);
&peak_count($flariat);
&compute_pvalue($fragment,$readlength);
print "Step select_peak_circRNA over!\n";
print "---------------------------------------\n";


__END__

=head1 NAME

select_peak_circRNA.pl

=head1 SYNOPSIS

select_peak_circRNA.pl [options]

=head1 OPTIONS

=over 

=item B<-fragment> <int>    

fragmentsize of paired-end reads (Default:300)

=item B<-unique> <int>

whether only keep uniquely mapped reads. 1 for uniquely mapped reads, other values for all mapped reads (Default:1)

=item B<-h|--help>   

Show help information.

=back

=head1 OUTPUT

This program will generate p-values and fdrs by call-peaking method. In file "peak_pvalue_qvalue".

=cut


