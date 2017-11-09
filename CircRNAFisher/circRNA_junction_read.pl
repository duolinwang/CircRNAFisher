#!/usr/bin/env perl
use strict;

use Getopt::Long;
use Pod::Usage;

use get_linear_seq;
use get_mapped_linear_scores;
use get_backsplice_seq;
use unmapped_perl_utils;
use pre_paired_filter qw(get_lariat_read merge_hg19_trans);
use change2fastq;
use paired_end_filter;
use candidate_circ_perl_utils qw(order_splice_read); 
use change_transpos2hg19pos qw(change2hg19);
use find_transcripts qw(find_transcripts);
use delete_lariat_read_backsplice;
use print_parameter; 
use order_linear_read;

#const
my $anchorl=10;
my $p = 8;
my $unique=1;
my $mismatch_thres=3;
my $MQ=10;
my $range=15;
my $flariat=1;

my $refseq="";
my $bt2_dix="";
my $help=0;
my $bowtie_path="";
my $samtools="";
my $argnum=scalar(@ARGV);
my $known_exons="";

print "---------------------------------------\n";
print "Running in step circRNA_jounction_read...\n";
GetOptions(
           "refseq=s" 		=> \$refseq,
           "ep=s" 		=> \$known_exons,
           "re-anchorl=i"  	=> \$anchorl,
	   "merge-range=i"	=> \$range,
	   "lariat=i"			=> \$flariat,
           "m=i"  		=> \$mismatch_thres,
           "p|threads=i"	=> \$p,
           "xg=s"  		=> \$bt2_dix,
           "mq=f"	=> \$MQ,
           "bowtie2=s" 		=> \$bowtie_path,
           "samtools=s" 	=> \$samtools,
           "unique=i"           => \$unique,
           "h|help" 		=> \$help) or pod2usage(-exitval => 2, -verbose => 2);
pod2usage(-exitval => 2,-verbose => 2) if ($help == 1);
pod2usage(-msg => "Invalid number of arguments!", -exitval => 2, -verbose => 2) if ($argnum < 6);
pod2usage(-msg => "bowtie2 genome index is missing", -exitval => 2, -verbose => 2) if ($bt2_dix eq '');
pod2usage(-msg => "Invalid reference sequence!", -exitval => 2, -verbose => 2) if ($refseq eq '');
pod2usage(-msg => "Invalid known exons!", -exitval => 2, -verbose => 2) if ($known_exons eq '');
pod2usage(-msg => "lariat must be 0 or 1!", -exitval => 2, -verbose => 2) if ($flariat != 1 and $flariat!=0);

#pod2usage(-msg => "Invalid merge-range, please set it more than 0", -exitval => 2, -verbose => 2) if ($range <=0);

&print_parameter("lariat",$flariat);

open(IN,"./output/parameters.txt")or die "ERROR: haven't product parameters.txt, please goto the first step!";
my $line;
my $f1=0;
my $f2=0;
my $f3=0;
my $f4=0;
my $anchorl_f;
my $readlength=0;
my $librarysize=0;
my $phred=0;
while($line=<IN>)
{
 chomp($line);
 my @paras=split(/\t/,$line);
 if($paras[0] eq "anchorl")
 {
  $anchorl_f=$paras[1];
  $f1=1;
 }
 if($paras[0] eq "readlength"){$readlength=$paras[1];$f2=1;}
 if($paras[0] eq "library_size"){$librarysize=$paras[1];$f3=1;}
 if($paras[0] eq "phred"){$phred=$paras[1];$f4=1;}
}
close IN;

if($f1==0){die "ERROR: no anchorl parameter, something wrong in the find_candidate_circRNAs step !";}
if($f2==0){die "ERROR: no readlength parameter, something wrong in the find_unmapped_reads step !";}
if($f3==0){die "ERROR: no librarysize parameter, something wrong in the find_unmapped_reads step !";}
if($f4==0){die "ERROR: no phred parameter, something wrong in the find_unmapped_reads step !";}

pod2usage(-msg => "Sequence length of an anchor make no sense !",-exitval => 2,-verbose => 2) if ($anchorl > $anchorl_f);

if ($refseq ne ""){$refseq .="/";}

&get_linear_seq($refseq,$anchorl,$readlength);
&get_backsplice_seq($refseq,$anchorl,$readlength);

if ($bowtie_path ne "") { $bowtie_path .= "/"; }
my $command=$bowtie_path."bowtie2-build -f --quiet ./output/ref_linear_seq.fa ./output/refseq_linear";
&runcommand($command);
my $command=$bowtie_path."bowtie2-build -f --quiet ./output/ref_backsplice_seq.fa ./output/refseq_back";
&runcommand($command);
&change2fastq();
###############linear realign####################################
$command = $bowtie_path."bowtie2 -p$p --very-sensitive --mm --score-min=C,-15,0 -k 20 -q --rfg 100,100 --rdg 100,100 -x ./output/refseq_linear -U ./output/unmapped.fastq -S ./output/sample_linear.sam 2> ./output/bowtie2_realign_linear.log";
  print "realign unmapped reads to reflinear...\n";
   &runcommand($command);
     system("cat ./output/bowtie2_realign_linear.log");
    if ($samtools ne "") { $samtools .= "/"; }
   my $samtoolspath="";  
   my $samtoolspath = $samtools."samtools";
   $command="$samtoolspath view -hubS ./output/sample_linear.sam | $samtoolspath sort -n - ./output/sample_linear";
   &runcommand($command);
   $command="rm ./output/sample_linear.sam";
   &runcommand($command);
   $command = "$samtoolspath view -F 4 ./output/sample_linear.bam > ./output/re_mapped_linear.sam";
   &runcommand($command);
   $command="rm ./output/sample_linear.bam";
   &runcommand($command);
 ###################################################
 ###################################################
 $command = $bowtie_path."bowtie2 -p$p --very-sensitive --mm --score-min=C,-15,0 -k 20 -q --rfg 100,100 --rdg 100,100 -x ./output/refseq_back -U ./output/unmapped.fastq -S ./output/sample.sam 2> ./output/bowtie2_realign.log";
 print "realign unmapped reads to genome...\n";
    &runcommand($command);
    system("cat ./output/bowtie2_realign.log");

   if ($samtools ne "") { $samtools .= "/"; }

   my $samtoolspath="";  
  
   my $samtoolspath = $samtools."samtools";
   $command="$samtoolspath view -hubS ./output/sample.sam | $samtoolspath sort -n - ./output/sample";
   &runcommand($command);
   $command="rm ./output/sample.sam";
   &runcommand($command);
   $command = "$samtoolspath view -F 4 ./output/sample.bam > ./output/re_mapped.sam";
   &runcommand($command);
   $command="rm ./output/sample.bam";
   &runcommand($command);
   ######################################
   ###if($range eq ""){$range=$readlength-2*$anchorl}
   ###&merge_read_backsplice($range); 
   &get_mapped_linear_scores($readlength,$anchorl);
   &order_linear_read();  
   &get_lariat_read($readlength,$anchorl);
   &order_splice_read("re_");  
   &change2hg19($known_exons);    
   &merge_hg19_trans($readlength);
   &paired_end_filter($known_exons,$unique,"re_",$MQ,$readlength); 
   ######&merge_backsplice($range);  
   &find_transcripts($known_exons);
   &delete_lariat_read_backsplice(); 
print "Step circRNA_jounction_read over!\n";
print "---------------------------------------\n";

__END__

=head1 NAME

circRNA_jounction_read.pl

=head1 SYNOPSIS

circRNA_jounction_read.pl [options]-xg <bt2-idx> -refseq <seq> -ep <knownexons>

=head1 ARGUMENTS

=over

=item B<-xg> <bt2-idx>

bowtie2 index filename prefix for genome alignment (minus trailing .X.bt2).

=item B<-refseq> <folder>

the folder contains fasta files of reference genome, one file for each chromosome, filename prefix is its corresponding chromosome, end with ".fa", for instance chr1.fa

=item B<-ep> <knownexons>

File with annotated exon information.

=back

=head1 OPTIONS

=over 

=item B<-re-anchorl> <int>   

Sequence length of an anchor for circRNA_jounction_read procedure (Default:5, don't set it longer than the length of an anchor set in last step, recommand 5-10)

=item B<-merge-range> <int>   

CircRNA merging range (Default:15)

=item B<-m> <int>

maximum acceptable mismatch number of a read (Default:3)

=item B<-unique> <int>

whether only keep uniquely mapped reads. 1 for uniquely mapped reads, other values for all mapped reads (Default:1)

=item B<-lariat> <int>

whether delete intron lariat. 1 for delete, 0 for not delete (Default:1)

=item B<-bowtie2> <path> 

The path to the Bowtie executables.(Default: the path to Bowtie executables is assumed to be in the user's PATH environment variable)

=item B<-samtools> <path> 

The path to the Samtools executables.(the path to Bowtie executables is assumed to be in the user's PATH environment variable)

=item B<-mq> <float>

Mapping quality threshold to filter circular junctions, ranges from 0 to 255 (Default:10)

=item B<-p/--threads> <int> 

number of alignment threads to launch (Default: 8)

=item B<-h|--help>   

Show help information.

=back

=head1 OUTPUT

This program will generate more exact counts for junction reads in file "re_unique_spliceplot_withstand.txt".

=cut

