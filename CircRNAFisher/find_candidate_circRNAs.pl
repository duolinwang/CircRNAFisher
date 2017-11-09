#!/usr/bin/env perl
use strict;
use Getopt::Long;
use Pod::Usage;
use candidate_circ_perl_utils qw(unmapped_2_anchors find_circ extension extension_line order_splice_read);
use unmapped_perl_utils;
use print_parameter;
use get_splice_site;

#const
my $anchorl = 20;
my $anchormis =3;
my $mismatch_thres = 3;#
my $misextend = 3;
my $find_range = 2000000;#circRNA detected range 
my $p = 8;
my $refseq = "";
my $bt2_dix = "";
my $help = 0;
my $bowtie_path = "";
my $samtools = "";

my $argnum=scalar(@ARGV);

print "---------------------------------------\n";
print "Running in step find_candidate_circRNAs ...\n";

GetOptions("refseq=s" 		=> \$refseq,
           "anchorl=i"  	=> \$anchorl,
           "anchormis=i" 	=> \$anchormis,
           "m=i"  		=> \$mismatch_thres,
           "em=i"               => \$misextend,
           "p|threads=i"	=> \$p,
           "xg=s"  		=> \$bt2_dix,
	   "bowtie2=s" 		=> \$bowtie_path,
           "samtools=s" 	=> \$samtools,
           "circ_range=i"       => \$find_range,
           "h|help" 		=> \$help) or pod2usage(-exitval => 2, -verbose => 2);
pod2usage(-exitval => 2,-verbose => 2) if ($help == 1);
pod2usage(-msg => "Invalid number of arguments!", -exitval => 2, -verbose => 2) if ($argnum < 4);
pod2usage(-msg => "Invalid reference sequence!", -exitval => 2, -verbose => 2) if ($refseq eq '');
pod2usage(-msg => "bowtie2 genome index is missing", -exitval => 2, -verbose => 2) if ($bt2_dix eq '');
pod2usage(-msg => "Beyond the acceptable mismatch number,please set it no more than $mismatch_thres.", -exitval => 2, -verbose => 2) if ($misextend>$mismatch_thres);
pod2usage(-msg => "Invalid length of anchor: too short, please set it no more than 10!", -exitval => 2, -verbose => 2) if ($anchorl <=10);

&print_parameter("anchorl",$anchorl);#write anchorl into parameters.txt

my $stat=&unmapped_2_anchors($anchorl);

pod2usage(-msg => "Invalid length of anchor: beyond the length of read!", -exitval => 2, -verbose => 2) if ($stat == 2);

pod2usage(-msg => "Invalid acceptable mismatch number of anchors, too big! The anchor is only $anchorl length,defualt 3 for 20 length anchors", -exitval => 2, -verbose => 2) if ($anchormis>=$anchorl/5);

pod2usage(-msg => "Invalid acceptable mismatch number of anchors, too big, please set it no more than $mismatch_thres", -exitval => 2, -verbose => 2) if ($anchormis>$mismatch_thres);

my $command = "";

if ($bowtie_path ne "") { $bowtie_path .= "/"; }
    my $score=-6*$anchormis;
    $command = $bowtie_path."bowtie2 -p$p --reorder --mm --score-min=C,$score,0 --rdg 20,20 --rfg 20,20 -k 20 -q --ignore-quals -x $bt2_dix -U ./output/sample_anchors.fastq -S ./output/realign_anchors_k20.sam 2> ./output/bowtie2_realign_anchors_k20.log";
   print "align anchors to genome...\n";
   &runcommand($command);
   system("cat ./output/bowtie2_realign_anchors_k20.log");

   if ($samtools ne "") { $samtools .= "/"; }

   my $samtoolspath="";  
  
   my $samtoolspath = $samtools."samtools";
   $command = "$samtoolspath view -F 4 -S ./output/realign_anchors_k20.sam > ./output/mapped_realign_anchors_k20.sam";
   &runcommand($command);
   $command = "$samtoolspath view -f 4 -S ./output/realign_anchors_k20.sam > ./output/unmapped_realign_anchors_k20.sam";
   &runcommand($command);

   $command="rm ./output/realign_anchors_k20.sam";
   # &runcommand($command);

  &find_circ($find_range,$anchorl);
  if ($refseq ne ""){ $refseq .="/"; }
  &extension_line($refseq,$misextend,$mismatch_thres,$anchorl);
  &extension($refseq,$misextend,$mismatch_thres,$anchorl);
  &order_splice_read("");
  &get_splice_site();

print "Step find_candidate_circRNAs over!\n";
print "---------------------------------------\n";

__END__

=head1 NAME

find_candidate_circRNAs.pl

=head1 SYNOPSIS

./find_candidate_circRNAs.pl [options] -xg <bt2-idx> -refseq <seq>

=head1 ARGUMENTS

=over

=item B<-xg> <bt2-idx>

bowtie2 index filename prefix for genome alignment (minus trailing .X.bt2).

=item B<-refseq> <folder>

the folder contains fasta files of reference genome, one file for each chromosome, filename prefix is its corresponding chromosome, end with ".fa", for instance chr1.fa

=back

=head1 OPTIONS

=over 

=item B<-anchorl> <int>   

Sequence length of an anchor (Default:20, better leave as default)

=item B<-m> <int>

maximum acceptable mismatch number of a read (Default:3)

=item B<-anchormis> <int> 

maximum acceptable mismatch number of anchors (Default:3) 

=item B<-em> <int>

maximum acceptable mismatch number during extension procedure in each direction (Default:3, please don't set it more than -m <int>)

=item B<-circ_range> <int>

largest acceptable distance on genome between the two splice sites of circRNA (Default:2000000)

=item B<-unique> <int>

whether only keep uniquely mapped reads. 1 for uniquely mapped reads, other values for all mapped reads (Default:1)

=item B<-bowtie2> <path> 

The path to the Bowtie executables.(Default: the path to Bowtie executables is assumed to be in the user's PATH environment variable)

=item B<-samtools> <path> 

The path to the Samtools executables.(the path to Bowtie executables is assumed to be in the user's PATH environment variable)

=item B<-p/--threads> <int> 

number of alignment threads to launch (Default: 8)

=item B<-h|--help>   

Show help information.

=back

=head1 OUTPUT

This program will generate candidate circRNA splice sites in file "unique_spliceplot_withstand.txt".

=cut






