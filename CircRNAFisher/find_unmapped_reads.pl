#!/usr/bin/env perl
use strict;
use Getopt::Long;
use Pod::Usage;
use unmapped_perl_utils;
use print_parameter; 
#const
my $p=8;

my $readfile1="";
my $readfile2="";
my $bt2_dix="";
my $bt2_trans="";
my $help=0;
my $bowtie_path="";
my $samtools="";
my $argnum=scalar(@ARGV);
my $phred;
print "---------------------------------------\n";
print "Running in step find_unmapped_reads...\n";

GetOptions("1=s" 		=> \$readfile1,
           "2=s" 		=> \$readfile2,
           "xg=s"  		=> \$bt2_dix,
           "xt=s"  		=> \$bt2_trans,
           "bowtie2=s" 		=> \$bowtie_path,
           "samtools=s" 	=> \$samtools,
           "p|threads=i"	=> \$p,
	   "phred=s"            => \$phred,
           "h|help" 		=> \$help) or pod2usage(-exitval => 2, -verbose => 2);
pod2usage(-exitval => 2,-verbose => 2) if ($help == 1);

pod2usage(-msg => "Invalid number of arguments!", -exitval => 2, -verbose => 2) if ($argnum < 10);
pod2usage(-msg => "read1 is missing", -exitval => 2, -verbose => 2) if ($readfile1 eq '');
pod2usage(-msg => "read2 is missing", -exitval => 2, -verbose => 2) if ($readfile2 eq '');
pod2usage(-msg => "Must specify the Phred encoding scheme as --phred33 or --phred64", -exitval => 2, -verbose => 2) if ($phred eq '');
pod2usage(-msg => "bowtie2 genome index is missing", -exitval => 2, -verbose => 2) if ($bt2_dix eq '');
pod2usage(-msg => "bowtie2 transcripts index is missing", -exitval => 2, -verbose => 2) if ($bt2_trans eq '');

if(system("ls -l output 2>/dev/null >/dev/null"))
{
print "mkdir output\n";
if(system("mkdir output")!=0){die "$!";}
}

else{print "Has output already ! Now cover it !\n";}

if(system("touch ./output/parameters.txt")!=0){die "$!";}

&print_parameter("phred",$phred);

if($readfile1 ne '' and $readfile2 ne '')
{
 &add_num($readfile1,$readfile2);
}

my $command = "";
if ($bowtie_path ne "") { $bowtie_path .= "/"; }
    $command = $bowtie_path."bowtie2 -p$p $phred --very-sensitive --mm --score-min=C,-15,0 -k 20 -q -x $bt2_dix -U ./output/merged_seq.fastq -S ./output/sample_sam_k20.sam 2> ./output/bowtie2_k20.log";
    print "align to genome...\n";
   &runcommand($command);
   system("cat ./output/bowtie2_k20.log");

if ($samtools ne "") { $samtools .= "/"; }
my $samtoolspath="";    
my $samtoolspath = $samtools."samtools";
    $command = "$samtoolspath view -hubS ./output/sample_sam_k20.sam | $samtoolspath sort -n - ./output/sample_bam_k20" ;
    &runcommand($command);
    $command="rm ./output/sample_sam_k20.sam";#delete sam
    &runcommand($command);
    $command = "$samtoolspath view -F 4 ./output/sample_bam_k20.bam > ./output/mapped_sample_k20.sam";
    &runcommand($command);
    $command = "$samtoolspath view -f 4 ./output/sample_bam_k20.bam > ./output/unmapped_sample_k20.sam";
    &runcommand($command); #delete bam
    $command="rm ./output/sample_bam_k20.bam";
    &runcommand($command);

    $command = $bowtie_path."bowtie2 -p$p $phred --very-sensitive --mm --score-min=C,-15,0 -k 20 -q -x $bt2_trans -U ./output/merged_seq.fastq -S ./output/sample_vs_transcript_temp.sam 2> ./output/bowtie2_trans_12.log";
    print "align to transcripts...\n";
    &runcommand($command);
    system("cat ./output/bowtie2_trans_12.log");
    $command = "$samtoolspath view -hbuS ./output/sample_vs_transcript_temp.sam| $samtoolspath sort -n - ./output/sample_vs_transcript";
    &runcommand($command);
    $command="rm ./output/sample_vs_transcript_temp.sam";#delete sam
    &runcommand($command);
    $command = "$samtoolspath view -f 4 ./output/sample_vs_transcript.bam > ./output/unmapped_sample_transcript.sam";
    &runcommand($command);
    $command = "$samtoolspath view -F 4 ./output/sample_vs_transcript.bam > ./output/mapped_sample_transcript.sam";
    &runcommand($command);
    $command="rm ./output/sample_vs_transcript.bam";#delete bam 
    &runcommand($command);
    &merge_unmapped_reads("./output/unmapped_sample_k20.sam","./output/unmapped_sample_transcript.sam");

    $command="rm ./output/unmapped_sample_k20.sam"; #delete unmapped
    &runcommand($command);
    $command="rm ./output/unmapped_sample_transcript.sam";#delete unmapped
    &runcommand($command);

print "---------------------------------------\n";
print "Step find_unmapped_reads over!\n";

__END__

=head1 NAME

find_unmapped_reads.pl

=head1 SYNOPSIS

./find_unmapped_reads.pl [options] -1 <read1> -2 <read2> -phred <--phred33|--phred64> -xg <bt2-idx> -xt <bt2-idx-trans>

=head1 ARGUMENTS

=over

=item B<-1> <read1>

file with #1 mates, paired with filte in <read2>, must be in fastq format, and same length or will cause issues

=item B<-2> <read2>

file with #2 mates, paired with filte in <read1>, must be in fastq format, and same length or will cause issues

=item B<-phred> <--phred33|phred64>

it only accepts "phred+33" and "phred+64" encoding qualities, same as bowtie2 --phred33 and --phred64, must specify which encoding system to be used. 

=item B<-xg> <bt2-idx>

bowtie2 index filename prefix for genome alignment (minus trailing .X.bt2).

=item B<-xt> <bt2-idx-trans>

bowtie2 index filename prefix for transcripts alignment (minus trailling .X.bt2).

=back

=head1 OPTIONS

=over 
	   
=item B<-bowtie2> <bowtie2_path> 

The path to the Bowtie executables.(Default: the path to Bowtie executables is assumed to be in the user's PATH environment variable)

=item B<-samtools> <samtools_path> 

The path to the Samtools executables.(the path to Bowtie executables is assumed to be in the user's PATH environment variable)

=item B<-p/--threads> <int> 

number of alignment threads to launch (Default: 8)

=item B<-h|--help>   

Show help information.

=back

=head1 OUTPUT

This program will generate 'unmapped_sample.sam'.

=cut



