#!/usr/bin/env perl
use strict;
use Getopt::Long;
use Pod::Usage;
use call_peak qw(padjust_and_filter);

print "---------------------------------------\n";
print "Running in step filter_by_fdr ...\n";

#const
my $cutoff=0.01; 
my $help=0;
my $padjust_method="bonferroni";

GetOptions("cut-offs=f"		=> \$cutoff,
	   "padjust_method=s"  => \$padjust_method,
           "h|help" 		=> \$help) or pod2usage(-exitval => 2, -verbose => 2);

&padjust_and_filter($cutoff,$padjust_method);

print "Step filter_by_fdr over!\n";
print "---------------------------------------\n";


__END__

=head1 NAME

filter_by_fdr.pl

=head1 SYNOPSIS

filter_by_fdr.pl [options]

=head1 OPTIONS

=over 

=item B<-cut-offs> <float>    

FDR cut-offs (Default:0.05)

=item B<-padjust_method><string>

method of adjust p-values for multiple comparisons same as in R function p.adjust can be "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none" (default:"bonferroni")

=item B<-h|--help>   

Show help information.

=back

=head1 OUTPUT

This program will generate final circRNA junctions which satisfying the FDR cut-offs. In file "./output/peak_fdr_passcutoff.txt".

=cut


