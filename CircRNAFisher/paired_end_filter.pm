#!/usr/bin/env perl

package paired_end_filter;
use strict;

use filter_splice_site_by_pairedend;
use get_unique_reads;
use get_splice_site_plot_strand;
use filter_candidate_junction_bymaq;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(paired_end_filter);
our @EXPORT_OK = qw(paired_end_filter);


sub paired_end_filter
{
	print "Paired end reads filter...\n";
	my ($known_exons,$unique,$reflag,$maq,$readlength)=@_;
	&filter_splice_site_by_pairedend($reflag);
        &filter_candidate_junction_bymaq($maq);	
	if($unique and $reflag eq "re_")
        {
	 &get_unique_reads($reflag);
	}
	&get_splice_site_plot_strand($unique,$readlength);
}

1;
