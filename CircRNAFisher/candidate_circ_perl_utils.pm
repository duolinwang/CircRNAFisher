#!/usr/bin/env perl

package candidate_circ_perl_utils;
use strict;
use unmapped_2_anchors;
use find_circ;
use extension;
use extension_line;
use order_splice_read;
use merge_backsplice;
use get_mapped_linear_scores;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(unmapped_2_anchors find_circ extension extension_line order_splice_read merge_backsplice get_mapped_linear_scores);
our @EXPORT_OK = qw(unmapped_2_anchors find_circ extension extension_line order_splice_read merge_backsplice get_mapped_linear_scores);

sub unmapped_2_anchors{
my $anchorl=$_[0];
return &unmapped_2_anchors::unmapped_2_anchors($anchorl);
}

sub find_circ
{
my ($find_range,$anchorl)=@_;
return &find_circ::find_circ($find_range,$anchorl);
}

sub extension_line
{
my ($refseq,$misextend,$mismatch_thres,$anchorl)=@_;
return &extension_line::extension_line($refseq,$misextend,$mismatch_thres,$anchorl);
}

sub extension
{
my ($refseq,$misextend,$mismatch_thres,$anchorl)=@_;
return &extension::extension($refseq,$misextend,$mismatch_thres,$anchorl);
}

sub order_splice_read
{
my $reflag=$_[0];
return &order_splice_read::order_splice_read($reflag);
}


sub merge_backsplice
{
my $range=$_[0];
return &merge_backsplice::merge_backsplice($range);
} 

sub get_mapped_linear_scores
{
return &get_mapped_linear_scores::get_mapped_linear_scores();
}
1;
