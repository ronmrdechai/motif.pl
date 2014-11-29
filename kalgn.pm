#!/usr/bin/perl -w
#
# kalgn.pm: Implements the 4algn lengthening and alignment algorithm. K-mers
#   are lengthened according to smaller 4-mers indexed within them
#
# 2013 (C) Ron Mordechai
#
# License: The Perl Artistic License
# <http://dev.perl.org/licenses/artistic.html>
#
# THIS PACKAGE IS PROVIDED "AS IS"  AND WITHOUT ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
# MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#

use strict;
use Text::LevenshteinXS qw(distance);

# allow using Smartmatch without those annoying warnings
no if $] >= 5.017011, warnings => 'experimental::smartmatch';

# index with levenshtein distance of n or less
sub lev_index($$$) {
    my ($str, $s, $n) = @_;
    my ($lenstr, $lens) = (length($str), length($s));
    my $temp = '';

    for (my $i = 0; $i < $lenstr - $lens + 1; $i++) {
	$temp = substr $str, $i, $lens;
	return $i if ( $n >= distance $temp, $s);
    }
    return -1;
}

# find the index all k-mers should be aligned to
sub get_avg_index {
    my $tbl_ref = shift;
    my ($main_k_mer, $ret, $c) = (get_main_k_mer($tbl_ref), 0, 0);
    foreach my $key (keys %$tbl_ref) {
    	if (index($key, $main_k_mer) > -1) {
    	    $c++;
    	    $ret = $ret + index $key, $main_k_mer;
    	}
    }
    return int $ret / $c;
}

# delete false hits in the table before lengthening
sub remove_false_hits {
    my ($tbl_ref, $avg_index) = @_;
    my $main_k_mer = get_main_k_mer $tbl_ref;
    my ($c, $d) = (0, scalar(keys %$tbl_ref) );
    my %ret = %$tbl_ref;

    foreach my $key (keys %ret) {
	delete $ret{$key} if (
	    abs(lev_index($key, $main_k_mer, 1) - $avg_index) > 4);
    }

    foreach my $key (keys %ret) {
	$c++ if (lev_index($key, $main_k_mer, 1) == $avg_index);
    }
    if ( ($d - $c) < 4) {
	foreach my $key (keys %ret) {
	    delete $ret{$key} if (
		lev_index($key, $main_k_mer, 1) != $avg_index);
	}
    }
    return %ret;
}

# lengthen k-mers on either side to align them
sub ka_lengthen_k_mers {
    my ($tbl_ref, $avg_index, $file) = @_;
    my ($main_k_mer, %ret) = (get_main_k_mer($tbl_ref), %$tbl_ref);

    foreach my $key (keys %ret) {
	if (lev_index($key, $main_k_mer, 1) <= $avg_index) {
	    $ret{add_base_on LEFT, $key, $file} = delete $ret{$key};
	}
	else {
	    $ret{add_base_on RGHT, $key, $file} = delete $ret{$key};
	}
    }
    return (scalar(keys %ret) > 1) ? %ret : %$tbl_ref;
}

# push k-mers one shift into alignment
sub ka_align_once {
    my ($tbl_ref, $file) = @_;
    my $main_k_mer = get_main_k_mer $tbl_ref;

    my $avg_index = get_avg_index $tbl_ref;
    my %ret = remove_false_hits $tbl_ref, $avg_index;

    foreach my $key (keys %ret) {
	if (lev_index($key, $main_k_mer, 1) != $avg_index) {
	    %ret = ka_lengthen_k_mers \%ret, $avg_index, $file;
	    last;
	}
    }
    return (scalar(keys %ret) > 1) ? %ret : %$tbl_ref; # return %ret?
}

# align k-mers according to the 4algn algorithm
sub ka_align_k_mers {
    my ($tbl_ref, $file) = @_;
    my %aligned_k_mers = ka_align_once $tbl_ref, $file;
    my %aligned_old;
# until the two hashes match, or until we're left with one k-mer:
    until ( %aligned_k_mers ~~ %aligned_old and
	    [sort values %aligned_k_mers] ~~ [sort values %aligned_old] or
	    scalar(keys %aligned_k_mers) <=1 ) {
	%aligned_old = %aligned_k_mers;
	%aligned_k_mers = ka_align_once \%aligned_old, $file;
    }
    return %aligned_k_mers;
}

1
