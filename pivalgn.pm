#!/usr/bin/perl -w
#
# pivalgn.pm: Implements the pivalign lengthening and aligning algorithm. 
#   K-mers are aligned using the least probable as a pivot.
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
use List::Util qw(min);
use Text::LevenshteinXS qw(distance);

# get index from pivot arrays
sub get_index {
    my ($valsl, $valsi, $valsr) = @_;
    my $min = min @$valsl, @$valsi, @$valsr;

    for my $i (0..$#{$valsi}) {
	return $i if ($valsi->[$i] == $min);
    }

    for my $i (0..$#{$valsr}) {
	return $#{$valsi} + 1 + $i if ($valsr->[$i] == $min);
    }

    for my $i (0..$#{$valsl}) {
	return -$i -1 if ($valsl->[$i] == $min);
    }    
}

# find pivot for pivalgn 
sub find_pivot {
    my ($left, $right) = @_; # !! left must be LONGER or EQUAL !!
    my ($lenl, $lenr) = (length $left, length $right);
    my ($mdev, $end) = ($lenr - 3, 0);
    my (@valsl, @valsi, @valsr);

    for (my $dev = 1; $dev <= $mdev; $dev++) {
	my $cmpl = substr $left, 0, $lenr - $dev;
	my $cmpr = substr $right, $dev;
	push @valsl, distance $cmpl, $cmpr;
    }
    for (my $diff = $lenl - $lenr; $diff >= 0; $diff--, $end++) {
	my $comp = substr $left, $end, $lenr;
	push @valsi, distance $comp, $right;
    } 
    for (my $dev = 1; $dev <= $mdev; $dev++) {
	my $cmpl = substr $left, -($lenr - $dev);
	my $cmpr = substr $right, 0, $lenr - $dev;
	push @valsr, distance $cmpl, $cmpr;
    }
    return get_index \@valsl, \@valsi, \@valsr;
}

# lengthen k-mers on either side to align them
sub pa_lengthen_k_mers {
    my ($left, $right, $file) = @_;
    my $pivot = find_pivot $left, $right;

    if ($pivot > 0) {
	if ($pivot + length $right > length $left) {
	    for (my $i = $pivot; $i--; ) {
		$right = add_base_on LEFT, $right, $file;
	    }
	    for (my $i = length($right) - length($left); $i > 0 ; $i--) {
		$left = add_base_on RGHT, $left, $file;
	    }
	}
	else {
	    for (my $i = $pivot; $i--; ) {
		$right = add_base_on LEFT, $right, $file;
	    }
	    for (my $i = length($left) - length($right); $i > 0 ; $i--) {
		$right = add_base_on RGHT, $right, $file;
	    }
	}
    } 
    else {
	for (my $i = $pivot; $i++; ) {
	    $left = add_base_on LEFT, $left, $file;
	}
	for (my $i = length($left) - length($right); $i > 0 ; $i--) {
		$right = add_base_on RGHT, $right, $file;
	}
    }
    return ($left, $right, $pivot);
}

# align k-mers according to the pivalgn algorithm
sub pa_align_k_mers {
    my ($tbl_ref, $file) = @_;
    my %ret = %$tbl_ref;
    my (@lengthen, @temp);

    my @k_mers = sort {$tbl_ref->{$a} <=> $tbl_ref->{$b} } keys %$tbl_ref;
    my $main_k_mer = shift @k_mers;
    foreach my $k_mer (@k_mers) {
	my ($left, $right, $pivot) = 
	    pa_lengthen_k_mers $main_k_mer, $k_mer, $file;
	$ret{$left}  = delete $ret{$main_k_mer};
	$ret{$right} = delete $ret{$k_mer};

	if (length $main_k_mer < length $left) {
	    my $diff = length($left) - length($main_k_mer);
	    foreach my $k_mer (@lengthen) {
		my $new = $k_mer;
		for (my $i = $diff; $i--; ) {
		    $new = add_base_on( ($pivot > 0)?RGHT:LEFT, $new, $file);
		}
		$ret{$new} = delete $ret{$k_mer};
		push @temp, $new;
	    }
	    @lengthen = @temp;
	    undef @temp;
	}
	push @lengthen, $right;
	$main_k_mer = $left;
    }
    return %ret;
}

1
