#!/usr/bin/perl -w
#
# gxw_to_seq.pl: turn gxw file to regex sequence
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

print (<DATA>) && exit 1 if (@ARGV != 2);

my $reverse = $ARGV[1];

open my $file, "<", "$ARGV[0]"
    or die "Couldn't open $ARGV[0]: '$!'";

# extract A, C, G, and T values from the gxw file
my @arr;
my $i = 0;
while (<$file>) {
    chomp;
    if (/<Position Weights=/) {
	my %tmp;
	$_ =~ s/<Position Weights="//;
	$_ =~ s/;"><\/Position>//;
	( $tmp{A}, $tmp{C}, $tmp{G}, $tmp{T} ) = split(/;/);
	$arr[$i++] = \%tmp;
    }
}
close $file;

# generate regex sequence
my $seq = "";
foreach my $hash (@arr) {
    my @bases;
    while ( my ($key, $val) = each(%$hash) ) {
	push @bases, $key if ($val > 0);
    }

    if (scalar(@bases) == 1) {
	$seq = $seq . $bases[0];
	next;
    }

    $seq = $seq . "[";
    $seq = $seq . "$_" foreach (@bases);
    $seq = $seq . "]";
}

# reverse the sequence pairs (optional)
my $rev = "";
if ($reverse) {
    foreach my $char ( split(//, $seq ) ) {
	$rev = $rev . "A"   if ($char eq "T");
	$rev = $rev . "C"   if ($char eq "G");
	$rev = $rev . "G"   if ($char eq "C");
	$rev = $rev . "T"   if ($char eq "A");
	$rev = $rev . $char if ($char =~ /[^TGCA]/); # else
    }
}
else {
    $rev = $seq;
}

print $rev, "\n";
__DATA__
Usage: gxw_to_seq <file.gxw> <reverse: (1|0)>
