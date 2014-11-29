#!/usr/bin/perl -w
#
# gen_cncr_pos.pl: parse large cancer tables
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

print (<DATA>) && exit 1 if (@ARGV != 3);

my $N_name = "N_$ARGV[2]";
my $T_name = "T_$ARGV[2]";

open my $data, "<", $ARGV[0] or
    die "Couldn't open $ARGV[0]: '$?'";
open my $utr, "<", $ARGV[1] or
    die "Couldn't open $ARGV[1]: '$?'";
open my $out, ">", "$ARGV[2].stab" or
    die "Couldn't open $ARGV[2].stab: '$?'";

my $line = <$data>;
my @cancer = split "\t", $line;

# get indecies for cell type
my (@N_index, @T_index);
for (my $i = 0; $i < @cancer; $i++) {
    push @N_index, $i if ($cancer[$i] eq $N_name);
    push @T_index, $i if ($cancer[$i] eq $T_name);
}

# list all irregularly regulated genes
my @genes;
while (my $line = <$data>) {
    my @tmp = split "\t", $line;
    my $name = $tmp[0];
    my $next = $tmp[0];

    my ($N_avg, $N_div) = (0, 0);
    my ($T_avg, $T_div) = (0, 0);

    for (;;) {
	seek($data, -length($line), 1) && last if ($next ne $name);

	foreach my $i (@N_index) {
	    $N_avg = $N_avg + $tmp[$i];
	    $N_div++;
	}
	foreach my $i (@T_index) {
	    $T_avg = $T_avg + $tmp[$i];
	    $T_div++;
	}

	$line = <$data>;
	last unless($line);
	@tmp = split "\t", $line;
	$next = $tmp[0];
    }
    $N_avg = $N_avg / $N_div;
    $T_avg = $T_avg / $T_div;

    push @genes, $name if ( ($T_avg - $N_avg) <= -1 );
}

# save gene 3'UTRs to file
while (<$utr>) {
    my @tmp  = split "\t";
    my $name = $tmp[0];
    my $seq  = $tmp[1];
    foreach my $gene (@genes) {
	print {$out} $seq if ($gene eq $name);
    }
}
print "Generated positive file at $ARGV[2].stab\n"
__DATA__
Usage: gen_cncr_pos.sh <data.stab> <3utr.stab> <cell type>
