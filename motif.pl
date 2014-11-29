#!/usr/bin/env perl -w
#
# motif.pl: discriminatively search for a motif of a given length using
# given nucleotide sequences
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
use File::Path;
use Fcntl qw(:seek);
use GD::Graph::bars;
use List::Util qw(max min shuffle);
use Text::LevenshteinXS qw(distance);
use Getopt::Long qw(:config bundling);
use Math::GSL::Fit qw(gsl_fit_linear);
use Math::GSL::CDF qw(:hypergeometric);
use Math::GSL::Statistics qw(gsl_stats_tss);

sub LEFT () { 0 };
sub RGHT () { 1 };

our $VERSION = "0.1";

# get standard deviation of number array
sub stdev {
	my $data = shift;

	if (@$data == 1) {
		return 0;
	}

	my $avg = 0;
	foreach (@$data) {
		$avg = $avg + $_;
	}
	$avg = $avg / @$data;

	my $sqtotal = 0;
	foreach (@$data) {
		$sqtotal += ($avg - $_) ** 2;
	}
	my $std = ($sqtotal / (@$data - 1) ) ** 0.5;
	return $std;
}

# generate negative data from positive data by shuffling each line
sub generate_neg_file {
	my $pos = shift;
	my $rand = int rand 100;
	open my $pos_file, "<", "$pos";
	open my $neg_file, ">", ".$rand.dat";

	while (<$pos_file>) {
		my @line = split //;
		print {$neg_file} join('', shuffle(shuffle(shuffle(@line))) ), "\n";
	}
	close $neg_file;
	close $pos_file;
	return ".$rand.dat";
}

# generate a list of all words of length k
sub gen_k_mer_list {
	my $k = shift;
	my $bases = [qw(A T C G)];

	return @$bases if $k eq 1;

	my @prev = gen_k_mer_list ($k - 1);
	my @ret;

	for my $base (@$bases) {
		push @ret, map { $base . $_ } @prev;
	}
	return @ret;
}

# count all k-mers in the specified file
sub count_k_mers {
	my ($k, $file) = @_;
	my @k_mer_list = gen_k_mer_list $k;
	my %ret;
	foreach my $k_mer (@k_mer_list) {
		seek $file, SEEK_SET, 0;
		my $cnt = 0;
		while (<$file>) {
			if (/$k_mer/) {
				$cnt++;
			}
		}
		$ret{$k_mer} = $cnt;
	}
	return %ret;
}

# set up a table of k-mers and their p-values
sub compute_hypergeometric_p_value {
	my ($pos_tbl, $neg_tbl, $pop_size) = @_;
	my %ret;
	foreach my $k_mer (keys %$pos_tbl) {
		my $p_value =
		  gsl_cdf_hypergeometric_Q($pos_tbl->{$k_mer},
								   $pop_size, $pop_size,
								   $pos_tbl->{$k_mer} + $neg_tbl->{$k_mer});
		$ret{$k_mer} = $p_value;
	}
	return %ret;
}

# generate a -log10 histogram of p-values to determine threshold
sub gen_histogram {
	my $tbl_ref = shift;
	my (@values, @ret);
	foreach (values %$tbl_ref) {
		if ($_ == 0) {			# can't compute a log10 of 0
			push @values, 100;
		} else {
			push @values, int( -log($_) / log(10) );
		}
	}
	for (my $i = 0; $i < 20; $i++) {
		my @temp;
		foreach (@values) {
			push @temp, $_ if ($_ == $i + 1);
		}
		$ret[$i] = \@temp;
	}
	return @ret;
}

# REMOVE: print histogram data to file
sub gen_histogram_plaintext {
	my ($arr_ref, $out_dir) = @_;
	my @E = ( "01", "02", "03", "04", "05",
			  "06", "07", "08", "09", "10",
			  "11", "12", "13", "14", "15",
			  "16", "17", "18", "19", "20" );
	my @values;
	foreach (@$arr_ref) {
		push @values, scalar(@$_);
	}
	open my $file, ">", "$out_dir/histogram.txt" or
	  die "Couldn't open '$out_dir/histogram.txt': '$?'";
	for (my $i = 0; $i < 20; $i++) {
		print {$file} "$E[$i]\t$values[$i]\n";
	}
	print "Printed histogram data to $out_dir/histogram.txt\n";
	close $file;
}

# DEBUG: generate bar graph from histogram
sub gen_histogram_graph {
	my ($arr_ref, $out_dir) = @_;
	my @values;
	foreach (@$arr_ref) {
		push @values, scalar(@$_);
	}
	my @data = ( [ "e-01", "e-02", "e-03", "e-04", "e-05",
				   "e-06", "e-07", "e-08", "e-09", "e-10",
				   "e-11", "e-12", "e-13", "e-14", "e-15",
				   "e-16", "e-17", "e-18", "e-19", "e-20" ],
				 \@values
			   );
	my $graph = GD::Graph::bars->new(301,188);
	$graph->set(
				x_label           => 'Value',
				y_label           => 'Count',
				title             => 'P-Value Histogram',
				transparent       => 0,
				box_axis          => 0,
			   ) or die $graph->error;
    my $out = $graph->plot(\@data) or die $graph->error;
    open my $file, ">", "$out_dir/histogram.png" or
	  die "Couldn't open '$out_dir/histogram.png': '$?'";
    binmode $file;
    print {$file} $out->png;
    print "Created histogram at $out_dir/histogram.png\n";
    close $file;
}

# get p-value threshold using R squared - faster, less accurate?
sub get_coef_pval_threshold {
	my $arr_ref = shift;
	my ($x, $y);
	my @R2;
	for (my $i = 0; $i < scalar(@$arr_ref); $i++) {
		push @$y, scalar($arr_ref->[$i]);
		push @$x, $i + 1;
	}
	while (scalar(@$y) > 1) {
		my $ss_resid = ( gsl_fit_linear($x, 1, $y, 1, scalar @$y) )[-1];
		my $ss_total = gsl_stats_tss($y, 1, scalar @$y);
		push @R2, 1 - $ss_resid/$ss_total;
		shift @$y;
		shift @$x;
	}
	for (my $i = 0; $i < @R2; $i++) {
		return $i + 3 if ($R2[$i] > 0.7); # +3 because R2 cuts a little early
	}
	return -1;
}

# REMOVE: get p-value threshold using histogram diffrence - slow as fuck
sub get_hist_pval_threshold {
	my ($hista, $k, $pos, $neg) = @_;
	my @diff;

	my $pmt = generate_neg_file $pos;
	open my $pmt_file, "<", $pmt or die "Couldn't open $pmt: '$?'";
	open my $neg_file, "<", $neg or die "Couldn't open $neg: '$?'";

	my $lines; ++$lines while <$pmt_file>;
	my %pmt_tbl = count_k_mers $k, $pmt_file;
	my %neg_tbl = count_k_mers $k, $neg_file;

	close $pmt_file;
	close $neg_file;
	unlink $pmt;

	my %p_value_tbl =
	  compute_hypergeometric_p_value \%pmt_tbl, \%neg_tbl, $lines;
	print scalar(keys %p_value_tbl), "\n";
	my @histb = gen_histogram \%p_value_tbl;

	for (my $i = 0; $i < @$hista; $i++) {
		$diff[$i] = scalar( @{ $hista->[$i]} ) - scalar( $histb[$i] );
	}
	for (my $i = 0; $i < @diff; $i++) {
		return $i + 4 if ($diff[$i] < 0); # same as above
	}
	return -1;
}

# remove any k-mers over the p-value threshold
sub clean_up {
	my ($tbl_ref, $p) = @_;
	my %ret;
	while ( my ($key, $val) = each %$tbl_ref ) {
		$ret{$key} = $val if ($val < 10 ** -$p);
	}
	return %ret;
}

# DEBUG: keep only top n k-mers
sub no_thres_clean_up {
	my ($tbl_ref, $n) = @_;
	my %ret;
	my @sorted = sort {$tbl_ref->{$a} <=> $tbl_ref->{$b} } keys %$tbl_ref;
	for (my $i = 0; $i < $n; $i++) {
		$ret{$sorted[$i]} = $tbl_ref->{$sorted[$i]};
	}
	return %ret
}


# find the most common 4-mer in table (part of 4clus and 4algn)
sub get_main_k_mer {
	my $tbl_ref = shift;
	my @by_val = sort { $tbl_ref->{$a} <=> $tbl_ref->{$b} } keys %$tbl_ref;
	my ($main_k_mer, $k_mer, %ret) = ($by_val[0]);

	my ($len_full, $len) =
	  (length($main_k_mer), int(length($main_k_mer) * 2/3));

	for (my $i = 0; $i < $len_full - $len; $i++) {
		$k_mer = substr $main_k_mer, $i, $len;
		$ret{$k_mer} = 0;
		foreach (@by_val) {
			$ret{$k_mer}++ if ( -1 < index $_, $k_mer );
		}
	}
	my @ret_by_val = sort { $ret{$b} <=> $ret{$a} } keys %ret;
	return $ret_by_val[0];
}

# turn the hash table into a double hash matrix
sub table_to_matrix {
	my $tbl_ref = shift;
	my %ret;
	foreach my $row (keys %$tbl_ref) {
		my %temp;
		foreach my $col (keys %$tbl_ref) {
			$temp{$col} = 0;
		}
		$ret{$row} = \%temp;
	}
	return \%ret;
}

# sort k-mers according to levenshtein distance
sub sort_dist {
	my $tbl_ref = shift;
	my @ret;
	foreach my $row (keys %$tbl_ref) {
		foreach my $col (keys %{ $tbl_ref->{$row} } ) {
			$tbl_ref->{$row}->{$col} = distance $row, $col;
		}
	}
	while (%$tbl_ref) {
		my $row = (keys %$tbl_ref)[0];
		my %temp;
		foreach my $col (keys %{ $tbl_ref->{$row} } ) {
			if ($tbl_ref->{$row}->{$col} && $tbl_ref->{$row}->{$col} <= 2) {
				$temp{$col} = 1;
				delete $tbl_ref->{$col};
			}
		}
		$temp{$row} = 0;
		delete $tbl_ref->{$row};
		push @ret, \%temp;
	}
	return @ret;
}

# cluster sorted k-mers according to smaller 4-mers
sub cluster {
	my $tbl_ref = shift;
	my @groups = sort_dist table_to_matrix $tbl_ref;
	foreach my $left (@groups) {
		my $left_k_mer = get_main_k_mer $left;
		foreach my $right (@groups) {
			next if ($left == $right);
			my $right_k_mer = get_main_k_mer $right;
			my $dist = distance($left_k_mer, $right_k_mer);
			if ($dist <= 2) {
				@{$left}{keys %$right} = values %$right;
				my $i = 0; $i++ until ($groups[$i] == $right);
				splice @groups, $i, 1;
			}
		}
	}
	foreach (@groups) {			# copy the p-values over
		@{$_}{ keys %$_ } = @{$tbl_ref}{ keys %$_ };
	}
	return @groups;
}

# add base to k-mer according to the supplied file on either left or right
sub add_base_on {
	my ($side, $k_mer, $file) = @_;
	my ($A, $T, $C, $G, $max);
	my @lines;

	open my $fh, "<", "$file";
	while (<$fh>) {
		chomp;
		push @lines, $_;
	}
	close $fh;

	if ($side eq LEFT) {
		$A = grep /A($k_mer)/, @lines;
		$T = grep /T($k_mer)/, @lines;
		$C = grep /C($k_mer)/, @lines;
		$G = grep /G($k_mer)/, @lines;
		$max = max $A, $T, $C, $G;
		return ($max eq $A ? 'A' : $max eq $T ? 'T' : $max eq $C ? 'C' : 'G') .
		  $k_mer;
	} else {
		$A = grep /($k_mer)A/, @lines;
		$T = grep /($k_mer)T/, @lines;
		$C = grep /($k_mer)C/, @lines;
		$G = grep /($k_mer)G/, @lines;
		$max = max $A, $T, $C, $G;
		return $k_mer .
		  ($max eq $A ? 'A' : $max eq $T ? 'T' : $max eq $C ? 'C' : 'G');
	}
}

# make a position specific matrix out of the k-mer list
sub make_pssm_array {
	my $tbl_ref = shift;
	my ($k, $c) = (length ~~(keys %$tbl_ref)[0],
				   scalar keys %$tbl_ref);
	my @PSSM;

	for (my $i = 0; $i < $k; $i++) {
		$PSSM[$i] = { 'A' => 0, 'C' => 0, 'G' => 0, 'T' => 0 };
	}
	for (my $i = 0; $i < $k; $i++) {
		foreach (keys %$tbl_ref) {
			$PSSM[$i]{ substr $_, $i, 1 }++;
		}
		foreach ( qw(A C G T) ) {
			$PSSM[$i]{$_} = $PSSM[$i]{$_} / $c;
		}
	}
	return @PSSM;
}

# test to see if a pssm column has easily predictable data
sub test_low_entropy {
	my $tbl_ref = shift;
	my @arr;
	foreach ( qw(A C G T) ) {
		push @arr, $tbl_ref->{$_};
	}
	my $dev = stdev \@arr;
	if ($dev > 0.15) {
		return 0;
	}
	return 1;
}

# remove junk values added to the PSSM as a result of lengthening
sub clean_pssm_array {
    my $arr_ref = shift;
    my %a = %{${$arr_ref}[0]};

    while (test_low_entropy(\%a) ) {
		shift @$arr_ref;
		%a = %{${$arr_ref}[0]};
    }
    my %b = %{${$arr_ref}[-1]};

    while (test_low_entropy(\%b) ) {
		pop @$arr_ref;
		%b = %{${$arr_ref}[-1]};
    }
    return $arr_ref;
}

# finds if pssm matches string of neucleotides with a prob. of 50% or more
sub match_pssm {
	my ($PSSM_ref, $string) = @_;
	my $plen = scalar @$PSSM_ref;
	my $slen = length $string;

	for (my $i = 0; $i < $slen - $plen; $i++) {
		my @cmp = split(//, substr($string, $i, $plen));
		my $prob = 0;
		for (my $j = 0; $j < $plen; $j++) {
			$prob = $prob + $PSSM_ref->[$j]->{$cmp[$j]};
		}
		$prob = $prob / $plen;
		return 1 if ($prob > 0.5);
	}
	return 0;
}

# validate the cleaned positive data against the pssm to increase accuracy
sub validate_pos_dat {
	my ($PSSM_ref, $pos) = @_;
	my $rand = int rand 100;
	open my $pos_file, "<", "$pos";
	open my $neg_file, ">", ".$rand.dat";

	while (my $line = <$pos_file>) {
		if (match_pssm $PSSM_ref, $line) {
			print {$neg_file} $line;
		}
	}
	close $neg_file;
	close $pos_file;
	return ".$rand.dat";
}

# returns a textual representation of the pssm
sub pssm_to_label {
	my $arr_ref = shift;
	my $ret = "";
	foreach my $elmnt (@$arr_ref) {
		my $max_val = max(values %$elmnt);
		foreach my $key (keys %$elmnt) {
			if ($elmnt->{$key} == $max_val) {
				$ret = $ret . $key;
				last;
			}
		}
	}
	return $ret;
}

# make a logo out of a PSSM using gxw2logo
sub make_pssm_logo {
	my ($PSSM, $label, $out_dir) = @_;
	open my $tmp_file, ">", ".pssm.gxw" or
	  die "Couldn't open .pssm.gxw: '$?'";

	print { $tmp_file } "<WeightMatrix Name=\"$out_dir/$label\" " .
	  "Type=\"PositionSpecific\" Order=\"0\">\n";

	foreach (@$PSSM) {
		print { $tmp_file } "  <Position Weights=\"";
		foreach my $nuc ( qw(A C G T) ) {
			print { $tmp_file } "$_->{$nuc};";
		}
		print { $tmp_file } "\"></Position>\n";
	}
	print { $tmp_file } "</WeightMatrix>\n";
	qx(./lib/gxw2logo.pl -a A,C,G,T -B 2 < .pssm.gxw); # >/dev/null 2>&1);
	unlink ".pssm.gxw";
	print "Created PSSM logo at $out_dir/$label.png\n";
}

# load script arguments
sub load_args {
	my %ret = (
			   help        => '0',     histogram   => '0',
			   k           => 6,	lengthening => '4algn',
			   validate    => 0,	pos_file    => '',
			   neg_file    => '',	no_cluster  => 0,
			   out_dir     => '.',	delete_neg  => 0,
			   threshold   => undef,	no_thres    => undef,
			   version     => '',
			  );

	GetOptions (
				'help|h|?'            => \$ret{help},
				'histogram|H'         => \$ret{histogram},
				'start-count|k=i'     => \$ret{k},
				'lengthening|l=s'     => \$ret{lengthening},
				'validate|L'          => \$ret{validate},
				'neg-file|n=s'        => \$ret{neg_file},
				'no-cluster|N'        => \$ret{no_cluster},
				'out-dir|o=s'         => \$ret{out_dir},
				'force-threshold|t=i' => \$ret{threshold},
				'no-threshold|T=i'    => \$ret{no_thres},
				'version|V'           => \$ret{version} );

	if ($#ARGV > 0) {
		print "Recived ", join(', ', @ARGV),
		  " as positive data, please supply only one positive file.\n";
		exit 1;
	}
	$ret{pos_file} = ($ARGV[0]) ? $ARGV[0] : '';
	return %ret;
}

# MAIN BLOCK
use kalgn;
use pivalgn;

my %args = load_args;

# make sure no temp files are left when killed
use sigtrap 'handler' => sub {
	print "\r\033KRemoving temporary files...\n" if ($_[0] eq 'INT');
	unlink ".pssm.gxw" if ( -f ".pssm.gxw");
	unlink $args{neg_file} if $args{delete_neg};
	exit 1;
}, 'normal-signals';

print (<DATA>) and
  exit 0 if ($args{help});
print ("motif.pl version $VERSION\n2013 (C) Ron Mordechai\n") and
  exit 0 if ($args{version});

# check variable arguments
unless ( $args{lengthening} eq "pivalgn" || $args{lengthening} eq "4algn" ) {
	print "Invalid lengthening algorithm. Options are `4algn' or `pivalgn'.\n";
	exit 1;
}

if ($args{pos_file} eq '') {
	print "You must supply positive data\n";
	exit 1;
}

if ($args{neg_file} eq '') {
	$args{neg_file} = generate_neg_file $args{pos_file};
	$args{delete_neg} = 1;
}

open my $pos_file, "<", $args{pos_file} or
  die "Couldn't open $args{pos_file}: '$?'";
open my $neg_file, "<", $args{neg_file} or
  die "Couldn't open $args{neg_file}: '$?'";

$args{out_dir} =~ s/\/$//;

mkpath $args{out_dir} unless -d $args{out_dir};

# count lines and k-mers in each file
my $lines; ++$lines while <$pos_file>;
my %pos_tbl = count_k_mers $args{k}, $pos_file;
my %neg_tbl = count_k_mers $args{k}, $neg_file;

close $pos_file;
close $neg_file;

# calculate hypergeometric p-values for each k-mer
my %p_value_tbl =
  compute_hypergeometric_p_value \%pos_tbl, \%neg_tbl, $lines;

# if the user asks for one, generate histogram and exit
if ( $args{histogram} ) {
	my @hist = gen_histogram \%p_value_tbl;
	gen_histogram_graph \@hist, $args{out_dir};
	exit 0
}

# remove any k-mers under a calculated or given threshold
my %cleaned_tbl;
if ( $args{no_thres} ) {
	%cleaned_tbl = no_thres_clean_up \%p_value_tbl, $args{no_thres};
} elsif ( $args{threshold} ) {
	%cleaned_tbl = clean_up \%p_value_tbl, $args{threshold};
} else {
	my @hist = gen_histogram \%p_value_tbl;
	my $thres = get_coef_pval_threshold \@hist;
	%cleaned_tbl = clean_up \%p_value_tbl, $thres;
}

unlink $args{neg_file} if $args{delete_neg};

if ( !keys %cleaned_tbl ) {
	print "No motif found!\n";
	print "This is usually a  result of poorly permuted negative data.\n";
	print "Try supplying your own negative data set. If you supplied a\n";
	print "negative set, recheck it or let the script generate one for\n";
	print "you (run without `-n').\n";
	exit 1;
}

# cluster k-mers, then align them and lengthen the motif
my $align_k_mers =
  ( $args{lengthening} eq "4algn"   ) ? \&ka_align_k_mers :
  ( $args{lengthening} eq "pivalgn" ) ? \&pa_align_k_mers : undef;
my @clusters =
  ( $args{no_cluster} ) ? ( \%cleaned_tbl ) : cluster \%cleaned_tbl;

foreach my $cluster (@clusters) {
	my %aligned_k_mers = &$align_k_mers($cluster, $args{pos_file});

	my @PSSM = make_pssm_array \%aligned_k_mers;
	my $cPSSM = clean_pssm_array(\@PSSM);

	if ($args{validate}) {
		my $pos = validate_pos_dat $cPSSM, $args{pos_file};
		system("$0", "$pos", "-N", "-l", "$args{lengthening}",
			   "-o", "$args{out_dir}");
		unlink $pos;
	} else {
		my $label = pssm_to_label $cPSSM;
		make_pssm_logo $cPSSM, $label, $args{out_dir};
	}
}
__DATA__
Usage: motif.pl [options] <positive data>
Discriminatively search for a motif of length k within the provided file.

  -h, --help                   Display this help and exit.
  -H, --histogram              Generate a -log10(p-value) histogram instead of
                               a PSSM.
  -k, --start-count=K          Starts with K-mers of size K. K=6 by default.
  -l, --lengthening=ALGORITHM  Set lengthening algorithm, options are 4algn
                               and pivalgn. Uses 4algn by default.
  -L, --validate               Validate positive data against the PSSM.
                               Slow and broken as of now.
  -n, --neg-file=FILE          Set a negative data file instead of generating
                               one. A completely random negative file may
                               increase accuracy.
  -N, --no-cluster             Disable clustering. Outputs only one motif but
                               may increase accuracy.
  -o, --out-dir=DIR            Set an output directory, defaults to `./'.
  -t, --force-threshold=T      Force a p-value threshold of 10^-T. All k-mers
                               with a p-value higher than 10^-T are  removed.
  -T, --no-threshold=N         Instead of keeping all k-mers with a p-value
                               under a certain threshold, keep the top N
                               k-mers with the lowest p-values.
  -V, --version                Display version and exit.

2013 (C) Ron Mordechai
Released under the Perl Artistic License
