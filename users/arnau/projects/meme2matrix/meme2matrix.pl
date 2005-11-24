#!/usr/local/bin/perl -w

use strict;

my $_debug = 0;

my $in_file = shift;
# score matrix or probability matrix mode
my $mode    = shift;

# mode eq "score" || "probability"

if (not ($mode =~ /^score$/i || $mode =~ /^probability$/i)) {
    print STDERR "you have to select a mode, 'score' or 'probability'\n";
}

my $matrices = [];
my $motif_index = 1;

open FILE, "<$in_file" or die "can't open file, $in_file!\n";
while (<FILE>) {
    my $line = $_;
    if ($mode =~ /^score$/i) {
	if ($line =~ /position-specific scoring matrix/) {
	    my $index = 1;
	    # Jump the next line
	    $line = <FILE>;
	    # parse from this line upto "----" line
	    $line = <FILE>;
	    my $matrix = "MEME_Motif_" . $motif_index . " $line";
	    $motif_index++;
	    $line = <FILE>;
	    while (not ($line =~ /^-/)) {
		$matrix .= "$index   $line";
		$line = <FILE>;
                $index++;
	    }
	    # add //
	    $matrix .= "//";
	    push (@$matrices, $matrix);
	}
    }
    elsif ($mode =~ /^probability$/i) {
	if ($line =~ /position-specific probability matrix/) {
	    my $index = 1;
	    # Jump the next line
	    $line = <FILE>;
	    # parse from this line upto "----" line
	    $line = <FILE>;
	    my $matrix = "MEME_Motif_" . $motif_index . " $line";
	    $motif_index++;
	    $line = <FILE>;
	    while (not ($line =~ /^-/)) {
		$matrix .= "$index   $line";
		$line = <FILE>;
		$index++;
	    }
	    # add //
	    $matrix .= "//";
	    push (@$matrices, $matrix);
	}
    }
}
close FILE;

if ($_debug) {
    print STDERR "matrices:\n";
}
print join ("\n", @$matrices);
