#!/usr/local/bin/perl -w

use strict;

my $in_file = shift;
# score matrix or probability matrix mode
my $mode    = shift;

# mode eq "score" || "probability"

if (not ($mode =~ /^score$/i || $mode =~ /^probability$/i)) {
    print STDERR "you have to select a mode, 'score' or 'probability'\n";
}

my $matrices = [];

open FILE, "<$in_file" or die "can't open file, $in_file!\n";
while (<FILE>) {
    my $line = $_;
    if ($mode =~ /^score$/i) {
	if ($line =~ /position-specific scoring matrix/) {
	    # Jump the next line
	    $line = <FILE>;
	    # parse from this line upto "----" line
	    my $matrix = "";
	    $line = <FILE>;
	    while (not ($line =~ /^-/)) {
		$matrix .= $line;
		$line = <FILE>;
	    }
	    push (@$matrices, $matrix);
	}
    }
    elsif ($mode =~ /^probability$/i) {
	if ($line =~ /position-specific probability matrix/) {
	    # Jump the next line
	    $line = <FILE>;
	    # parse from this line upto "----" line
	    my $matrix = "";
	    $line = <FILE>;
	    while (not ($line =~ /^-/)) {
		$matrix .= $line;
		$line = <FILE>;
	    }
	    push (@$matrices, $matrix);
	}
    }
}
close FILE;

print STDERR "matrices:\n";
print join ("\n", @$matrices);
