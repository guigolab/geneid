#!/usr/local/bin/perl -w

# Generate a score matrix, given various input data
# so far, you can give:
# * meta-alignment format output
# * GFF format output

# e.g. generateScoreMatrice.pl meta-alignment sota < data/runMultiMetaAlignment.1000.out > scores.matrice

# you can get different types of formatting for your matrices, e.g. for fneighbor, it's got to be space delimited, in other cases it's got to be tab delimited.

use strict;

my $_debug      = 0;
my $_delimitor  = "\t";
my $_equal_mark = "-";

my $inputformat  = shift || "meta-alignment";
my $outputformat = shift || "SOTA";

my $scores_per_sequences;
my $sequence_identifiers;

# Parsing input file

if ($inputformat eq "meta-alignment") {
    ($sequence_identifiers, $scores_per_sequences) = parseMeta ();
}
else {
    print STDERR "unknown input format, $inputformat\n";
    exit 0;
}

# Format the matrix, as an array of arrays

my @matrix       = ();

for (my $i=0; $i < @$sequence_identifiers; $i++) {
    my @scores_per_row = ();
    my $seq_identifier_1 = $sequence_identifiers->[$i];
    push (@scores_per_row, $seq_identifier_1);
    
    for (my $j=0; $j < @$sequence_identifiers; $j++) {

	if ($i == $j) {
	    push (@scores_per_row, $_equal_mark);
	}

	if ($i < $j) {
	    my $seq_identifier_2 = $sequence_identifiers->[$j];
	    my $id    = $seq_identifier_1 . "_" . $seq_identifier_2;
	    my $score = $scores_per_sequences->{$id} || $_equal_mark;
	    
	    push (@scores_per_row, $score);
	}
	if ($i > $j) {
	    my $seq_identifier_2 = $sequence_identifiers->[$j];
	    my $id    = $seq_identifier_2 . "_" . $seq_identifier_1;
	    my $score = $scores_per_sequences->{$id} || $_equal_mark;
	    
	    push (@scores_per_row, $score);
	}
    }
    push (@matrix, \@scores_per_row);
}

# Format the matrix for output

# Ximo Dopazo compliant formatting

if ($outputformat eq "SOTA") {
    
    print "#$_delimitor" . join ("$_delimitor", @$sequence_identifiers) . "\n";
    
    foreach my $row (@matrix) {
	my @row = @$row;
	print join ("$_delimitor", @$row) . "\n";
    }
}
else {
    print STDERR "unknown output format, $outputformat\n";
    exit 0;
}

#################
#
# The End...
#
#################

sub parseMeta {
    
    my $scores_per_sequences;
    my $sequence_identifiers = [];

    while (my $line = <STDIN>) {
	my $map1;
	my $map2;
	my $score;
	
	if ($line =~ /MAP1 ([^\s]+) .+/) {
	    # parse map1 identifier
	    $map1 = $1;
	    if (not defined $map1) {
		print STDERR "Error, can not parse map1!!\n";
		print STDERR "line, $line\n";
		exit 0;
	    }
	    else {

		# Add the sequence identifier in the seqids array if not in there yet.
		if ((not defined $sequence_identifiers->[0]) || ($sequence_identifiers->[@$sequence_identifiers-1] ne $map1)) {
		    push (@$sequence_identifiers, $map1);
		}
		
		# parse map2
		
		# one line further...
		$line = <STDIN>;
		if ($line =~ /MAP2 ([^\s]+) .+/) {
		    # parse map2 identifier
		    $map2 = $1;
		}
		else {
		    print STDERR "Error, can not parse map2!!\n";
		    print STDERR "line, $line\n";
		    exit 0;
		}
		
		# parse score
		
		# two lines further...
		<STDIN>;
		$line = <STDIN>;
		if ($line =~ /Maximum similarity: (.+)/) {
		    # parse map2 identifier
		    $score = $1;
		    chomp $score;
		}
		else {
		    print STDERR "Error, can not parse the similarity score!!\n";
		    print STDERR "line, $line\n";
		    exit 0;
		}
		
		if ($_debug) {
		    print STDERR "got $map1, $map2, $score.\n";
		}

		my $id = $map1 . "_" . $map2;
		$scores_per_sequences->{$id} = $score;
	    }
	}
    }

    if ($_debug) {
	print STDERR "parsed " . @$sequence_identifiers . " sequences.\n";
	print STDERR "parsed " . keys (%$scores_per_sequences) . " scores.\n";
    }
    
    return ($sequence_identifiers, $scores_per_sequences);
}
