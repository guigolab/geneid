#!/usr/local/bin/perl -w

# Generate a score matrix, given various input data
# so far, you can give:
# * meta-alignment format output
# * GFF format output

# e.g. generateScoreMatrice.pl meta-alignment sota < data/runMultiMetaAlignment.1000.out > scores.matrice

# you can get different types of formatting for your matrices, e.g. for fneighbor, it's got to be space delimited, in other cases it's got to be tab delimited.

use strict;
use Data::Dumper;

my $_debug      = 0;
my $_delimitor  = "\t";

# inbHierarchicalCluster not working with dashes!!
my $_equal_mark    = "-";
$_equal_mark    = "";

my $_plus_infinite_mark  = "1000000";
my $_minus_infinite_mark = "-1000000";
$_plus_infinite_mark  = "-";
$_minus_infinite_mark = "-";

$_plus_infinite_mark  = "";
$_minus_infinite_mark = "";

# SOTA is compliant with valencia and uam

my $inputformat  = shift || "meta-alignment";
my $outputformat = shift || "SOTA";

my $scores_per_sequences;
my $sequence_identifiers;

# Parsing input file

if ($inputformat eq "meta-alignment") {

    if ($_debug) {
	print STDERR "parsing meta-alignment input data...\n";
    }
    
    ($sequence_identifiers, $scores_per_sequences) = parseMeta ();
}
else {
    print STDERR "unknown input format, $inputformat\n";
    exit 0;
}

# Generate the matrix, as an array of arrays

my @matrix = ();

for (my $i=0; $i < @$sequence_identifiers; $i++) {
    my @scores_per_row = ();
    my $seq_identifier_1 = $sequence_identifiers->[$i];
    push (@scores_per_row, $seq_identifier_1);
    
    for (my $j=0; $j < @$sequence_identifiers; $j++) {

	if ($i == $j) {
	    push (@scores_per_row, $_plus_infinite_mark);
	}

	if ($i < $j) {
	    my $seq_identifier_2 = $sequence_identifiers->[$j];
	    my $id    = $seq_identifier_1 . "_" . $seq_identifier_2;
	    my $score;
	    if (defined $scores_per_sequences->{$id}) {
		$score = $scores_per_sequences->{$id};
	    }
	    else {
		# it's the other way around !!
		$id    = $seq_identifier_2 . "_" . $seq_identifier_1;
		$score = $scores_per_sequences->{$id} || $_minus_infinite_mark;
	    }
	    
	    push (@scores_per_row, $score);
	}
	if ($i > $j) {
	    my $seq_identifier_2 = $sequence_identifiers->[$j];
	    my $id    = $seq_identifier_2 . "_" . $seq_identifier_1;
	    my $score;
	    if (defined $scores_per_sequences->{$id}) {
		$score = $scores_per_sequences->{$id};
	    }
	    else {
		# it's the other way around !!
		$id    = $seq_identifier_1 . "_" . $seq_identifier_2;
		$score = $scores_per_sequences->{$id} || $_minus_infinite_mark;
	    }
	    push (@scores_per_row, $score);
	}
    }
    push (@matrix, \@scores_per_row);
}

# Format the matrix for output

# SOTA compliant formatting

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
    
    my $scores_per_sequences = {};
    my $sequence_identifiers = {};

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

		# Add the sequence identifier in the seqids hash if not in there yet.
		if ((not defined $sequence_identifiers->{$map1}) || (not $sequence_identifiers->{$map1})) {
		    $sequence_identifiers->{$map1} = 1;
		}
		
		# parse map2
		
		# one line further...
		$line = <STDIN>;
		if ($line =~ /MAP2 ([^\s]+) .+/) {
		    # parse map2 identifier
		    $map2 = $1;

		    # Add the sequence identifier in the seqids hash if not in there yet.
		    if ((not defined $sequence_identifiers->{$map2}) || (not $sequence_identifiers->{$map2})) {
			$sequence_identifiers->{$map2} = 1;
		    }

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
	print STDERR "parsed " . keys (%$sequence_identifiers) . " sequences.\n";
	print STDERR "parsed " . keys (%$scores_per_sequences) . " scores.\n";
    }
    
    my @seqids = keys (%$sequence_identifiers);
    return (\@seqids, $scores_per_sequences);
}
