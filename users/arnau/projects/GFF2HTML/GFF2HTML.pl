#!/usr/local/bin/perl -w

# Replace 'motif' by 'feature' or $feature_key !

use strict;
use POSIX qw(ceil floor);

my $_debug = 0;
my @input_files = @ARGV;

print STDERR "files:@input_files.\n";

# could be 'Lowest p-value' if MEME
my $score_type = "Score";
my $seq_start;
my $seq_end;
my $seq_length;
# Default
my $gap   = 25;
my $gap_width = 50;
my $gap_width_adjusted = $gap_width - 2;

my $index = 1;

my %sequences;

# parse GFF3 input file

foreach my $input_file (@input_files) {
    
    print STDERR "parsing $input_file...\n";
    
    my %features;
    my $seq_id;
    
    open FILE, "<$input_file" or die "can't open file, $input_file!\n";
    while (<FILE>) {
	my $line = $_;
	
	# parse the sequence length
	
	if ($line =~ /\# Sequence.region [^\s]+\s([^\s]+)\s(\d+)/) {
	    print STDERR "parsing sequence length...\n";
	    $seq_start   = $1;
	    $seq_end     = $2;
	    $seq_length = $seq_end - $seq_start + 1;
	}
	
	if ($line =~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t*([^\t]*)/) {
	    
	    if ($_debug) {
		print STDERR "Feature information\n";
	    }
	    
	    $seq_id    = $1;
	    my $method = $2;
	    my $feature_type = $3;
	    my $start  = $4;
	    my $end    = $5;
	    my $score  = $6;
	    my $feature_id;
	    
	    if ($method =~ /geneid/i) {
		# keep only exon features
		if ($feature_type ne "exon") {
		    next;
		}
	    }
	    
	    my $feature = {
		feature_type => $feature_type,
		start => $start,
		end   => $end,
		score => $score,
	    };
	    
	    $features{$index} = $feature;
	    
	    $index++;
	    
	}
	
    }
    close FILE;
    
    if (defined $seq_length) {
	print STDERR "sequence length, $seq_length\n";
    }
    else {
	print STDERR "no sequence length defined!\n";
	exit 1;
    }
    
    my $nb_features = keys (%features);
    print STDERR "nb features, $nb_features\n";

    my $sequence = {
	seq_id => $seq_id,
	seq_length => $seq_length,
	features => \%features,
	};
    
    $sequences{$seq_id} = $sequence;
    
}

# $gap = 25 * $seq_length / 500;
$gap = floor ($seq_length / 20);

print STDERR "gap, $gap\n";

# HTML generation

print "<html>\n<body BGCOLOR='#D5F0FF'>\n";
print "<CENTER><BIG><B>Features block diagrams</B></BIG></CENTER><HR>\n";

print "<TABLE SUMMARY='feature diagrams' BORDER=1 ALIGN=CENTER>\n";
print "<TR><TH>Name<TH>Lowest<BR>$score_type<TH ALIGN=LEFT>&nbsp;&nbsp; Features\n";

my @seq_ids = keys ( %sequences);
foreach my $seq_id (@seq_ids) {

    my $sequence_href = $sequences{$seq_id};
    
    print "<TR>\n";
    
    my $score = $sequence_href->{score} || "&nbsp;";
    
    print "<TD>$seq_id\n";
    print "<TD ALIGN=RIGHT NOWRAP>$score\n";
    
    my $features_href  = $sequence_href->{features};

    print "<TD>\n";

    my $width = 1000 + 32;
    
    print "<TABLE SUMMARY='diagram $seq_id' WIDTH=$width BORDER=0 ALIGN=LEFT CELLSPACING=0 CELLPADDING=0><TR ALIGN=CENTER>\n";
    if (keys %$features_href > 0) {

	# in the meantime
	# print "<TD WIDTH=$width><HR SIZE=4 NOSHADE>\n";
	
	my @feature_indexes = keys (%$features_href);
	foreach my $feature_index (@feature_indexes) {
	    my $feature = $features_href->{$feature_index};
	    my $start   = $feature->{start};
	    my $end     = $feature->{end};
	    my $before_feature_width = $seq_start * 2;
	    print "<TD WIDTH=$before_feature_width><HR SIZE=4 NOSHADE>\n";
	    print "<TD CLASS='c0' WIDTH=30>+$feature_index\n";
	    
	    # temporary
	    last;
	}
    }
    else {
	print "<TD WIDTH=$width><HR SIZE=4 NOSHADE>\n";
    }

    print "</TABLE>\n";
    print "<TR><TH CLASS='blue' COLSPAN=2 ROWSPAN=2 ALIGN=LEFT>SCALE\n";
    print "<TD><TABLE SUMMARY='scale' WIDTH=$width BORDER=0 ALIGN=LEFT CELLSPACING=0 CELLPADDING=0><TR ALIGN=CENTER>\n";

    print "<TD CLASS='blue' WIDTH=$gap_width_adjusted ALIGN=LEFT>|</TD>\n";
    for (my $i =1; $i < 20; $i++) {
	print "<TD CLASS='blue' WIDTH=$gap_width ALIGN=LEFT>|</TD>\n";
    }
    
    print "<TR>\n";
    
    print "<TD CLASS='blue' WIDTH=$gap_width_adjusted ALIGN=LEFT>$seq_start</TD>\n";
    
    for (my $i =1; $i < 20; $i++) {
	my $coordinate = $seq_start + ($i * $gap);
	print "<TD CLASS='blue' WIDTH=$gap_width ALIGN=LEFT>$coordinate</TD>\n";
    }
    print "</TABLE>\n";

}

print "</TABLE>\n";

print "</body>\n</html>\n";
