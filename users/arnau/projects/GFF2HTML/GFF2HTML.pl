#!/usr/local/bin/perl -w

# Replace 'motif' by 'feature' or $feature_key !

use strict;
use POSIX qw(ceil floor);

my $_debug = 0;
my @input_files = @ARGV;

if ($_debug) {
    print STDERR "files: @input_files.\n";
}

# could be 'Lowest p-value' if MEME
my $score_type = "Score";

# Get this from the command line in case of meta-alignment or MatScan !!!

my $seq_start = 1;
my $seq_end = 500;
my $seq_length = 500;

# Default
my $gap   = 25;
my $gap_width = 50;
my $gap_width_adjusted = $gap_width - 2;

my %sequences;

# parse GFF3 input file

foreach my $input_file (@input_files) {

    my $index = 1;

    if ($_debug) {
	print STDERR "parsing $input_file...\n";
    }
    
    my %features;
    my @seq_starts = ();
    my %feature_index_by_seq_starts;
    my $seq_id;
    my $seq_score = "&nbsp";
    
    open FILE, "<$input_file" or die "can't open file, $input_file!\n";
    while (<FILE>) {
	my $line = $_;
	
	# Parse the sequence length
	
	if ($line =~ /\# Sequence.region [^\s]+\s([^\s]+)\s(\d+)/) {
	    if ($_debug) {
		print STDERR "parsing sequence length...\n";
	    }
	    $seq_start   = $1;
	    $seq_end     = $2;
	    $seq_length = $seq_end - $seq_start + 1;
	}
	
	# Parse the meta-alignment sequence similarity score
	
	if ($line =~ /^\# Maximum similarity\: (.+)/) {
	    $seq_score = $1;
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
	    
	    if ($_debug) {
		print STDERR "start, $start\n";
		print STDERR "end, $end\n";
	    }
	    
	    if ($method =~ /geneid/i) {
		# keep only exon features
		if ($feature_type ne "CDS") {
		    
		    if ($_debug) {
			print STDERR "filtering $feature_type\n";
		    }
		    
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
	    if (! is_in ($seq_start, @seq_starts)) {
		push (@seq_starts, $seq_start);
	    }
	    if (defined $feature_index_by_seq_starts{$seq_start}) {
		my $feature_indexes_tmp =  $feature_index_by_seq_starts{$seq_start};
		push (@$feature_indexes_tmp, $index);
		$feature_index_by_seq_starts{$seq_start} = $feature_indexes_tmp;
	    }
	    else {
		$feature_index_by_seq_starts{$seq_start} = [$index];
	    }
	    
	    $index++;
	    
	}
	
    }
    close FILE;
    
    if (defined $seq_length) {
	if ($_debug) {
	    print STDERR "sequence length, $seq_length\n";
	}
    }
    else {
	print STDERR "no sequence length defined!\n";
	exit 1;
    }
    
    my $nb_features = keys (%features);
    if ($_debug) {
	print STDERR "nb features, $nb_features\n";
    }
    
    # Sort the sequence starts
    @seq_starts = sort {$a <=> $b} @seq_starts;
    
    my $sequence = {
	seq_id     => $seq_id,
	seq_length => $seq_length,
	features   => \%features,
	seq_starts => \@seq_starts,
	feature_index_by_seq_starts => \%feature_index_by_seq_starts,
	seq_score  => $seq_score,
    };
    
    $sequences{$seq_id} = $sequence;
    
}

$gap = floor ($seq_length / 20);

if ($_debug) {
    print STDERR "gap, $gap\n";
}

# HTML generation

my $width = 1032;
my $length_feature = 1;

print "<html>\n<body BGCOLOR='#D5F0FF'>\n";
print "<CENTER><BIG><B>Features block diagrams</B></BIG></CENTER><HR>\n";

print "<TABLE SUMMARY='feature diagrams' BORDER=1 ALIGN=CENTER>\n";
print "<TR><TH>Name<TH>Lowest<BR>$score_type<TH ALIGN=LEFT>&nbsp;&nbsp; Features\n";

my @seq_ids = keys ( %sequences);
foreach my $seq_id (@seq_ids) {
    
    my $sequence_href = $sequences{$seq_id};
    
    my $seq_score = $sequence_href->{seq_score} || "&nbsp";
    
    print "<TR>\n";
    
    print "<TD>$seq_id\n";
    print "<TD ALIGN=RIGHT NOWRAP>$seq_score\n";
    
    my $features_href   = $sequence_href->{features};
    my $feature_index_by_seq_starts_href = $sequence_href->{feature_index_by_seq_starts};
    my $seq_starts_aref = $sequence_href->{seq_starts};
    
    print "<TD>\n";
    
    my $before_feature_width;
    
    print "<TABLE SUMMARY='diagram $seq_id' WIDTH=$width BORDER=0 ALIGN=LEFT CELLSPACING=0 CELLPADDING=0><TR ALIGN=CENTER>\n";
    if (keys %$features_href > 0) {
	
	my $start_reference = 1;
	
	# Loop over the features but in the order of position onto the sequence from 5' to 3' !!!!!!!!!!!
	
	foreach my $seq_start (@$seq_starts_aref) {
	    my $feature_indexes_aref = $feature_index_by_seq_starts_href->{$seq_start};
	    
	    foreach my $feature_index (@$feature_indexes_aref) {
		
		my $feature = $features_href->{$feature_index};
		my $start   = $feature->{start};
		my $end     = $feature->{end};
		$before_feature_width = floor (($start - $start_reference) * $width / $seq_length);
		
		if ($_debug) {
		    print STDERR "feature index, $feature_index\n";
		    print STDERR "exon start, $start\n";
		    print STDERR "before_feature_width, $before_feature_width\n";
		}
		
		print "<TD WIDTH=$before_feature_width><HR SIZE=4 NOSHADE>\n";
		print "<TD CLASS='c0' WIDTH=$length_feature>+$feature_index\n";
		
		$start_reference = $start;
	    }
	    
	    my $remaining = $width - $before_feature_width - $length_feature;
	    
	    if ($_debug) {
		print STDERR "remaining, $remaining\n";
	    }
	    print "<TD WIDTH=$remaining><HR SIZE=4 NOSHADE>\n";
	}
    }
    else {
	print "<TD WIDTH=$width><HR SIZE=4 NOSHADE>\n";
    }

    print "</TABLE>\n";
    
   # score not here anymore !!!

}

# Scale here now !!
print "<TR><TH CLASS='blue' COLSPAN=2 ROWSPAN=2 ALIGN=LEFT>SCALE\n";
print "<TD><TABLE SUMMARY='scale' WIDTH=$width BORDER=0 ALIGN=LEFT CELLSPACING=0 CELLPADDING=0><TR ALIGN=CENTER>\n";

# The scale

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


####

print "</TABLE>\n";

print "</body>\n</html>\n";


sub is_in {
    my ($element, @elements) = @_;
    
    foreach my $element_test (@elements) {
	if ($element_test == $element) {
	    return 1;
	}
    }
    
    return 0;
}
