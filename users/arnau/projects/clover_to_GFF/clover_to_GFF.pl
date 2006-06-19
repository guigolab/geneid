#!/usr/local/bin/perl -w

use strict;

my $_algorithm = "clover";
my $_debug     = 0;
my $filename   = shift;

# Parsing clover output file

my %motifs_features_per_sequence;
my $features = [];
my $seqId;
my $start_parsing = 0;

open FILE, $filename or die "can't open input file, $filename!\n";
while (my $line = <FILE>) {
    if ($line =~ /^>([^\s]+)\s.+/) {
	# sequence header, which means, new set of features
	$seqId = $1;
	$start_parsing++;
	
	if (@$features > 0) {
	    # Store previous putative binding sites dataset
	    
	    $motifs_features_per_sequence{$seqId} = $features;
	    
	    # Reinitialisation
	    $features = [];
	}
    }
    elsif ($start_parsing) {
	# get rid of the spaces at the beginning of the line
	$line =~ s/^\s+//;
	
	if ($_debug) {
	    print STDERR "line after space removal, $line\n";
	}
	
	if ($line =~ /^(.+)\s+(\w+)\W+(\w+)\s+(\S)\s+([^\s]+)\s+(.+)/) {
	    my $motif_id  = $1;
	    my $start     = $2;
	    my $end       = $3;
	    my $strand    = $4;
	    my $motif_seq = $5;
	    my $score     = $6;
	    
	    $motif_id =~ s/\s+$//;
	    
	    if ($_debug) {
		print STDERR "motif_id, $motif_id, start, $start, end, $end, strand, $strand, sequence, $motif_seq, score, $score\n";
	    }
	    
	    if ($motif_id =~ /MA/) {
		# it is a jaspar motif
		# ...
	    }
	    
	    # Store the information
	    
	    my %feature = (
			   seqId     => $seqId,
			   algorithm => $_algorithm,
			   motif_id  => $motif_id,
			   start     => $start,
			   end       => $end,
			   score     => $score,
			   strand    => $strand,
			   frame     => ".",
			   motif_seq => $motif_seq,
			   );
	    push (@$features, \%feature);
	}
    }
}
close FILE;

# Printing out in GFF

foreach my $seqId (keys (%motifs_features_per_sequence)) {
    my $features = $motifs_features_per_sequence{$seqId};
    foreach my $feature_href (@$features) {
	my $algorithm = $feature_href->{algorithm} || die "algorithm not defined!\n";
	my $motif_id  = $feature_href->{motif_id}  || die "motif identifier not defined!\n";
	my $start     = $feature_href->{start}     || die "start not defined!\n";
	my $end       = $feature_href->{end}       || die "end not defined!\n";
	my $score     = $feature_href->{score}     || die "score not defined!\n";
	my $strand    = $feature_href->{strand}    || die "strand not defined!\n";
	my $frame     = $feature_href->{frame}     || die "frame not defined!\n";
	my $motif_seq = $feature_href->{motif_seq} || die "motif sequence not defined!\n";
	
	print "$seqId\t$algorithm\t$motif_id\t$start\t$end\t$score\t$strand\t$frame\t# $motif_seq\n";
    }
}

