#!/usr/local/bin/perl -w

# e.g. perl GFF2_to_GFF3.pl features.gff2 > features.gff3
# Can deal with GeneID, MatScan, Meta-alignment GFF output

###############################################
#
# Still to do - support partial CDSs if it the last reported exon !!
#
###############################################

use strict;
use Getopt::Std;

sub help {
my $usage = "
Description: Converts a set of features in GFF2 format into GFF3 format
Usage:

    GFF2_to_GFF3.pl [-h] -f {input file} -p {Gene Identifier Prefix} {Gene Identifier Index Start}
	-h help
	-f Input file with the set of features in GFF3 format
        -p Gene Identifier Prefix (Optional, no prefix by default)
        -i Gene Identifier Index start, Default is 1

Examples using some combinations:
    GFF2_to_GFF3.pl -f AC005155.geneid.gff2 > AC005155.geneid.gff3
    GFF2_to_GFF3.pl -f AC005155.geneid.gff2 -p Dm -i 26000000 > AC005155.geneid.gff3
	
";

print STDERR "$usage\n";
exit 1;

}

BEGIN {
	
    # Global variables definition
    use vars qw /$_debug %opts $in_file $id_prefix $id_index_start/;
    
    # these are switches taking an argument (a value)
    my $switches = 'hf:p:i:';
    
    # Get the switches
    getopts($switches, \%opts);
    
    # If the user does not write anything, skip to help
    if (defined($opts{h}) || !defined($opts{f})){
	print STDERR "Specify an input file name!\n";
	print STDERR help;
	exit 0;
    }
    
}

# Default values
$id_prefix      = "";
$id_index_start = 0;

defined $opts{f} and $in_file        = $opts{f};
defined $opts{p} and $id_prefix      = $opts{p};
defined $opts{i} and $id_index_start = $opts{i};

$_debug = 0;

# Default => no filtering
my $_score_filtering = -1;
# meta-alignment filtering
$_score_filtering = 7;

my @gff_features = ("## gff-version 3");

my $seqId;
my $algorithm;
my $geneId;
my $mRNAId;
my $geneStrand;
my $geneScore;
my $geneStart;
my $geneEnd;

my $parsed_first_exon = 0;
my $cds_features      = [];
my $has_genes         = 0;

open FILE, "<$in_file" or die "can't open input file, $in_file!\n";
while (<FILE>) {
    my $line = $_;
    
    # Meta-alignment Global score
    
    if ($line =~ /^# Maximum similarity/) {

	if ($_debug) {
	    print STDERR "Sequence score similarity information\n";
	}
	
	chomp $line;
	
	push (@gff_features, $line);
    }
    
    # Sequence information
    
    if ($line =~ /^# Sequence ([^\s]+)\s\D+(\d+) bps/) {
	
	if ($has_genes) {
	    
	    # The last gene of the previous sequence if the GFF file contains more than one sequence
	    
	    # Need to process an extra information, re. the last exon if it is partial gene or not !!!
	    
	    # ...
	    
	    my $gene_feature = "$seqId\t$algorithm\tgene\t$geneStart\t$geneEnd\t$geneScore\t$geneStrand\t.\tID=$geneId";
	    my $mRNA_feature = "$seqId\t$algorithm\tmRNA\t$geneStart\t$geneEnd\t$geneScore\t$geneStrand\t.\tID=$mRNAId;Parent=$geneId";
	    push (@gff_features, $gene_feature, $mRNA_feature, @$cds_features);
	}
	
	if ($_debug) {
	    print STDERR "Sequence information\n";
	}
	
	$seqId  = $1;
	my $length = $2;
	my $gff3_seq_line = "# Sequence-region $seqId 1 $length";
	push (@gff_features, $gff3_seq_line);
    }
    
    # Gene information
    
    if ($line =~ /^# Gene (\d+) \(([^\)]+)\)[^S]+Score\s*=\s*(.+)/) {
	
	if ($_debug) {
	    print STDERR "Gene information\n";
	}
	
	$has_genes++;

	if ($parsed_first_exon) {
	    $parsed_first_exon = 0;
	    
	    # Need to process an extra information, re. the last exon if it is partial gene or not !!!
	    
	    # ...
	    
	    my $gene_feature = "$seqId\t$algorithm\tgene\t$geneStart\t$geneEnd\t$geneScore\t$geneStrand\t.\tID=$geneId";
	    my $mRNA_feature = "$seqId\t$algorithm\tmRNA\t$geneStart\t$geneEnd\t$geneScore\t$geneStrand\t.\tID=$mRNAId;Parent=$geneId";
	    push (@gff_features, $gene_feature, $mRNA_feature, @$cds_features);
	    
	    $cds_features = [];
	}
	
	my $gene_index = $1;
	$geneStrand    = $2;
	$geneScore     = $3;
	
	# Set up the gene and mRNA identifiers
	
	my $index = $id_index_start + $gene_index;
	$geneId = $id_prefix . "_geneid_" . $index;
	$mRNAId = $id_prefix . "_geneid_mRNA_" . $index;

	# $geneId    = $seqId . "_gene_" . $gene_index;
	# $mRNAId    = $seqId . "_mRNA_" . $gene_index;
	
	if ($geneStrand =~ /forward/i) {
	    $geneStrand  = "+";
	}
	else { $geneStrand = "-"; }
    }
    
    if ($line =~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t*([^\t]*)/) {
	
	if ($_debug) {
	    print STDERR "Feature information\n";
	}
	
	$seqId     = $1;
	$algorithm = $2;
	my $featureType = $3;
	my $start  = $4;
	my $end    = $5;
	my $score  = $6;
	my $strand = $7;
	my $phase  = $8;
	chomp $phase;
	my $attributes = $9;
	
	# Clean the score
	$score =~ s/\s//g;
	
	# Filtering meta-alignment hits
	if ($algorithm =~ /meta_/) {
	    if ($score < $_score_filtering) {
		next;
	    }
	}
	
	# Feature type mapping if genes
	if ($featureType =~ /Internal|First|Terminal|Single/) {
	    
	    $attributes = undef;
	    
	    if ($_debug) {
		print STDERR "it is an exon feature\n";
	    }
	    
	    if (! $parsed_first_exon) {
		$parsed_first_exon = 1;
		$geneStart = $start;
		
		# if featureType is equal to Internal then it is a partial CDS
		if ($featureType eq "Internal") {
		    if ($strand eq "+") {
			$attributes .= "5_prime_partial=true";
		    }
		    else {
			$attributes .= "3_prime_partial=true";
		    }
		}
		
	    }
	    
	    $featureType = "CDS";
	    if (! defined $attributes) {
		$attributes  = "Parent=$mRNAId";
	    }
	    else {
		$attributes  = "Parent=$mRNAId;" . $attributes;
	    }
	    
	    my $cds_feature = "$seqId\t$algorithm\t$featureType\t$start\t$end\t$score\t$strand\t$phase\t$attributes";
	    push (@$cds_features, $cds_feature);
	    
	    $geneEnd = $end;
	}
	elsif ($algorithm =~ /MatScan|meta/) {
	    
	    if ($_debug) {
		print STDERR "it is an binding_site feature\n";
	    }
	    
	    my $featureId   = $featureType;
	    $featureType    = "binding_site";
	    
	    if ($algorithm eq "MatScan") {
		
		# Store also the binding_site sequence
		
		if ($attributes =~ /^# (\w+)/) {
		    my $sequenceAttribute = "Seq=" . $1;
		    $attributes  = "ID=$featureId;$sequenceAttribute";
		}
		else {
		    $attributes  = "ID=$featureId";
		}
	    }
	    else {
		$attributes     = "ID=$featureId";
	    }
	    
	    # Figure out the Dbxref attribute
	    
	    if ($featureId =~ /^\w\$/) {
		# It is Transfac identifier
		my $Dbxref = "Dbxref=Transfac:$featureId";
		$attributes  .= ";" . $Dbxref;
	    }
	    elsif ($featureId =~ /^MEME/) {
		# MEME motif
	    }
	    else {
		# Must be a Jaspar motif !!!
		my $Dbxref = "Dbxref=Jaspar:$featureId";
		$attributes  .= ";" . $Dbxref;
	    }
	    
	    my $gff_feature = "$seqId\t$algorithm\t$featureType\t$start\t$end\t$score\t$strand\t$phase\t$attributes";
	    push (@gff_features, $gff_feature);
	}
	else {
	    
	    if ($_debug) {
		print STDERR "it is a feature of unprocessed type\n";
	    }
	    
	    chomp $line;
	    push (@gff_features, $line);
	}
	
    }    
}
close FILE;

if ($has_genes) {
    # The last gene
    
    # Need to process an extra information, re. the last exon if it is partial gene or not !!!
    
    # ...
    
    my $gene_feature = "$seqId\t$algorithm\tgene\t$geneStart\t$geneEnd\t$geneScore\t$geneStrand\t.\tID=$geneId";
    my $mRNA_feature = "$seqId\t$algorithm\tmRNA\t$geneStart\t$geneEnd\t$geneScore\t$geneStrand\t.\tID=$mRNAId;Parent=$geneId";
    push (@gff_features, $gene_feature, $mRNA_feature, @$cds_features);
}

print join ("\n", @gff_features);
print "\n";
