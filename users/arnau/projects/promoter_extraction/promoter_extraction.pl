#!/usr/local/bin/perl -w

# NB: Ensembl start of gene definition: TSS or if not known, start codon
# Do we want ot include it in the output sequence ?

# ..........................

# Make sure we don't go over the chromosome start/end !!!
# Perldoc to make

# ..........................

# Issue warnings about suspicious programming.
use warnings 'all';

# Must declare and initialize all variables
use strict;

# be prepare for command-line options/arguments
use Getopt::Std;

use Data::Dumper;
use Benchmark;

# Bioperl
use Bio::SeqIO;
use Bio::PrimarySeq;

sub help {
return <<"END_HELP";
Description: Extract Promoter regions for a set of given genes
Usage:

    promoter_extraction.pl [-h] -f {gene ids} -s {Species} -r {Genome Release} -u {upstream} -d {downstream} -i {Flag intergenic regions only} -g {Flag report attached features in GFF format}
	-h help
	-f gene identifiers input file
	-s gender species, e.g. homo_sapiens
	-r database release, e.g 29_35b
	-u length of the upstream region   (Default 2000)
	-d length of the downstream region (Default 500)
	-i intergenic sequence only
	-g report overlapping gene features in GFF format

Examples using some combinations:
	promoter_extraction.pl -f geneIds.lst -s homo_sapiens -r 29_35b -u 2000 -d 500 -i -g
	promoter_extraction.pl -f refseqs.lst -s rattus_norvegicus -r 29_3f -u 2000 -d 500 -i -g

	promoter_extraction.pl -f geneIds.lst -s homo_sapiens -r 31_35d -u 2000 -d 500 -i -g
	promoter_extraction.pl -f geneIds.lst -s homo_sapiens -r 30_35c -u 2000 -d 500 -i -g

END_HELP

}

BEGIN {
	
    # Global variables definition
    use vars qw /$_debug %opts $report_features $intergenic_only $upstream_length $downstream_length $geneIds_file $dbhost $dbuser $dbensembl/;

    if (-f "/home/ug/gmaster/projects/promoter_extraction/config.pl") {
        require "/home/ug/gmaster/projects/promoter_extraction/config.pl";
    }
    else {
        print STDERR "There is a problem with the configuration, contact gmaster\@imim.es for help\n";
        exit 1;
    }
    
    # these are switches taking an argument (a value)
    my $switches = 'hf:s:r:u:d:ig';
    
    # Get the switches
    getopts($switches, \%opts);
    
    # If the user does not write nothing, skip to help
    if (defined($opts{h}) || !defined($opts{f}) || !defined($opts{s}) || !defined($opts{r})){
	print help;
	exit 0;
    }    
    
}

## Benchmarking code
my $time0 = new Benchmark;

my $species;
my $release;

$_debug            = 0;
$upstream_length   = 2000;
$downstream_length = 500;
$report_features   = 0;
$intergenic_only   = 0;

defined $opts{f} and $geneIds_file      = $opts{f};
defined $opts{s} and $species           = $opts{s};
defined $opts{r} and $release           = $opts{r};
defined $opts{u} and $upstream_length   = $opts{u};
defined $opts{d} and $downstream_length = $opts{d};
defined $opts{g} and $report_features++;
defined $opts{i} and $intergenic_only++;

# Ensembl database configuration

$dbhost      = "ensembldb.ensembl.org";
$dbuser      = "anonymous";
$dbensembl = $species . "_core_" . $release;

$release =~ /^(\d+)_\w+$/;
my $releaseNumber = $1;

my $ensembl_API_path = "/home/ug/gmaster/projects/promoter_extraction/lib/ensembl-$releaseNumber/modules";

if ($_debug) {
    print STDERR "ensembl_API_path, $ensembl_API_path.\n";
}

if (-d $ensembl_API_path) {
    # use lib "$ensembl_API_path";
    unshift (@INC, $ensembl_API_path);
}
else {
    print STDERR "release, $release, not supported, contact gmaster\@imim.es for help\n";
    exit 1;
}

use Bio::EnsEMBL::DBSQL::DBAdaptor;

if ($_debug) {
    print STDERR "connecting to Ensembl database...\n";
}

my $dbh = new Bio::EnsEMBL::DBSQL::DBAdaptor (
					      -host   => $dbhost,
					      -user   => $dbuser,
					      -dbname => $dbensembl
					      )
    or die "can't connect to Ensembl database, $dbensembl, contact gmaster\@imim.es for help!\n";

my $dbEntry_adaptor = $dbh->get_DBEntryAdaptor ();
my $gene_adaptor    = $dbh->get_GeneAdaptor ();
my $slice_adaptor   = $dbh->get_SliceAdaptor ();

# Parse input file

# Mapping : RefSeq DNA XM_346823 or AFFY_Rat230_2: 1374933_at, 1369793_a_at
# AFFY_RG_U34C: rc_AA859350_at

# to

# ensembl identifier: ENSRNOG00000007726

# my @geneIds = ("ENSG00000179403", "ENSG00000197185", "ENSG00000197785");
# my @geneIds = ("ENSRNOG00000007726");
my @geneIds = ();

open FILE, "<$geneIds_file" or die "can't open file, $geneIds_file!\n";
while (my $geneId = <FILE>) {
    chomp ($geneId);
    push (@geneIds, $geneId);
}
close FILE;

my $fasta_output_file = "upstream_regions.fa";
my $gff_output_file   = "upstream_regions.gff";

my $sout = Bio::SeqIO->new (
			    -file   => ">$fasta_output_file",
			    -format => 'fasta'
			    );

if ($report_features) {
    open GFF, ">$gff_output_file" or die "can't open file, $gff_output_file!\n";
}

if ($_debug) {
    print STDERR "processing the " . @geneIds . " gene identifiers...\n";
}

foreach my $geneId (@geneIds) {

    if ($_debug) {
	print STDERR "\nprocessing gene identifier, $geneId...\n";
    }
    
    # Get the stable identifier if not it
    
    # Mapping
    
    # RefSeq DNA XM_346823 or AFFY_Rat230_2: 1374933_at, 1369793_a_at
    # AFFY_RG_U34C: rc_AA859350_at
    # This is display_id or primary_id, different from dbid!
    # Only using dbid, 16078 ???
    
    # to
    
    # ensembl identifier: ENSRNOG00000007726

    my $genes = $gene_adaptor->fetch_all_by_external_name ($geneId);
    # my $genes = $gene_adaptor->fetch_all_by_external_name ("NM_023983");
    # my $genes = $gene_adaptor->fetch_all_by_external_name ("XM_346823");
    # Affy Id
    # $genes = $gene_adaptor->fetch_all_by_external_name ("1369793_a_at");

    # We don't know if the geneId is an external identifier or an Ensembl identifier

    # It is an external identifier
    
    if (@$genes > 1) {
	print STDERR "Found more than one Ensembl gene (" . @$genes . ") associated with external identifier, $geneId!!!\n";
	print STDERR "Will return the upsteam region of all Ensembl genes\n";
	
    }
    
    # It must be an Ensembl identifier
    
    if ((not defined $genes) || (@$genes == 0)) {
	my $gene = $gene_adaptor->fetch_by_stable_id ($geneId) or warn "can't instanciate gene object with identifier, $geneId!\n";
	if (not defined $gene) {
	    print STDERR "can't find any Ensembl gene for this given identifier, $geneId!\n";
	}
	else {
	    if ($_debug) {
		print STDERR "gene stable identifier found: " . $gene->stable_id . ", on strand, " . $gene->strand . "\n";
	    }
	    $genes = [ $gene ];
	}
    }

    foreach my $gene (@$genes) {

	my $gene_stable_id    = $gene->stable_id;
	my $coordinate_system = $gene->coord_system_name();
	my $chromosome        = $gene->seq_region_name();
	my $strand            = $gene->strand;
	my $start             = "";
	my $end               = "";
	my $features = [];
	
	my $tss               = $gene->start;
	if ($strand == -1) {
	    $tss = $gene->end;
	}

	if ($_debug) {
	    print STDERR "tss, $tss\n";
	}

	# Make sure we don't go over the chromosome start/end !!!
	
	if ($strand > 0) {
	    $start = $gene->start - $upstream_length;
	    $end   = $gene->start + $downstream_length;
	}
	else {
	    $end   = $gene->end + $upstream_length;
	    $start = $gene->end - $downstream_length;
	}
	
	if ($_debug) {
	    print STDERR "gene coordinates, " . $gene->start . ".." . $gene->end . " on strand $strand\n";
	}
	
	# Get rid of the base 0 if we don't want any downstream region, this way we won't get any overlapping feature associated with the given gene (its first exon for example)
	
	if (($downstream_length == 0) && ($strand == 1)) {
	    $end -= 1;
	}
	elsif (($downstream_length == 0) && ($strand == -1)) {
	    $start += 1;
	}
	
	if ($_debug) {
	    print STDERR "start and end of the slice region, $start, $end\n";
	}
	
	# Get the Slice object corresponding to the region of the genome we want to extract
	
	my $slice_region = $slice_adaptor->fetch_by_region ($coordinate_system,$chromosome,$start,$end,$strand);
	
	if ($_debug) {
	    print STDERR "length of the extracted region, " . $slice_region->length . "\n";
	}
	
	if ($intergenic_only) {
	    if (has_gene_upstream ($slice_region, $strand, $downstream_length, $tss)) {
		
		if ($_debug) {
		    print STDERR "has gene upstream\n";
		}
		
		my ($subslice_start, $subslice_end) = getIntergenicSequence ($slice_region, $strand, $downstream_length, $tss);
		
		if ($_debug) {
		    print STDERR "subslice start, end: $subslice_start, $subslice_end\n";
		}
		
		# Update the slice object so it includes only the intergenic sequence
		
		$slice_region = $slice_adaptor->fetch_by_region ($coordinate_system,$chromosome,$subslice_start,$subslice_end,$strand);
		
		if ($_debug) {
		    print STDERR "length of the extracted region after intergenic processing, " . $slice_region->length . "\n";
		}
		
	    }
	    else {
		if ($_debug) {
		    print STDERR "no gene upstream.\n";
		}
	    }
	}
	
	if ($_debug) {
	    print STDERR "writing down the sequence...\n";
	}
	
	my $upstream_sequence = $slice_region->seq;
	my $seqobj = Bio::PrimarySeq->new (
					   -id    => $gene_stable_id,
					   -desc  => $geneId,
					   -seq   => $upstream_sequence,
					   );
	$sout->write_seq ($seqobj);
	
	if ($_debug) {
	    print STDERR "sequence done.\n";
	}
	
	if ($report_features) {
	    
	    if ($_debug) {
		print STDERR "getting the features...\n";
	    }
	    
	    my $gff_output = getFeatures ($slice_region, $gene_stable_id);
	    
	    if ($gff_output ne "") {
		
		if ($_debug) {
		    print STDERR "writing GFF output...\n";
		}
		    
		print GFF "$gff_output";
	    }
	    else {
		if ($_debug) {
		    print STDERR "no overlapping features\n";
		}
	    }
	} # End reporting Features
	
    } # End processing Ensembl gene objects
}

if ($report_features) {
    close GFF;
}

if ($_debug) {
    print STDERR "processing done.\n";
}

## Benchmarking code
my $time1 = new Benchmark;

if ($_debug) {
   print STDERR "Script Benchmark: ".timestr(timediff($time1,$time0))."\n";
}

##
# End
##

# Deprecated !!

sub convertEnsembl2GFF {
    my ($seqId, $ensembl_features, $comments) = @_;
    my $gff_output = "";

    foreach my $feature (@$ensembl_features) {

	my $feature_module_name = ref ($feature);
	$feature_module_name    =~ /\w+\:+\w+\:+(\w+)/;
	my $feature_type        = $1;

	if ($_debug) {
	    print STDERR "processing feature, " . $feature->stable_id . ", of type, $feature_type...\n";
	}

	my $feature_id  = $feature->stable_id;
	my $comments = "EnsemblIdentifier=$feature_id";

	# if Exon
	if ($feature_type eq "exon") {
	    
	    # my $transcript_id = ?;
	    # my $gene_id       = ?;
	    
	    my $evidences = $feature->get_all_supporting_features();
	    $evidences = [];
	    
	    if (@$evidences > 0) {
		$feature_module_name = ref ($evidences->[0]);
		$feature_module_name =~ /\w+\:+\w+\:+(\w+)/;
		$feature_type        = $1;
		
		foreach my $evidence (@$evidences) {
		    $gff_output .= parse_ensembl_feature ($seqId, $evidence, $feature_type, "", "ensembl");
		}	    
	    }
	}
	
	$gff_output .= parse_ensembl_feature ($seqId, $feature, $feature_type, $comments, "ensembl");
    }
    
    return $gff_output;
}


sub parse_ensembl_feature {
    my ($seqId, $feature, $feature_type, $comments, $source) = @_;
    my $gff_output = "";

    my $start  = $feature->start;
    my $end    = $feature->end;
    my $strand = $feature->strand;
    if ($strand == 1) {
	$strand = "+";
    }
    else {
	$strand = "-";
    }
    
    # Analysis data ?????
    
    my $score = "";
    my $frame = "";

    $gff_output .= "$seqId\t$source\t$feature_type\t$start\t$end\t$score\t$strand\t$frame\t$comments\n";

    return $gff_output;
    
}

sub parse_ensembl_sequence {
    my ($seqId, $start, $end, $sequence_type, $comments, $strand, $source) = @_;
    my $gff_output = "";

    if ($strand == 1) {
	$strand = "+";
    }
    else {
	$strand = "-";
    }
    
    my $score = "";
    my $frame = "";

    $gff_output .= "$seqId\t$source\t$sequence_type\t$start\t$end\t$score\t$strand\t$frame\t$comments\n";

    return $gff_output;
    
}


sub has_gene_upstream {
    my ($slice_region, $strand, $downstream_length, $tss) = @_;
    
    my $genes = $slice_region->get_all_Genes();
    # print STDERR @$genes . " genes upstream before processing\n";
    # foreach my $gene (@$genes) {
	# print STDERR "stable id, " . $gene->stable_id . "\n";
    # }
    
    if ($downstream_length > 0) {
	
	# In that case include the input gene, so don't take it into account
	
	if (@$genes > 1) {
	    return 1;
	}
    }
    else {
	if (@$genes > 0) {
	    return 1;
	}
    }
    
    return 0;
}

sub getIntergenicSequence {
    my ($slice_region, $strand, $downstream_length, $tss) = @_;
    
    # The slice coordinates are related to chromosome coordinates therefore slice->start is not equal to 1
    # But genes associated with the slice region are attached to it with slice->start = 1 therefore we need to shift the coordinate of slice->start
    
    my $subslice_start = 1;
    my $subslice_end   = $slice_region->end - $slice_region->start;
    # Shift also the tss coordinate
    $tss = $tss - $slice_region->start;
    
    if ($_debug) {
	print STDERR "slice start and end, " . $slice_region->start . ", " . $slice_region->end . "\n";
	print STDERR "sub slice start and end shifted, $subslice_start, $subslice_end\n";
	print STDERR "tss shifted, $tss\n";
    }

    # Get the end of the most downstream of the upstream genes...
    
    my $latest_gene;
    # forward strand
    my $previous_gene_end = 1;
    # reverse strand
    if ($strand != 1) {
	$previous_gene_end = $subslice_end;
    }

    my @upstream_genes = @{$slice_region->get_all_Genes};
    my @upstream_exons = @{$slice_region->get_all_Exons};

    if ($_debug) {
	if (@upstream_genes > 0) {
	    print STDERR @upstream_genes . " genes upstream...\n";
	}
	
	if (@upstream_exons > 0) {
	    print STDERR @upstream_exons . " exons upstream...\n";
	}
    }

    foreach my $gene (@upstream_genes) {

	if ($_debug) {
	    print STDERR "previous_gene_end, $previous_gene_end\n";
	}

	my $gene_start = $gene->start;
	my $gene_end   = $gene->end;

	if ($_debug) {
	    print STDERR "processing gene, " . $gene->stable_id . " at " . $gene_start . ".." . $gene_end . " on strand " . $gene->strand . "\n";
	}

	if ($strand == -1) {
	    if ($_debug) {
		print STDERR "reversing gene coordinates!!\n";
	    }

	    my $gene_end_tmp = $subslice_end - $gene_start;
	    $gene_start = $subslice_end - $gene_end;
	    $gene_end   = $gene_end_tmp;
	}

	if ($_debug) {
	    print STDERR "processing gene after reversing its coordinates, " . $gene->stable_id . " at " . $gene_start . ".." . $gene_end . " on strand " . $gene->strand . "\n";
	}

	if ($strand == 1) {
	    if (($gene_end > $previous_gene_end) && ($gene_end < ($tss - 1))) {

		if ($_debug) {
		    print STDERR "\tupdating previous_gene_end - strand 1 processing\n";
		}

		$previous_gene_end = $gene_end;
		$latest_gene = $gene;

		if ($_debug) {
		    print STDERR "previous_gene_end updated, $previous_gene_end\n";
		}
	    }
	}
	else {
	    if (($gene_start < $previous_gene_end) && ($gene_start > ($tss + 1))) {
		
		if ($_debug) {
		    print STDERR "\tupdating previous_gene_end - strand -1 processing\n";
		}

		$previous_gene_end = $gene_start;
		$latest_gene = $gene;
	    }
	}
	
    }

    if (not defined $latest_gene) {
	if ($_debug) {
	    print STDERR "problem, suppose to have a gene in the upstream region but can't find it !!\n";
	}
    }
    
    if ($_debug) {
	print STDERR "before, got subslice, start, end, $subslice_start, $subslice_end\n";
    }
    
    if ($strand == 1) {
	$subslice_start = $previous_gene_end + 1;
    }
    else {
	$subslice_end   = $previous_gene_end - 1;
    }

    if ($_debug) {
	print STDERR "after, got subslice, start, end, $subslice_start, $subslice_end\n";
    }

    # Now we shift back to chromosome coordinates

    if ($strand == 1) {
	$subslice_start = $subslice_start + $slice_region->start;
	$subslice_end   = $slice_region->end;
    }
    else {
	$subslice_start = $slice_region->start;
	$subslice_end   = $subslice_end + $slice_region->start;
    }

    if ($_debug) {
	print STDERR "back to chromosome coordinates, this is giving, $subslice_start..$subslice_end\n";
    }

    return ($subslice_start, $subslice_end);
}


sub getFeatures {
    my ($slice_region, $input_gene_id) = @_;
    my $gff_output = "";
    
    # Right now, only report gene and exon features...
    # Report also introns, and UTRs

    # Only those on the slice !! - Make sure that's the case...
    
    my @genes = ();
    my @genes_tmp = @{$slice_region->get_all_Genes()};
    
    # Filter the input gene
    
    if ($downstream_length > 0) {
	foreach my $gene (@genes_tmp) {
	    if ($gene->stable_id ne $input_gene_id) {
		push (@genes, $gene);
	    }
	}
    }
    else {
	@genes = @genes_tmp;
    }

    # Get the transcripts, exons, introns, UTRs information for each gene
    
    foreach my $gene (@genes) {

	# Process the Gene Feature

	# Check the coordinates... (no negative ones)

	my $geneId      = $gene->stable_id;
	# source is only available since Ensembl release 31
	# Default source is "ensembl"
	my $source      = "ensembl";
	
	if ($releaseNumber >= 31) {
	    $source = $gene->source;
	}

	my $gene_strand = $gene->strand;
	my $comments    = "EnsemblIdentifier=$geneId";
	my $biotype     = $gene->type;
	if ($releaseNumber >= 31) {
	    # Since Ensembl release 31, type is deprecated and has been replaced by biotype
	    $biotype = $gene->biotype;
	}
	my $RNA_type;
        SWITCH: {
	    if ($biotype eq "protein_coding") { $RNA_type = "mRNA";   last SWITCH; }
	    if ($biotype =~ "rRNA")           { $RNA_type = "rRNA";   last SWITCH; }
	    if ($biotype =~ "tRNA")           { $RNA_type = "tRNA";   last SWITCH; }
	    if ($biotype =~ "ncRNA")          { $RNA_type = "ncRNA";  last SWITCH; }
	    if ($biotype =~ "snRNA")          { $RNA_type = "snRNA";  last SWITCH; }
	    if ($biotype =~ "snoRNA")         { $RNA_type = "snoRNA"; last SWITCH; }
	    if ($biotype =~ "scRNA")          { $RNA_type = "scRNA";  last SWITCH; }
	    if ($biotype =~ "miRNA")          { $RNA_type = "miRNA";  last SWITCH; }
	    # ensembl type returned for protein coding transcript for release before v31 !!!
	    if ($biotype =~ "ensembl")          { $RNA_type = "mRNA";  last SWITCH; }
	    # Default is to map biotype as a RNA type
	    print STDERR "don't know anything about this gene biotype, $biotype!\n";
	    
	    if ($biotype =~ /^(\w+)_pseudogene/) {
		$RNA_type = $biotype;
	    }
	    else {
		$RNA_type = $biotype;
	    }
	}

	my $strand_info;
	if ($gene_strand == 1) {
	    $strand_info = "Forward";
	}
	else {
	    $strand_info = "Reverse";
	}
	
	my $GFF_feature = parse_ensembl_feature ($input_gene_id, $gene, "gene", $comments ,$source);
	$gff_output    .= $GFF_feature;

        # Process the Transcript Features
	
	my $transcripts = $gene->get_all_Transcripts;

	if ($RNA_type eq "mRNA") {
	    # Report the number of transcripts of the protein coding genes
	    $gff_output    .= "# Gene $geneId ($strand_info). " . @$transcripts . " transcripts.\n";
	}
	
	foreach my $transcript (@$transcripts) {
	    my $transcriptId = $transcript->stable_id;
	    $comments        = "EnsemblIdentifier=$transcriptId; GeneIdentifier=$geneId";
	    $GFF_feature     = parse_ensembl_feature ($input_gene_id, $transcript, $RNA_type, $comments, $source);
	    $gff_output     .= $GFF_feature;
	    
	    # Process the Transcripts subfeatures
	    
	    my $exons = $transcript->get_all_translateable_Exons;
	    if (@$exons > 0) {

		# Gene 1 (Forward). 8 exons. 360 aa.
	    
		my $pep_length  = $transcript->translate->length;
		$gff_output    .= "# Gene $transcriptId ($strand_info). " . @$exons . " exons. $pep_length aa\n";

		# Deal with protein coding genes
		
		# Check we have the UTRs (something they're not predicted !! therefore the gene start end = CDS start end)
		# Get the UTRs start/end
		
		my $five_UTR_start;
		my $five_UTR_end;

		my $three_UTR_start;
		my $three_UTR_end;
		
		my $has_UTRs = 0;
		if (($gene->start < $transcript->coding_region_start) && ($gene->end > $transcript->coding_region_end)) {
		    $has_UTRs++;
		}

		if ($has_UTRs) {
		    if ($gene_strand == 1) {
			$five_UTR_start  = $gene->start;
			$five_UTR_end    = $transcript->coding_region_start - 1;
			
			$three_UTR_start = $transcript->coding_region_end + 1;
			$three_UTR_end   = $gene->end;
		    }
		    else {
			$five_UTR_end    = $gene->end;
			$five_UTR_start  = $transcript->coding_region_end + 1;
			
			$three_UTR_end   = $transcript->coding_region_start - 1;
			$three_UTR_start = $gene->start;
		    }
		    
		    my $five_UTR     = $transcript->five_prime_utr;
		    $comments        = "TranscriptIdentifier=$transcriptId; GeneIdentifier=$geneId";
		    $GFF_feature     = parse_ensembl_sequence ($input_gene_id, $five_UTR_start, $five_UTR_end, "5'UTR", $comments, $gene_strand, $source);
		    $gff_output     .= $GFF_feature;
		}

		foreach my $exon (@$exons) {
		    my $exonId   = $transcript->stable_id;
		    $comments    = "EnsemblIdentifier=$exonId; TranscriptIdentifier=$transcriptId; GeneIdentifier=$geneId";
		    $GFF_feature = parse_ensembl_feature ($input_gene_id, $exon, "exon", $comments, $source);
		    $gff_output .= $GFF_feature;
		}
		
		if (@$exons > 1) {
		    # More than one exon, we've got at least one intron
		    my $introns      = $transcript->get_all_Introns;
		    foreach my $intron (@$introns) {
			$comments    = "TranscriptIdentifier=$transcriptId; GeneIdentifier=$geneId";
			$GFF_feature = parse_ensembl_feature ($input_gene_id, $intron, "intron", $comments, $source);
			$gff_output .= $GFF_feature;
		    }
		}

		if ($has_UTRs) {
		    my $three_UTR    = $transcript->three_prime_utr;
		    $comments        = "TranscriptIdentifier=$transcriptId; GeneIdentifier=$geneId";
		    $GFF_feature     = parse_ensembl_sequence ($input_gene_id, $three_UTR_start, $three_UTR_end, "3'UTR", $comments, $gene_strand, $source);
		    $gff_output     .= $GFF_feature;
		}

	    }
	}
    }

    return $gff_output;
}

