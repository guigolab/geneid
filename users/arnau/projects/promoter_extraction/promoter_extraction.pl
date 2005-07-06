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

# Bioperl
use Bio::SeqIO;
use Bio::PrimarySeq;

sub help {
return <<"END_HELP";
Description: Extract Promoter regions for a set of given genes
Usage:

    promoter_extraction.pl [-h] -f {gene ids} -s {Species} -r {Genome Release} -u {upstream} -d {downstream} -i {Flag intergenic regions only} -g {Flag report attached features in GFF format}
	-h help
	-f Gene Identifiers input file
	-s Gender species, e.g. homo_sapiens
	-r Database release, e.g 29_35b
	-u length of the upstream region   (Default 2000)
	-d length of the downstream region (Default 500)
	-i 
	-g

Examples using some combinations:
	promoter_extraction.pl -f geneIds.lst -s homo_sapiens -r 29_35b -u 2000 -d 500 -i -g
	promoter_extraction.pl -f refseqs.lst -s rattus_norvegicus -r 29_3f -u 2000 -d 500 -i -g

	promoter_extraction.pl -f geneIds.lst -s homo_sapiens -r 31_35d -u 2000 -d 500 -i -g

END_HELP

}

BEGIN {
	
    # Global variables definition
    use vars qw /$_debug %opts $report_features $intergenic_only $upstream_length $downstream_length $geneIds_file $dbhost $dbuser $dbensembl/;

    if (-f "config.pl") {
        require "config.pl";
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
my $branch = $1;

my $ensembl_API_path = "./lib/ensembl-$branch/modules";

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
	    
	    my $upstream_features = getFeatures ($slice_region, $gene_stable_id);
	    
	    if (@$upstream_features > 0) {
		
		my $gff_output = convertEnsembl2GFF ($geneId, $upstream_features);
		
		if ($gff_output ne "") {
		    
		    if ($_debug) {
			print STDERR "writing GFF output...\n";
		    }
		    
		    print GFF "$gff_output";
		}
		else {
		    if ($_debug) {
			print STDERR "Warning: no GFF output returned!!\n";
		    }
		}
		
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

##
# End
##

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
		    $gff_output .= parse_ensembl_feature ($seqId, $evidence, $feature_type, "");
		}	    
	    }
	}
	
	$gff_output .= parse_ensembl_feature ($seqId, $feature, $feature_type, $comments);
    }
    
    return $gff_output;
}


sub parse_ensembl_feature {
    my ($seqId, $feature, $feature_type, $comments) = @_;
    my $gff_output = "";

    # Get the gene name for each exon and report it.

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
    
    # algorithm
    my $algorithm = "Ensembl";
    # score
    my $score = "";

    $gff_output .= "$seqId\t$algorithm\t$feature_type\t$start\t$end\t$score\t$strand\t$comments\n";

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
    my $features = [];
    
    # Right now, only report gene and exon features...
    # What else ?

    my @genes = @{$slice_region->get_all_Genes()};
    
    # Filter the input gene
    
    if ($downstream_length > 0) {
	foreach my $gene (@genes) {
	    if ($gene->stable_id ne $input_gene_id) {
		push (@$features, $gene);
		# Get the exons
		my $exons = $gene->get_all_Exons;
		push (@$features, @$exons);
	    }
	}
    }
    else {
	push (@$features, @genes);
	foreach my $gene (@genes) {
	    my $exons = $gene->get_all_Exons;
	    push (@$features, @$exons);
	}
    }

    # What other features do we want to report ??
    
    return $features;
}

