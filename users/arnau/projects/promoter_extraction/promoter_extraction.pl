#!/usr/local/bin/perl -w

# NB: Start of gene definition: TSS or if not known, start codon

# Issue warnings about suspicious programming.
use warnings 'all';

# Must declare and initialize all variables
use strict;

# be prepare for command-line options/arguments
use Getopt::Std;

use Data::Dumper;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

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

END_HELP

}

BEGIN {
	
    # Determines the options with values from program
    # use vars qw/$opt_h $opt_f $opt_s $opt_r $opt_u $opt_d $opt_i $opt_g/;
    use vars qw /%opts/;
    
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

my $_debug = 1;

my $geneIds_file;
my $species;
my $release;
my $upstream_length   = 2000;
my $downstream_length = 500;
my $report_features   = 0;
my $intergenic_only   = 0;

defined $opts{f} and $geneIds_file      = $opts{f};
defined $opts{s} and $species           = $opts{s};
defined $opts{r} and $release           = $opts{r};
defined $opts{u} and $upstream_length   = $opts{u};
defined $opts{d} and $downstream_length = $opts{d};
defined $opts{g} and $report_features++;
defined $opts{i} and $intergenic_only++;

print STDERR "report features, $report_features, intergenic only, $intergenic_only\n";

# Ensembl database configuration

my $host      = "ensembldb.ensembl.org";
my $user      = "anonymous";
my $dbensembl = $species . "_core_" . $release;

$release =~ /^(\d+)_\w+$/;
my $branch = $1;

print STDERR "branch, $branch\n";

my $ensembl_API_path = "/home/ug/arnau/cvs/ensembl-$branch/modules";

print STDERR "ensembl_API_path, $ensembl_API_path.\n";

if (-d $ensembl_API_path) {
    use lib "$ensembl_API_path";
}
else {
    print STDERR "release, $release, not supported, contact gmaster\@imim.es for help\n";
    exit 1;
}
print STDERR "connecting to Ensembl database...\n";

my $dbh = new Bio::EnsEMBL::DBSQL::DBAdaptor (
					      -host => $host,
					      -user => $user,
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

print STDERR "processing the " . @geneIds . " gene identifiers...\n";

foreach my $geneId (@geneIds) {

    print STDERR "\nprocessing gene identifier, $geneId...\n";

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
    
    if (@$genes > 0) {
	if (@$genes > 1) {
	    print STDERR "Found more than one Ensembl gene (" . @$genes . ") associated with external identifier, $geneId!!!\n";
	    print STDERR "Will return the upsteam region of all Ensembl genes\n";
	    
	}
    }
    
    # It must be an Ensembl identifier

    if ((not defined $genes) || (@$genes == 0)) {
	my $gene = $gene_adaptor->fetch_by_stable_id ($geneId) or warn "can't instanciate gene object with identifier, $geneId!\n";
	if (not defined $gene) {
	    print STDERR "can't find any gene for this given identifier, $geneId!\n";
	}
	else {
	    print STDERR "gene stable identifier found: " . $gene->stable_id . ", on strand, " . $gene->strand . "\n";
	    $genes = [ $gene ];
	}
    }

    foreach my $gene (@$genes) {

	my $gene_stable_id = $gene->stable_id;
	my $dblinks = $gene->get_all_DBLinks ();
	
	# print STDERR "db links, " . Dumper ($dblinks) . "\n";
	
	my $coordinate_system = $gene->coord_system_name();
	my $chromosome        = $gene->seq_region_name();
	my $strand            = $gene->strand;
	my $start             = "";
	my $end               = "";
	my $features = [];
	
	if ($strand > 0) {
	    $start = $gene->start - $upstream_length;
	    $end   = $gene->start + $downstream_length;
	}
	else {
	    $end   = $gene->end + $upstream_length;
	    $start = $gene->end - $downstream_length;
	}
	
	# Get rid of the base 0 if we don't want any downstream region, this way we won't get any overlapping feature associated with the given gene (its first exon for example)
	
	if (($downstream_length == 0) && ($strand == 1)) {
	    $end -= 1;
	}
	elsif (($downstream_length == 0) && ($strand == -1)) {
	    $start += 1;
	}
	
	# Get the Slice object corresponding to the region of the genome we want to extract
	
	my $slice_region = $slice_adaptor->fetch_by_region ($coordinate_system,$chromosome,$start,$end,$strand);
	
	if ($intergenic_only) {
	    if (has_gene_upstream ($slice_region, $strand, $downstream_length)) {
		
		print STDERR "has gene upstream\n";
		
		my ($subslice_start, $subslice_end) = getIntergenicSequence ($slice_region, $strand, $downstream_length);
		
		print STDERR "subslice start, end: $subslice_start, $subslice_end\n";

		# Update the slice object so it inlcudes only the intergenic sequence
		
		$slice_region = $slice_adaptor->fetch_by_region ($coordinate_system,$chromosome,$subslice_start,$subslice_end,$strand);
	    }
	}
	
	my $upstream_sequence = $slice_region->seq;
	my $seqobj = Bio::PrimarySeq->new (
					   -id    => $gene_stable_id,
					   -desc  => $geneId,
					   -seq   => $upstream_sequence,
					   );
	$sout->write_seq ($seqobj);
	
	if ($report_features) {
	    
	    my $upstream_features = getFeatures ($slice_region, $strand, $downstream_length);
	    
	    if (@$upstream_features > 0) {
		
		my $gff_output = convertEnsembl2GFF ($geneId, $upstream_features);
		
		if ($gff_output ne "") {
		
		    print STDERR "writing GFF output...\n";
		    
		    print GFF "$gff_output";
		}
		else {
		    print STDERR "Warning: no GFF output returned!!\n";
		}
		
	    }
	    else {
		print STDERR "no overlapping features\n";
	    }
	} # End reporting Features

    } # End processing Ensembl gene objects
}

if ($report_features) {
    close GFF;
}

print STDERR "processing done.\n";

##
# End
##

sub convertEnsembl2GFF {
    my ($seqId, $ensembl_features) = @_;
    my $gff_output = "";

    foreach my $feature (@$ensembl_features) {

	# if Exon
	my $evidences = $feature->get_all_supporting_features();

	my $feature_module_name = ref ($evidences->[0]);
	$feature_module_name    =~ /\w+\:+\w+\:+(\w+)/;
	my $feature_type        = $1;

	foreach my $evidence (@$evidences) {
	    $gff_output .= parse_ensembl_feature ($seqId, $evidence, $feature_type);
	}

	# print STDERR "processing feature, " . Dumper ($feature) . "...\n";
	
	$feature_module_name = ref ($feature);
	$feature_module_name    =~ /\w+\:+\w+\:+(\w+)/;
	$feature_type        = $1;

	$gff_output .= parse_ensembl_feature ($seqId, $feature, $feature_type);
    }

    return $gff_output;
}


sub parse_ensembl_feature {
    my ($seqId, $feature, $feature_type) = @_;
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
    
    # algorithm
    my $algorithm = "Ensembl";
    # score
    my $score = "";
    # Comments
    my $comments = "";
    
    $gff_output .= "$seqId\t$algorithm\t$feature_type\t$start\t$end\t$score\t$strand\t$comments\n";

    return $gff_output;
    
}


sub has_gene_upstream {
    my ($slice_region, $strand, $downstream_length) = @_;
    
    print STDERR "length slice region, " . $slice_region->length()  . "\n";
    
    my $tss;
    if ($strand == 1) {
	$tss = $slice_region->end   - $downstream_length;
    }
    else {
	$tss    = $slice_region->start + $downstream_length;
    }

    my $exons1 = $slice_region->get_all_Exons();
    print STDERR @$exons1 . " exons upstream before processing\n";
    
    # Get rid of the downstream region if any
    
    my $slice_subregion = $slice_region;
    if ($downstream_length > 0) {
	if ($strand == 1) {
	    $slice_subregion = $slice_region->sub_Slice ($slice_region->start, $tss - 1, $strand);
	    $slice_subregion = $slice_region->sub_Slice ($tss + 1, $slice_region->end, $strand);
	}
    }
    
    print STDERR "length sub slice region, " . $slice_subregion->length()  . "\n";

    my $exons2 = $slice_subregion->get_all_Exons();
    print STDERR @$exons2 . " exons upstream after processing\n";

    if (@{$slice_subregion->get_all_Genes} > 0) {
	return 1;
    }

    return 0;
}

sub getIntergenicSequence {
    my ($slice_region, $strand, $downstream_length) = @_;

    print STDERR "strand, $strand\n";

    # Initialize the downstream coordinate - kept unchanged

    my $subslice_start;
    my $subslice_end;
    if ($strand == 1) {
	$subslice_start = $slice_region->start;
	$subslice_end   = $slice_region->end;
    }
    else {
	$subslice_start = $slice_region->start;
	$subslice_end   = $slice_region->end;
    }

    print STDERR "slice start and end, " . $slice_region->start . ", " . $slice_region->end . "\n";
    
    # Get the end of the most downstream of the upstream genes...
    
    my $tss;
    if ($strand == 1) {
	$tss = $slice_region->end - $downstream_length;
    }
    else {
	$tss    = $slice_region->start + $downstream_length;
    }

    my $latest_gene;
    # forward strand
    my $previous_gene_end = 0;
    # reverse strand
    if ($strand != 1) {
	$previous_gene_end = $slice_region->end;
    }

    my @upstream_genes = @{$slice_region->get_all_Genes};
    my @upstream_exons = @{$slice_region->get_all_Exons};
    
    if (@upstream_genes > 0) {
	print STDERR "genes upstream...\n";
    }
    
    if (@upstream_exons > 0) {
	print STDERR "exons upstream...\n";
    }
    
    foreach my $gene (@upstream_genes) {
	
	if ($strand == 1) {
	    if (($gene->end > $previous_gene_end) && ($gene->end < $tss)) {
		$previous_gene_end = $gene->end;
		$latest_gene = $gene;
	    }
	}
	else {
	    if (($gene->start < $previous_gene_end) && ($gene->start > $tss)) {
		$previous_gene_end = $gene->start;
		$latest_gene = $gene;
	    }
	}
	
    }

    if (not defined $latest_gene) {
	print STDERR "problem, suppose to have a gene in the upstream region but can't find it !!\n";
    }
    
    print STDERR "before, got subslice, start, end, $subslice_start, $subslice_end\n";

    if ($strand == 1) {
	$subslice_start = $previous_gene_end + 1;
    }
    else {
	$subslice_end   = $previous_gene_end - 1;
    }

    print STDERR "after, got subslice, start, end, $subslice_start, $subslice_end\n";

    return ($subslice_start, $subslice_end);
}


sub getFeatures {
    my ($slice_region, $strand, $downstream_length) = @_;
    my $features = [];
    
    if (not $intergenic_only) {
	my $exon_features = $slice_region->get_all_Exons();

	# Check we don't get the exons of the input gene itself !

	# ...

	push (@$features, @$exon_features);
    }

    # What other features do we want to report ??

    return $features;
}

