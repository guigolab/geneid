#!/usr/local/bin/perl -w

=head1 AUTHOR

Arnaud Kerhornou, akerhornou@imim.es

=head1 COPYRIGHT

Copyright (c) 2006, Arnaud Kerhornou and INB - Nodo Vertical 1 INB/CRG.
 All Rights Reserved.

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.


=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

##################################################################
#
# GenericSequence to FASTA conversion Moby Service Client
#
##################################################################

use strict;
use Data::Dumper;

# be prepare for command-line options/arguments
use Getopt::Std;

# BioMoby and SOAP libraries

use MOBY::Client::Central;
use MOBY::Client::Service;
# use SOAP::Lite + 'trace';

# XML DOM parsing library
use XML::LibXML;

# MIME encoding/decoding
use MIME::Base64;

# Benchmark Module
use Benchmark;

# HTML encoding/decoding module
use HTML::Entities;

# GEPAS CGI calls support
use LWP::UserAgent;
use HTTP::Request::Common;

##################################################################

sub help {
return <<"END_HELP";
Description: Execute a gene clustering workflow, based on patterns found in the upstream regions of a given set of genes. This workflow takes a set of gene upstream sequences in FASTA format and return in STDOUT a clustering tree picture in PNG format.
Usage:

    GenesClustering_FASTA.pl [-h] -x {Moby Central} -f {sequence FASTA file} -t {MatScan threshold} -d {MatScan database} -a {meta-alignment alpha penalty} -l {meta-alignment lambda penalty} -u {meta-alignment mu penalty} -g {multiple meta-alignment gamma penalty} -r {non-collinearity penalty} -m {Hierarchical clustering method} -n {Number of K-means clusters} -i {Number of K-means iterations} -o {Output directory}
	-h help
	-x MOBY Central: Inab, BioMoby, Mobydev (optional - Default is BioMoby registry)
		<1> or Inab
		<2> or BioMoby
		<3> or Mobydev
	-f Sequence(s) input file, in FASTA format
	-t MatScan probability threshold (Default is 0.85)
        -d MatScan Motifs database [Jaspar, Transfac] (Default is Transfac)
	-a Meta-alignment alpha penalty (Default is 0.5)
	-l Meta-alignment lambda penalty (Default is 0.1)
	-u Meta-alignment mu penalty (Default is 0.1)
	-g Multiple meta-alignment gamma penalty (Default is -10)
	-r Multiple meta-alignment non-collinearity penalty (Default is 100)
        -m HierarchicalCluster method, e.g nearest neighbour joining or furthest neighbour joining [nearest, furthest] (Default is nearest)
	-n Number of clusters returned by the k-means algorithm (Default is 10)
	-i Number of k-means iterations (Default is 200)
        -o Output directory name, if not specified, the output is turned off, the script will just return a tree clustering picture in STDOUT.
	-c workflow configuration file (default is \$HOME/.workflow.config)

Examples using some combinations:
	perl GenesClustering_FASTA.pl -x 2 -f /home/ug/arnau/data/ENSRNOG00000007726.orthoMode.withRat.1000.fa -c \$HOME/.workflow.config -t 0.80 -d jaspar -a 0.5 -l 0.1 -u 0.1 -g -10 -r 100 -m nearest -n 4 -i 200 -o output

END_HELP

}

BEGIN {
	
    # Determines the options with values from program
    use vars qw/$opt_h $opt_x $opt_f $opt_c $opt_t $opt_d $opt_a $opt_l $opt_u $opt_m $opt_g $opt_r $opt_n $opt_i $opt_o $opt_s $output_dir/;
    
    # these are switches taking an argument (a value)
    my $switches = 'x:shf:c:t:d:a:l:u:m:g:r:n:i:o:';
    
    # Get the switches
    getopts($switches);
    
    # If the user does not write nothing, skip to help
    if (defined($opt_h) || !defined ($opt_f)){
	print STDERR help;
	if (defined $output_dir) {
	    print 1;
	}
	exit 1;
    }
    
}

my $t1 = Benchmark->new ();

my $_debug = 0;
# Need meta-alignment software because it is run locally, in case there are too many input sequences!
my $_meta_dir = "/usr/local/molbio/bin";
my $_meta_bin = "meta";

##################################################################
#
# Setup sequences input data and moby parameters
#
##################################################################

############################
#
# Hidden option - a hacky one !!
#
# use this (-s) flag when you have too many sequences, and you've done manually the multi pairwise alignments, and you want to start the workflow at the score matrix generation step!!
#
############################

my $shortcut = 0;
defined ($opt_s) && $shortcut++;

if ($shortcut) {
    print STDERR "shortcut activated, will start the workflow at 'fromMetaAlignmentsToTextScoreMatrix' step!\n";
}

# input file

my $in_file = $opt_f;
if (! defined $in_file) {
    print STDERR "Error, please specify an input file with promoter FASTA sequences\n";
    if (defined $output_dir) {
	print 1;
    }
    exit 1;
}
if (not (-f $in_file)) {
    print STDERR "Error, can't find input file, $in_file\n";
    if (defined $output_dir) {
	print 1;
    }
    exit 1;
}

# Command-line parameters

# NJ parameters

my $method    = $opt_m || "nearest";
if (lc ($method) eq "nearest") {
  $method = "Nearest neighbor (single linkage)";
}
elsif (lc ($method) eq "furthest") {
  $method = "Furthest neighbor (complete linkage)";
}
else {
  print STDERR "don't know about the value for method parameter, $method!\n";
  print STDERR help;
  if (defined $output_dir) {
      print 1;
  }
  exit 1;
}

my $method_xml = "<Value>$method</Value>";

# MatScan parameters

my $threshold = $opt_t || 0.85;
my $database  = $opt_d || "transfac";

my $threshold_xml   = "<Value>$threshold</Value>";
my $matrix_xml      = "<Value>$database</Value>";

# Pairwise Meta-alignment parameters

my $alpha  = $opt_a || 0.5;
my $lambda = $opt_l || 0.1;
my $mu     = $opt_u || 0.1;

my $alpha_xml  = "<Value>$alpha</Value>";
my $lambda_xml = "<Value>$lambda</Value>";
my $mu_xml     = "<Value>$mu</Value>";

# K-means parameters

my $cluster_number       = $opt_n || 10;
my $cluster_number_xml   = "<Value>$cluster_number</Value>";
my $iteration_number     = $opt_i || 200;
my $iteration_number_xml = "<Value>$iteration_number</Value>";

# Multiple Meta-alignment parameters

# Same than Pairwise Meta-alignment parameters, plus:
my $gamma            = $opt_g || -10;
my $non_collinearity = $opt_r || 100;

my $gamma_xml = "<Value>$gamma</Value>";
my $non_collinearity_xml = "<Value>$non_collinearity</Value>";

# parameters configuration file parsing

my $serviceName = undef;

my $config_file = $opt_c || $ENV{HOME} . "/.workflow.config";
if (not (-f $config_file)) {
    print STDERR "Error, can't find config file, $config_file\n";
    if (defined $output_dir) {
	print 1;
    }
    exit 1;
}

my %parameters = setConfigurationData ($config_file);

# Output
$output_dir = $opt_o || undef;
if (defined $output_dir) {
  if (!$shortcut) {
      if (-d $output_dir) {
	  qx/rm -rf $output_dir/;
      }
      qx/mkdir $output_dir/;
  }
}

##################################################################
#
# Setup Moby configuration parameters
#
##################################################################

my $namespace = $parameters{input}->{namespace};
if ($_debug) {
  print STDERR "namespace, $namespace\n";
}
my $input_type = $parameters{input}->{input_type};
# Hardcoded it here !
$input_type = "identifiers";

if (! ($input_type eq "FASTA" || $input_type eq "identifiers")) {
	print STDERR "don't know about this input type, $input_type!\n";
	print STDERR "must be either 'FASTA' or 'identifiers'\n";
	if (defined $output_dir) {
	    print 1;
	}
	exit 1;
}

##############################################
#
# ASSIGN THE MOBY URI AND MOBY SERVER
#
##############################################

# Default is Icapture production registry

my $URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://mobycentral.icapture.ubc.ca/cgi-bin/MOBY05/mobycentral.pl';
my $URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://mobycentral.icapture.ubc.ca/MOBY/Central';
my $PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

if (defined($opt_x)) {
    
	# Delete spaces
        $opt_x =~ s/\s//g;
    
	# Assign the MOBY Server and MOBY URI
	if (($opt_x == 1) || (lc ($opt_x) eq 'inab')) {

	    $URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://www.inab.org/cgi-bin/MOBY-Central.pl';
	    $URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://www.inab.org/MOBY/Central';
	    $PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

	}elsif (($opt_x == 2) || (lc ($opt_x) eq 'biomoby')) {

	    $URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://mobycentral.icapture.ubc.ca/cgi-bin/MOBY05/mobycentral.pl';
	    $URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://mobycentral.icapture.ubc.ca/MOBY/Central';
	    $PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';
		
	}elsif (($opt_x == 3) || (lc ($opt_x) eq 'mobydev')) {

	    $URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://moby-dev.inab.org/cgi-bin/MOBY-Central.pl';
	    $URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://moby-dev.inab.org/MOBY/Central';
	    $PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

	}else {
	    print STDERR "Don't anything about this registry, $opt_x!\n";
	    print help;
	    if (defined $output_dir) {
		print 1;
	    }
	    exit 1;
	}

}

##################################################################
#
# Moby Server Instanciation
#
##################################################################

if ($_debug) {
    print STDERR "Instanciating Central object...\n";
}

my $C = MOBY::Client::Central->new(
				   Registries => {mobycentral => {URL => $URL,URI => $URI}}
				   );

if (defined $C) {
    if ($_debug) {
	print STDERR "Instanciation done.\n";
	print STDERR "TESTING MOBY CLIENT with\n\tURL: $URL\n\tURI: $URI\n\tProxy: $PROXY\n\n";
    }
}
else {
    print STDERR "Error, Central could not be instanciated!\n";
    if (defined $output_dir) {
	print 1;
    }
    exit 1;
}

##################################################################
#
# Setup the moby objects
#
##################################################################

my $input_xml;

if ($input_type eq "FASTA") {

    # fromFASTAToDNASequenceCollection

    my $input = qx/cat $in_file/;
    
    $input_xml = <<PRT;
<FASTA_NA_multi namespace="$namespace" id="">
<String id="" namespace="" articleName="content">
<![CDATA[
$input
]]>
</String>
</FASTA_NA_multi>
PRT

}
else {
    
    # getUpstreamSeqFromEnsembl
    
    my $input = qx/cat $in_file/;
    
    $input_xml = <<PRT;
<List_Text namespace="$namespace" id="">
<String id="" namespace="" articleName="content">
<![CDATA[
$input
]]>
</String>
</List_Text>
PRT

}

if ($_debug) {
    print STDERR "input xml,\n $input_xml\n";
}

# 1/ fromFASTAtoDNASequenceCollection or getUpstreamSeqFromEnsembl

if ($input_type eq "FASTA") {
	# fromFASTAtoDNASequenceCollection
	$serviceName = "fromFASTAToDNASequenceCollection";
}
else {
	# getUpstreamSeqFromEnsembl
	$serviceName = "getUpstreamSeqFromEnsembl";
}

# Get the parameters from the configuration file

my $authURI           = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
my $articleName       = $parameters{$serviceName}->{articleName} || die "no article name for $serviceName\n";
my $output_datatype   = $parameters{$serviceName}->{outputType}  || die "no output type for $serviceName\n";
my $species_xml           = "";
my $upstream_length_xml   = "";
my $downstream_length_xml = "";
my $intergenic_only_xml   = "";

if ($input_type eq "identifiers") {
    my $species           = $parameters{$serviceName}->{species} || die "no species\n";
    my $upstream_length   = $parameters{$serviceName}->{upstream_length}   || die "no upstream length\n";
    my $downstream_length = $parameters{$serviceName}->{downstream_length} || 0;
    my $intergenic_only   = $parameters{$serviceName}->{intergenic_only}   || die "no intergenic only\n";
    
    $species_xml           = "<Value>$species</Value>";
    $upstream_length_xml   = "<Value>$upstream_length</Value>";
    $downstream_length_xml = "<Value>$downstream_length</Value>";
    $intergenic_only_xml   = "<Value>$intergenic_only</Value>";
}

##################################################################
#
# Moby Service instantiation
#
##################################################################

if ($_debug) {
    print STDERR "finding service, $serviceName...\n";
}

my $Service = getService ($C, $serviceName, $authURI);

##################################################################
#
# Service execution
#
##################################################################

my $moby_response;
if ($input_type eq "FASTA") {
    $moby_response = $Service->execute (XMLinputlist => [
							 ["$articleName", $input_xml],
							]);
}
else {

    print STDERR "Preliminary step, upstream sequence retrieval from Ensembl...\n\n";
    
    $moby_response = $Service->execute (XMLinputlist => [
							 ["$articleName", $input_xml, "species", $species_xml, "upstream length", $upstream_length_xml, "downstream length", $downstream_length_xml, "intergenic only", $intergenic_only_xml]
							]);
}

##################################################################
#
# Moby_Response processing
#
##################################################################

if ($_debug) {
	print STDERR "$serviceName results\n";
	print STDERR $moby_response;
	print STDERR "\n";
}

if (hasFailed ($moby_response)) {
    print STDERR "service, $serviceName, has failed!\n";
    my $moby_error_message = getExceptionMessage ($moby_response);
    print STDERR "reason is the following,\n$moby_error_message\n";
    if (defined $output_dir) {
	print 1;
    }
    exit 1;
}

$input_xml = parseResults ($moby_response, $output_datatype);

# Convert into FASTA if necessary

if (($input_type ne "FASTA") && (defined $output_dir)) {
    my $Conversion_Service = getService ($C, "fromGenericSequenceCollectionToFASTA", $authURI);
    my $moby_fasta_response = $Conversion_Service->execute (XMLinputlist => [
									     ["sequences", $input_xml]
									     ]);
    saveResults ($moby_fasta_response, "FASTA", "upstream_sequences", $output_dir);
}

if ($_debug) {
    print STDERR "input xml for next service:\n";
    print STDERR join (', ', @$input_xml);
    print STDERR ".\n";
}

my %matscan_results;
my $moby_matscan_response;

if (!$shortcut) {
    
    # runMatScanGFFCollection
    
    print STDERR "First step, $database binding sites predictions...\n";
    print STDERR "Executing MatScan...\n";
    
    $serviceName        = "runMatScanGFFCollection";
    $authURI            = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
    $articleName        = $parameters{$serviceName}->{articleName} || die "article name for $serviceName\n";
    
    $Service = getService ($C, $serviceName, $authURI);
    
    $moby_response = $Service->execute (XMLinputlist => [
							 ["$articleName", $input_xml, "threshold", $threshold_xml, "motif database", $matrix_xml]
							 ]);
    
    if ($_debug) {
	print STDERR "$serviceName results\n";
	print STDERR $moby_response;
	print STDERR "\n";
    }
    
    if (hasFailed ($moby_response)) {
	print STDERR "service, $serviceName, has failed!\n";
	my $moby_error_message = getExceptionMessage ($moby_response);
	print STDERR "reason is the following,\n$moby_error_message\n";
	if (defined $output_dir) {
	    print 1;
	}
	exit 1;
    }
    
    $input_xml = parseResults ($moby_response, "GFF");
    # Setup a hastable to fast access and keep it for later !
    %matscan_results = setMatScan_hash ($input_xml);
    $moby_matscan_response = $moby_response;
    my @sequence_matscan_ids = keys (%matscan_results);
    
    if (defined $output_dir) {
	saveResults ($moby_response, "GFF", "MatScan", $output_dir);
    }
    
    if ($_debug) {
	print STDERR "input xml for next service:\n";
	print STDERR join (', ', @$input_xml);
	print STDERR ".\n";
    }
    
    print STDERR "First step done\n\n";
    
    # runMultiPairwiseMetaAlignmentGFF & runMultiPairwiseMetaAlignment
    
    print STDERR "Second step, making the pairwise alignments of the binding site maps, using meta-alignment...\n";
    
    if (@$input_xml < 50) {
	
	print STDERR "Invoking meta-alignment remotely!\n";
	
	# runMultiPairwiseMetaAlignmentGFF first
	
	# Run this service just to save the results in GFF format
	
	if (defined $output_dir) {
	    
	    $serviceName = "runMultiPairwiseMetaAlignmentGFF";
	    $authURI     = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
	    $articleName = $parameters{$serviceName}->{articleName} || die "article name for $serviceName\n";
	    
	    $Service = getService ($C, $serviceName, $authURI);
	    
	    $moby_response = $Service->execute (XMLinputlist => [
								 ["$articleName", $input_xml, 'alpha', $alpha_xml, 'lambda', $lambda_xml, 'mu', $mu_xml]
								 ]);
	    
	    if ($_debug) {
		print STDERR "$serviceName result\n";
		print STDERR $moby_response;
		print STDERR "\n";
	    }
	    
	    if (hasFailed ($moby_response)) {
		print STDERR "service, $serviceName, has failed!\n";
		my $moby_error_message = getExceptionMessage ($moby_response);
		print STDERR "reason is the following,\n$moby_error_message\n";
		if (defined $output_dir) {
		    print 1;
		}
		exit 1;
	    }
	    
	    saveResults ($moby_response, "GFF", "Meta", $output_dir);
	}
	
	if ($_debug) {
	    print STDERR "input xml for next service:\n";
	    print STDERR join (', ', @$input_xml);
	    print STDERR ".\n";
	}
	
	# Then runMultiPairwiseMetaAlignment
	
	$serviceName = "runMultiPairwiseMetaAlignment";
	$authURI     = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
	$articleName = $parameters{$serviceName}->{articleName} || die "article name for $serviceName\n";
	
	$Service = getService ($C, $serviceName, $authURI);
	
	$moby_response = $Service->execute (XMLinputlist => [
							     ["$articleName", $input_xml, 'alpha', $alpha_xml, 'lambda', $lambda_xml, 'mu', $mu_xml]
							     ]);
	
	if ($_debug) {
	    print STDERR "$serviceName result\n";
	    print STDERR $moby_response;
	    print STDERR "\n";
	}
	
	if (hasFailed ($moby_response)) {
	    print STDERR "service, $serviceName, has failed!\n";
	    my $moby_error_message = getExceptionMessage ($moby_response);
	    print STDERR "reason is the following,\n$moby_error_message\n";
	    if (defined $output_dir) {
		print 1;
	    }
	    exit 1;
	}
	
	$input_xml = parseResults ($moby_response, "Meta_Alignment_Text");
	if (defined $output_dir) {
	    saveResults ($moby_response, "Meta_Alignment_Text", "Meta", $output_dir);
	}
	
	if ($_debug) {
	    print STDERR "input xml for next service:\n";
	    print STDERR join (', ', @$input_xml);
	    print STDERR ".\n";
	}
	
    } # End running Multi Pairwise alignments Web services
    else {
	
	print STDERR "Too many maps to process, thus invoking meta-alignment locally\n";
	
	my $matscan_maps = parseTextContent ($moby_matscan_response, "GFF");
	
	if ($_debug) {
	    print STDERR "got " , @$matscan_maps . " maps!\n";
	}
	
	# Make the pairwise alignments
	
	my $i = 0;
	my $j = 0;
	my $meta_index = 1;
	$input_xml     = [];
	qx/mkdir $output_dir\/Meta/;
	
	foreach my $id1 (@sequence_matscan_ids) {
	    my $file_map1 = "$output_dir/MatScan/$id1.MatScan.gff";
	    
	    if (! -f $file_map1) {
		print STDERR "Error, can't find map file, $file_map1!\n";
	    }
	    
	    $j = 0;
	    foreach my $id2 (@sequence_matscan_ids) {
		my $file_map2 = "$output_dir/MatScan/$id2.MatScan.gff";
		
		if (! -f $file_map2) {
		    print STDERR "Error, can't find map file, $file_map2!\n";
		}
		
		if ($i < $j) {
		    my $meta_data = qx/$_meta_dir\/$_meta_bin -a 1 -l 0.1 -m 1 $file_map1 $file_map2/;
		    my $metadata_xml = "<Meta_Alignment_Text namespace='' id=''><String namespace='' id='' articleName='Content'><![CDATA[\n" . $meta_data . "]]></String></Meta_Alignment_Text>";
		    push (@$input_xml, $metadata_xml);
		    
		    # Also Store it in a file
		    
		    my $meta_filename = "$output_dir/Meta/$meta_index.Meta.txt";
		    qx/echo "$meta_data" > $meta_filename/;
		    
		    $meta_index++;
		}
		
		$j++;
	    }
	    
	    $i++;
	} # End performing the pairwise alignments
	
    } # End running Multi Pairwise meta-alignments locally
    
    print STDERR "Second step done!\n\n";
    
}
else {
    
    # Shortcup path...
    
    # Get the results from the Meta output directory
    
    $input_xml = [];
    
    print STDERR "parsing meta-alignment data in $output_dir/Meta directory...\n";
    
    opendir METADIR, "$output_dir/Meta";
    my @metafiles = grep /\.meta|\.txt/, readdir METADIR;
    
    if (@metafiles < 1) {
	print STDERR "Error, can't parse any meta-alignment data files in $output_dir/Meta!\n";
	closedir METADIR;
	exit 1;
    }
    closedir METADIR;
    
    if ($_debug) {
	print STDERR "got " . @metafiles . " meta files\n";
	print STDERR "(first one, " . $metafiles[0] . ")\n";
    }
    
    foreach my $metafile (@metafiles) {
	if (-f ("$output_dir/Meta/$metafile")) {
	    my $meta_data = qx/cat $output_dir\/Meta\/$metafile/;
	    my $metadata_xml = "<Meta_Alignment_Text namespace='' id=''><String namespace='' id='' articleName='Content'><![CDATA[\n" . $meta_data . "]]></String></Meta_Alignment_Text>";
	    push (@$input_xml, $metadata_xml);
	}
	else {
	    print STDERR "metafile, $output_dir/Meta/$metafile, doesn't return any data!\n";
	}
    }
    
    print STDERR "shortcut processing done!\n";
    
}

# fromMetaAlignmentsToTextScoreMatrix

print STDERR "Third step, generating a score matrix by parsing meta-alignment data...\n";

$serviceName = "fromMetaAlignmentsToTextScoreMatrix";
$authURI     = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
$articleName = $parameters{$serviceName}->{articleName} || die "article name for $serviceName\n";

$Service = getService ($C, $serviceName, $authURI);

$moby_response = $Service->execute (XMLinputlist => [
						     ["$articleName", $input_xml]
						     ]);

if ($_debug) {
	print STDERR "$serviceName results\n";
	print STDERR $moby_response;
	print STDERR "\n";
}

if (hasFailed ($moby_response)) {
    print STDERR "service, $serviceName, has failed!\n";
    my $moby_error_message = getExceptionMessage ($moby_response);
    print STDERR "reason is the following,\n$moby_error_message\n";
    if (defined $output_dir) {
	print 1;
    }
    exit 1;
}

my $input_xml_aref = parseResults ($moby_response, "MicroArrayData_Text");
if (defined $output_dir) {
  saveResults ($moby_response, "MicroArrayData_Text", "score_matrix", $output_dir);
}

if ($_debug) {
	print STDERR "input xml for next service:\n";
	print STDERR join (', ', @$input_xml_aref);
	print STDERR ".\n";
}

# convert the input xml into a scalar
$input_xml =  $input_xml_aref->[0];

print STDERR "Third step done!\n\n";

# runKMeansClustering

print STDERR "Fourth step, gene clustering using a k-means clustering algorithm...\n";

if ($_debug) {
  print STDERR "\nExecuting runKMeansClustering...\n\n";
}

$serviceName   = "runKMeansClustering";
$authURI       = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
$articleName   = $parameters{$serviceName}->{articleName} || "";

$Service = getService ($C, $serviceName, $authURI);

$moby_response = $Service->execute (XMLinputlist => [
						     ["$articleName", "$input_xml\n", "clusters number", $cluster_number_xml, 'iterations number', $iteration_number_xml]
						     ]);

if ($_debug) {
	print STDERR "$serviceName results\n";
	print STDERR $moby_response;
	print STDERR "\n";
    }

if (hasFailed ($moby_response)) {
    print STDERR "service, $serviceName, has failed!\n";
    my $moby_error_message = getExceptionMessage ($moby_response);
    print STDERR "reason is the following,\n$moby_error_message\n";
    if (defined $output_dir) {
	print 1;
    }
    exit 1;
}

# Get the simples (collection of List_Text)

$input_xml_aref = parseResults ($moby_response, "List_Text");

print STDERR "Fourth step done!\n\n";

# runMultiMetaAlignmentGFF & runMultiMetaAlignment & runGFF2JPEG

print STDERR "Fifth step, foreach predicted gene cluster, make the multiple binding sites maps alignment...\n";

my $cluster_index = 1;
foreach my $gene_cluster_input_xml (@$input_xml_aref) {
    
    my $matscan_maps_xml = [];
    my $multimeta_map_xml;
    
    # Create a subdirectory which will contain cluster analysis results
    qx/mkdir $output_dir\/$cluster_index".cluster_results"/;
    
    # Filter MatScan results based on the Moby identifier
    # so prior, need to get the list of genes as an array
    
    my $gene_list_str_aref = parseTextContent ($gene_cluster_input_xml, "List_Text");
    my $gene_number = 0;
    
    # So far the list of genes in a formatted string
    my $gene_list_str = $gene_list_str_aref->[0];
    # Convert the list as an array of identifiers
    # and Filter empty lines
    my @cluster_genes = grep {$_ ne ""} split ("\n", $gene_list_str);
    
    if ($_debug) {
	print STDERR "cluster genes, " . join (', ', @cluster_genes) . ".\n";
    }
    
    # store the list of gene identifiers in a file
    open FILE, ">$output_dir/$cluster_index.cluster_results/$cluster_index.ids.lst" or die "can't open file, $output_dir/$cluster_index.cluster_results/$cluster_index.ids.lst!\n";
    print FILE join ("\n", @cluster_genes);
    print FILE "\n";
    close FILE;
    
    foreach my $gene_identifier (@cluster_genes) {
	
	if ($_debug) {
	    print STDERR "processing gene identifier, $gene_identifier.\n";
	}
	
	$gene_number++;
	my $map_xml = $matscan_results{$gene_identifier};
	push (@$matscan_maps_xml, $map_xml);
    }
    
    if ($gene_number > 1) {
	
	# runMultiMetaAlignmentGFF first - only when more than one gene in the cluster
	# Run this service just to save the results in GFF format and also for gff2ps visualization purposes
	
	print STDERR "Executing multiple meta-alignment...\n";
	
	$serviceName = "runMultiMetaAlignmentGFF";
	$authURI     = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
	$articleName = $parameters{$serviceName}->{articleName} || die "article name for $serviceName\n";
	
	$Service = getService ($C, $serviceName, $authURI);
	
	$moby_response = $Service->execute (XMLinputlist => [
							     ["$articleName", $matscan_maps_xml, 'alpha', $alpha_xml, 'lambda', $lambda_xml, 'mu', $mu_xml, 'gap penalty', $gamma_xml, 'NoN-colinear penalty', $non_collinearity_xml]
							     ]);
	
	if ($_debug) {
	    print STDERR "$serviceName result\n";
	    print STDERR $moby_response;
	    print STDERR "\n";
	}
	
	if (hasFailed ($moby_response)) {
	    print STDERR "service, $serviceName, has failed on cluster $cluster_index data!\n";
	    my $moby_error_message = getExceptionMessage ($moby_response);
	    print STDERR "reason is the following,\n$moby_error_message\n";
	}
	else {
	    my $multimeta_map_xml_aref = parseResults ($moby_response, "GFF");
	    $multimeta_map_xml = $multimeta_map_xml_aref->[0];
	    
	    if (defined $output_dir) {
		saveResults ($moby_response, "GFF", $cluster_index . ".MultiMeta", $output_dir . "/$cluster_index.cluster_results");
	    }
	    
	    # Then runMultiMetaAlignment
	    
	    $serviceName = "runMultiMetaAlignment";
	    $authURI     = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
	    $articleName = $parameters{$serviceName}->{articleName} || die "article name for $serviceName\n";
	    
	    $Service = getService ($C, $serviceName, $authURI);
	    
	    $moby_response = $Service->execute (XMLinputlist => [
								 ["$articleName", $matscan_maps_xml, 'alpha', $alpha_xml, 'lambda', $lambda_xml, 'mu', $mu_xml, 'gap penalty', $gamma_xml, 'NoN-colinear penalty', $non_collinearity_xml]
								 ]);
	    
	    if ($_debug) {
		print STDERR "$serviceName result\n";
		print STDERR $moby_response;
		print STDERR "\n";
	    }
	    
	    if (hasFailed ($moby_response)) {
		print STDERR "service, $serviceName, has failed!\n";
		my $moby_error_message = getExceptionMessage ($moby_response);
		print STDERR "reason is the following,\n$moby_error_message\n";
	    }
	    elsif (defined $output_dir) {
		saveResults ($moby_response, "Meta_Alignment_Text", $cluster_index . ".MultiMeta", $output_dir . "/$cluster_index.cluster_results");
	    }
	    
	    if ($_debug) {
		# print STDERR "input xml for next service:\n";
		# print STDERR join (', ', @$input_xml);
		# print STDERR ".\n";
	    }
	    
	    # runGFF2JPEG
	    
	    print STDERR "Executing gff2ps...\n";
	
	    $serviceName = "runGFF2JPEG";
	    $authURI     = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
	    $articleName = $parameters{$serviceName}->{articleName} || die "article name for $serviceName\n";
	    
	    $Service = getService ($C, $serviceName, $authURI);
	    
	    my $maps_input_xml = [];
	    push (@$maps_input_xml, $multimeta_map_xml);
	    push (@$maps_input_xml, @$matscan_maps_xml);
	    
	    if ($_debug) {
		print STDERR "runGFF2JPEG input, " . @$maps_input_xml . " maps\n";
	    }
	    
	    $moby_response = $Service->execute (XMLinputlist => [
								 ["$articleName", $maps_input_xml]
								 ]);
	    
	    if ($_debug) {
		print STDERR "$serviceName result\n";
		print STDERR $moby_response;
		print STDERR "\n";
	    }
	    
	    if (hasFailed ($moby_response)) {
		print STDERR "service, $serviceName, has failed!\n";
		my $moby_error_message = getExceptionMessage ($moby_response);
		print STDERR "reason is the following,\n$moby_error_message\n";
	    }
	    elsif (defined $output_dir) {
	    
		# Store the jpeg format image
		
		my $picture_b64_aref = parseTextContent ($moby_response, "b64_encoded_jpeg");
		
		if (! defined $picture_b64_aref) {
		    print STDERR "problem, can't parse anything from moby runGFF2JPEG response!\n";
		}
		else {
		    if ($_debug) {
			print STDERR "picture array has " . @$picture_b64_aref . " elements\n";
		    }
		    
		    # Store the image into a file
		    
		    my $picture_b64 = $picture_b64_aref->[0];
		    my $picture = decode_base64($picture_b64);
		    
		    open FILE, ">$output_dir/$cluster_index.cluster_results/$cluster_index.TFBSs_maps.jpg" or die "can't open file, $output_dir/$cluster_index.cluster_results/$cluster_index.TFBSs_maps.jpg!\n";
		    print FILE "$picture\n";
		    close FILE;
		}
	    }
	    
	    print STDERR "\n";
	}
    }
    
    $cluster_index++;
}

# Neighbour-joining (NJ) service execution disabled as it is not working !!

if (0) {
    
    print STDERR "Fifth step done!\n\n";
    
    # Request Inab registry for the two following ones
    $URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://www.inab.org/cgi-bin/MOBY-Central.pl';
    $URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://www.inab.org/MOBY/Central';
    $PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';
    
    $C = MOBY::Client::Central->new(
				    Registries => {mobycentral => {URL => $URL,URI => $URI}}
				    );
    
    if (defined $C) {
	if ($_debug) {
	    print STDERR "Instanciation done.\n";
	    print STDERR "TESTING MOBY CLIENT with\n\tURL: $URL\n\tURI: $URI\n\tProxy: $PROXY\n\n";
	}
    }
    else {
	print STDERR "Error, Central could not be instanciated!\n";
	if (defined $output_dir) {
	    print 1;
	}
	exit 1;
    }
    
    # runHierarchicalClustering
    
    print STDERR "Sixth step, gene clustering using a $method neighbour joining clustering algorithm...\n";
    
    if ($_debug) {
	print STDERR "\nExecuting runHierarchicalClustering...\n\n";
    }
    
    $serviceName   = "runHierarchicalClustering";
    $authURI       = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
    $articleName   = $parameters{$serviceName}->{articleName} || "";
    
    $Service = getService ($C, $serviceName, $authURI);
    
    $moby_response = $Service->execute (XMLinputlist => [
							 ["$articleName", "$input_xml\n", "method", $method_xml]
							 ]);
    if (! defined $moby_response) {
	# moby service is not working, so use GEPAS instead !!!
	
	# GEPAS doesn't return any results, no idea why because it works fine with GeneID web server !!! ????
	
	print STDERR "using GEPAS...\n";
	
	print STDERR "method, $method\n";
	
	if ($method =~ /nearest/i) {
	    $method = "single";
	}
	else {
	    $method = "complete";
	}
	
	# Assume the output flag was on, otherwise need to create a temporary file with the data which would be a pain !!!
	my $curdir = qx/pwd/;
	chomp $curdir;
	
	my $matrix_filename = "score_matrix.txt";
	
	if (! -f "$output_dir/$matrix_filename") {
	    print STDERR "problem, score matrix file not found, $output_dir/$matrix_filename!!\n";
	}
	
	my $cluster_url = "http://gepas.bioinfo.cipf.es/cgi-bin/cluster";
	my $agent_diag   = LWP::UserAgent->new(timeout => 30000);
	
	# Cluster CGI Call
	
	my $request_diag = POST($cluster_url,
				Content_Type => 'form-data',
				Content      => [
						 file     => ["$output_dir/$matrix_filename"],
						 method   => "$method",
						 ]
				);
	my $result_diag = $agent_diag->request($request_diag);
	
	# print STDERR "Dumping result_diag, " . Dumper ($result_diag) . "\n";
	
	my $html_results = $result_diag->content;
	
	# Get the data from GEPAS server !
	my $newick_tree = "";
	if (defined $output_dir) {
	    open FILE, ">$output_dir/newick.txt" or die "can't open in write access newick text file, $output_dir/newick.txt!\n";
	    print FILE $newick_tree;
	    close FILE;
	}
	
	# Call TreeView
	
	# ...
	
    }
    else {
	
	# Fine so carry on using BioMOBY infrastructure...
	
	if ($_debug) {
	    print STDERR "\n$serviceName results\n";
	    print STDERR $moby_response;
	    print STDERR "\n";
	}
	
	# Check if we got results
	my $newick_trees = parseTextContent ($moby_response, "Newick_Text");
	my $newick_tree  = $newick_trees->[0];
	
	if ($newick_tree eq "") {
	    print STDERR "Hierachical clustering failed!\n";
	    if (defined $output_dir) {
		print 1;
	    }
	    exit 1;
	}
	
	$input_xml_aref = parseResults ($moby_response, "HierarchicalClustering");
	if (defined $output_dir) {
	    saveResults ($moby_response, "Newick_Text", "newick", $output_dir);
	}
	
	if ($_debug) {
	    print STDERR "input xml for next service:\n";
	    print STDERR join (', ', @$input_xml_aref);
	    print STDERR ".\n";
	}
	
	# convert the input xml into a scalar
	$input_xml =  $input_xml_aref->[0];
	
	print STDERR "Clustering done!\n\n";
	
	# plotClustersTree
	
	print STDERR "Generating a picture of the clustering tree...\n";
	
	$serviceName   = "plotClustersTree";
	$authURI       = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
	$articleName   = $parameters{$serviceName}->{articleName} || "";
	
	$Service = getService ($C, $serviceName, $authURI);
	
	$moby_response = $Service->execute (XMLinputlist => [
							     ["$articleName", $input_xml]
							     ]);
	
	if ($_debug) {
	    print STDERR "$serviceName results\n";
	    print STDERR $moby_response;
	    print STDERR "\n";
	}
	
	my $picture_b64_aref = parseTextContent ($moby_response, "Image_Encoded");
	my $picture_b64 = $picture_b64_aref->[0];
	
	# Convert in to a picture and store into a file
	
	my $picture = decode_base64($picture_b64);
	
	print STDERR "Picture done\n\n";
	print STDERR "Workflow has terminated.\n";
	if (! defined $output_dir) {
	    print $picture;
	}
	else {
	    open FILE, ">$output_dir/clustering_tree.png" or die "can't open in write access file, '$output_dir/clustering_tree.png'!\n";
	    print FILE $picture;
	    close FILE;
	}
    }
    
} # End disabling execution of NJ service

my $t2 = Benchmark->new ();
print STDERR "\nTotal : ", timestr (timediff ($t2, $t1)), "\n";

exit 0;

##
# The end
##

sub setConfigurationData {
    my ($config_file) = @_;
    my %parameters;
    
    open CONFIG, "<$config_file" or die "can't open config file, $config_file\n";
    while (<CONFIG>) {
	my $line = $_;
	
	next if $line =~ /^#/;
	
	if ($line =~ /^>(.+)/) {
	    $serviceName = $1;
	    chomp $serviceName;
	    if (! defined $serviceName) {
		print STDERR "pb when parsing config file - can't get the service name\n";
		print STDERR "line , $line\n";
	    }
	    else {
	        if ($_debug) {
  		  print STDERR "parsed following article name, $serviceName\n";
                }
	    }
	}
	elsif ($line =~ /(.+)=(.+)/) {
	    
	    if ($_debug) {
  	      print STDERR "parsing a parameter...\n";
            }
	    
	    if (! defined $serviceName) {
		print STDERR "pb - no service name defined!\n";
	    }
	    
	    my $parameter_name  = $1;
	    my $parameter_value = $2;
	    
	    if ($_debug) {
              print STDERR "parsed name, value, $parameter_name, $parameter_value.\n";
            }
	    
	    my $parameters_tmp;
	    if (defined $parameters{$serviceName}) {
		$parameters_tmp = $parameters{$serviceName};
		$parameters_tmp->{$parameter_name} = $parameter_value;
	    }
	    else {
		$parameters_tmp = {
		    $parameter_name => $parameter_value,
		};
		$parameters{$serviceName} = $parameters_tmp;
	    }
	}
    }
    close CONFIG;

    return %parameters;
}

sub getService {
    my ($C, $serviceName, $authURI) = @_;
    my ($service_instances, $reg) = $C->findService (
						     serviceName  => $serviceName,
						     authURI      => $authURI
						     );
    
    if ($_debug) {
	print STDERR "Service instances references: " . @$service_instances . "\n";
    }
    
    if (@$service_instances == 0) {
	print STDERR "Error, can't find any service called $serviceName!\n";
	return 0;
    }
    
    my $service_instance = $service_instances->[0];
    
    my $wsdl = $C->retrieveService($service_instance);
    
    if ($_debug) {
	print STDERR "wsdl: $wsdl\n";
    }
    
    if (!$wsdl || ($wsdl !~ /\<definitions/)){
	if ($_debug) {
	    print STDERR "test \t\t[FAIL]\tWSDL was not retrieved\n\n";
	}
    }
    
    my $Service = MOBY::Client::Service->new(service => $wsdl);
    
    return $Service;
}


# parse the MobyData bloc to get the collection or simple bloc
# return an array of objects
sub parseResults {
    my ($XML, $object_type) = @_;
    my $inputs_xml = [];
    
    my $parser = XML::LibXML->new();
    my $doc    = $parser->parse_string( $XML );
    $XML       = $doc->getDocumentElement();
    my $elements = $XML->getElementsByTagName( "moby:$object_type" );
    
    my $size = $elements->size();
    
    if ($size == 0) {
	$elements = $XML->getElementsByTagName( "object_type" );
	$size = $elements->size();
	if ($size == 0) {
	    print STDERR "Error, can't parse the moby output from the moby XML...\n";
	}
    }
    
    my $i = 0;
    while ($i < $size) {
	my $element   = $elements->get_node ($i);
	my $input_xml = $element->toString();
	push (@$inputs_xml, $input_xml);
	
	$i++;
    }
    
    return $inputs_xml;
}

sub parseTextContent {
    my ($XML, $object_type) = @_;
    my $inputs_text = [];
    
    # If the XML is not embbeded in a MOBY element, add it because otherwise we won't be able to instanciate an XML document
    # I think it is compulsory to define a namespace attribute (xmlns)
    
    if (! ($XML =~ /MOBY/)) {
	$XML = "<?xml version='1.0' encoding='UTF-8'?><moby:MOBY xmlns:moby='http://www.biomoby.org/moby' xmlns='http://www.biomoby.org/moby'>\n" . $XML;
	$XML .= "</moby:MOBY>\n";
	
	if ($_debug) {
	    print STDERR "List_Text element, $XML\n";
	}
    }
    
    my $parser = XML::LibXML->new();
    my $doc    = $parser->parse_string( $XML );
    $XML       = $doc->getDocumentElement();
    my $elements = $XML->getElementsByTagName ( "moby:$object_type" );
    my $size   = $elements->size();
    
    if ($size == 0) {
      $elements = $XML->getElementsByTagName ("$object_type");
      $size = $elements->size();
      if ($size == 0) {
        print STDERR "Error, can't parse the moby output from the moby XML...\n";
      }
    }

    my $i = 0;
    while ($i < $size) {
	my $element    = $elements->get_node ($i);
	my $input_text = $element->textContent();
	push (@$inputs_text, $input_text);
	$i++;
    }

    return $inputs_text;
}

sub saveResults {
  my ($XML, $object_type, $softwareName, $output_dir) = @_;
  
  if ($_debug) {
    print STDERR "saving results for tool, $softwareName...\n";
  }
  
  my %inputs_text;
    
  my $parser = XML::LibXML->new();
  my $doc    = $parser->parse_string( $XML );
  $XML       = $doc->getDocumentElement();
  my $elements = $XML->getElementsByTagName ( "moby:$object_type" );
  my $size   = $elements->size();
    
  if ($size == 0) {
    $elements = $XML->getElementsByTagName ("$object_type");
    $size = $elements->size();
    if ($size == 0) {
      print STDERR "Error, can't parse the moby output from the moby XML...\n";
    }
  }

  if ($_debug) {
    print STDERR "found $size elements\n";
  }

  # parsing
  
  my $i = 0;
  while ($i < $size) {
    my $element = $elements->[$i];
    my $id = $element->getAttribute ("id");
    
    # if more than one simple object, create an identifier index
    if ((! defined $id) || $id eq "") {
      if ($size > 1) {
        $id = $i + 1;
      }
      else {
        $id = "";
      }
    }
    my $input_text = $element->textContent();
    if ($object_type eq "b64_Encoded_PNG") {
      $input_text = decode_base64($input_text);
    }
    
    $inputs_text{$id} = $input_text;
    $i++;
  }

  my $suffix = getFileNameSuffix ($object_type);

  # Save on the local disk the data
  my @ids = keys (%inputs_text);
  if (@ids == 1) {
    # simple => just one file
    my $filename = $softwareName . "." . $suffix;
    my $id = $ids[0];
    if ($id ne "") {
      $filename = $id . "." . $filename;
    }
    
    my $data = $inputs_text{$id};
    # clean data
    
    $data =~ s/^\n+//;
    $data =~ s/\n+$/\n/;
    
    open FILE, ">$output_dir/$filename" or die "can' open output file, $output_dir/$filename";
    print FILE "$data";
    close FILE;
  }
  else {
    # collection of simples => a sub directory
    my $subdir = "$softwareName";
    
    if (! -d "$output_dir/$subdir") {
      qx/mkdir $output_dir\/$subdir/;
    }
    
    foreach my $id (@ids) {
      my $filename = $id . "." . $softwareName . "." . $suffix;
      my $data = $inputs_text{$id};
      # clean data
      
      $data =~ s/^\n+//;
      $data =~ s/\n+$/\n/;
      
      open FILE, ">$output_dir/$subdir/$filename" or die "can' open output file, $output_dir/$subdir/$filename";
      print FILE "$data";
      close FILE;
    }
  }
}

# Set up a file name suffix giving the input object type
# e.g. if the object is GFF, the file suffix will be '.gff'

sub getFileNameSuffix {
  my ($object_type) = @_;
  my $suffix = "txt";
  
  if (lc ($object_type) eq "gff") {
    return "gff";
  }
  elsif (lc ($object_type) eq "text-formatted" || lc ($object_type) eq "text_formatted" || lc ($object_type) eq "microarraydata_text" || lc ($object_type) eq "newick_text") {
    return "txt";
  }
  elsif (lc ($object_type) eq "b64_encoded_png") {
    return "png";
  }
  elsif (lc ($object_type) eq "fasta") {
      return "fa";
  }
  elsif (lc ($object_type) eq "list_text" || lc ($object_type) eq "list-text") {
      return "lst";
  }
  
  if ($_debug) {
      print STDERR "Suffix for object_type, $object_type unknown, return txt!\n"; 
  }
  
  return $suffix;
}


# parse the moby exception status
sub hasFailed {
    my ($moby_response) = @_;
    
    # convert 
    # "&lt;Notes&gt;Service execution succeeded&lt;/Notes&gt;"
    # into
    # "<Notes>Service execution succeeded</Notes>"
    
    $moby_response = HTML::Entities::decode( $moby_response );
    
    if (($moby_response =~ /service execution succeeded/i) || ($moby_response =~ /exceptionCode>700/)) {
	# succeedded
	return 0;
    }
    else {
	# failed
	return 1;
    }
}

sub getExceptionMessage {
    my ($moby_response) = @_;
    
    # convert string like the following one
    # "&lt;Notes&gt;Service execution succeeded&lt;/Notes&gt;"
    # into
    # "<Notes>Service execution succeeded</Notes>"
    
    $moby_response = HTML::Entities::decode( $moby_response );
    my $exceptionMessages = parseTextContent ($moby_response, "exceptionMessage");
    
    my $exceptionMessage = $exceptionMessages->[0];
    if (defined $exceptionMessage) {
	# Clean the message
	while ($exceptionMessage =~ s/^\n//) {
	}
	while (chomp $exceptionMessage) { # just chomp... 
	}
    }
    else {
	print STDERR "exception message could not be parsed!\n";
	$exceptionMessage = "";
    }
    
    return $exceptionMessage;
}

sub setMatScan_hash {
  my ($matscan_results_xml_aref) = @_;
  my %matscan_results;
  
  foreach my $matscan_results_xml (@$matscan_results_xml_aref) {
    
    # print STDERR "parsed matscan results xml, $matscan_results_xml\n"; 
    
    $matscan_results_xml =~ /id="([^\"]+)"/;
    my $gene_identifier = $1;
    
    if ($_debug) {
      print STDERR "parsed gene identifier, $gene_identifier\n";
    }
    
    $matscan_results{$gene_identifier} = $matscan_results_xml;
  }
  
  return %matscan_results;
}

