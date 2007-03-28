#!/usr/local/bin/perl -w

##################################################################
#
#  INB Genomics Node Moby Services testing client
#
##################################################################

use strict;

# be prepare for command-line options/arguments
use Getopt::Std;
use Benchmark;
use Data::Dumper;
# BioMoby and SOAP libraries

use MOBY::Client::Central;
use MOBY::Client::Service;
# use SOAP::Lite + 'trace';

use File::Temp qw/tempfile/;

sub help {
return <<"END_HELP";
Description: Execute Moby services available from genome.imim.es
Usage:

    testServices.pl [-h] -x {Moby Central}
	-h help
	-x MOBY Central: Chirimoyo, Mobydev, Inab or BioMoby
		<1> or Chirimoyo
		<2> or Mobydev
		<3> or Inab
		<4> or BioMoby

Examples using some combinations:
	testServices -x 1

END_HELP

}

BEGIN {

	# Determines the options with values from program
	use vars qw/$opt_h $opt_x $URL $URI $PROXY $AUTH $_debug/;

	# these are switches taking an argument (a value)
	my $switches = 'hx';

	# Get the switches
	getopt($switches);

	# If the user does not write nothing, skip to help
	if (defined($opt_h) || !defined($opt_x)){
		print STDERR help;
		exit 0;
	}

}

my $t1 = Benchmark->new ();
$_debug = 0;
$AUTH = "genome.imim.es";

##################################################################
#
# Setup sequence input data
#
##################################################################

my $rootdir                   = $ENV{HOME} . "/cvs/GRIB/projects/moby/scripts/testServices";
my $input_data_dir            = "$rootdir/inputData";
my $control_data_dir          = "$rootdir/outputControl";

my $nucleotide_sequence_xml_file = "Hsap_BTK.msk.xml";
my $tblastx_output_xml_file   = "Hsap_BTK.tbx.xml";
my $GeneIDGFF_xml_file        = "Hsap_BTK.msk.GeneIDGFF.xml";
my $geneIds_lst_xml_file      = "geneIds.lst.xml";
my $gostat_regulated_xml_file = "mut1_downreg.fbgn.xml";
my $gostat_allArray_xml_file  = "allArray.fbgn.xml";
my $ENSG00000197785_upstream_sequence_xml_file = "ENSG00000197785.xml";
my $ENSG00000197785_matscan_xml_file           = "ENSG00000197785.runMatScanGFF.xml";
my $ENSPTRG00000000031_matscan_xml_file        = "ENSPTRG00000000031.runMatScanGFF.xml";
my $ENSG00000197785_fasta_xml_file             = "ENSG00000197785.fa.xml";
my $Dmel_MultiMeta_xml_file   = "Dmel.MultiMetaOutput.xml";
my $Lep_Gasdermin_score_matrix_xml_file        = "Lep_Gasdermin_score_matrix.xml";

# Check that the files exist !!!

if ((not -f "$input_data_dir/$nucleotide_sequence_xml_file") || (not -f "$input_data_dir/$tblastx_output_xml_file") || (not -f "$input_data_dir/$geneIds_lst_xml_file") || (not -f "$input_data_dir/$GeneIDGFF_xml_file") || (not -f "$input_data_dir/$gostat_regulated_xml_file") || (not -f "$input_data_dir/$gostat_allArray_xml_file") || (not -f "$input_data_dir/$ENSG00000197785_upstream_sequence_xml_file") || (not -f "$input_data_dir/$ENSG00000197785_matscan_xml_file") || (not -f "$input_data_dir/$ENSPTRG00000000031_matscan_xml_file") || (not -f "$input_data_dir/$ENSG00000197785_fasta_xml_file") || (not -f "$input_data_dir/$Dmel_MultiMeta_xml_file") || (not -f "$input_data_dir/$Lep_Gasdermin_score_matrix_xml_file")) {
    print STDERR "Error, can't find one of the input files in directory, $input_data_dir!\n";
    exit 1;
}

my $tblastx_output_xml   = qx/cat $input_data_dir\/$tblastx_output_xml_file/;
my $nucleotide_sequence_xml = qx/cat $input_data_dir\/$nucleotide_sequence_xml_file/;
my $GeneIDGFF_xml        = qx/cat $input_data_dir\/$GeneIDGFF_xml_file/;
my $geneIds_lst_xml      = qx/cat $input_data_dir\/$geneIds_lst_xml_file/;
my $gostat_regulated_xml = qx/cat $input_data_dir\/$gostat_regulated_xml_file/;
my $gostat_allArray_xml  = qx/cat $input_data_dir\/$gostat_allArray_xml_file/;
my $ENSG00000197785_upstream_sequence_xml = qx/cat $input_data_dir\/$ENSG00000197785_upstream_sequence_xml_file/;
my $ENSG00000197785_matscan_xml           = qx/cat $input_data_dir\/$ENSG00000197785_matscan_xml_file/;
my $ENSPTRG00000000031_matscan_xml        = qx/cat $input_data_dir\/$ENSPTRG00000000031_matscan_xml_file/;
my $ENSG00000197785_fasta_xml             = qx/cat $input_data_dir\/$ENSG00000197785_fasta_xml_file/;
my $Dmel_MultiMeta_xml   = qx/cat $input_data_dir\/$Dmel_MultiMeta_xml_file/;
my $Lep_Gasdermin_score_matrix_xml        = qx/cat $input_data_dir\/$Lep_Gasdermin_score_matrix_xml_file/;

my $runGeneID_control_file                  = "Hsap_BTK.msk.runGeneID.control";
my $runGeneIDGFF_control_file               = "Hsap_BTK.msk.runGeneIDGFF.control";
my $runSGP2GFF_control_file                 = "Hsap_BTK.msk.runSGP2GFF.control";
my $translateGeneIDGFFPredictions_control_file = "Hsap_BTK.msk.GeneIDGFF.translateGeneIDGFFPredictions.control";
my $getUpstreamSeqFromEnsembl_control_file  = "geneIds.lst.getUpstreamSeqFromEnsembl.control";
my $runGOstat_control_file                  = "mut1_downreg.fbgn.runGOstat.control";
my $fromGenericSequenceToFASTA_control_file = "Hsap_BTK.msk.fromGenericSequenceToFASTA.control";
my $fromGenericSequenceCollectionToFASTA_control_file = "Hsap_BTK.msk.fromGenericSequenceCollectionToFASTA.control";
my $runMatScanGFF_control_file              = "ENSG00000197785.runMatScanGFF.control";
my $runMatScanGFFCollection_control_file    = "ENSG00000197785.runMatScanGFFCollection.control";
my $runMetaAlignment_control_file           = "ENSG00000197785.runMetaAlignment.control";
my $runMetaAlignmentGFF_control_file        = "ENSG00000197785.runMetaAlignmentGFF.control";
my $fromFASTAToDNASequenceCollection_control_file    = "ENSG00000197785.fromFASTAToDNASequenceCollection.control";
my $fromMetaAlignmentsToTextScoreMatrix_control_file = "mut1_downreg.fbgn.ScoresMatrix.control";
my $runRepeatMasker_control_file = "ENSG00000197785.runRepeatMasker.control";
my $runDust_control_file = "ENSG00000197785.runDust.control";
my $runSOTAClustering_control_file = "Lep_Gasdermin.runSOTAClustering.control";

##################################################################
#
# Setup Moby configuration parameters
#
##################################################################

if (defined($opt_x)) {

    # Assign the MOBY Server and MOBY URI
    if (($opt_x == 1) || ($opt_x eq 'Chirimoyo')) {

	$URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';
	$URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://chirimoyo.ac.uma.es/MOBY/Central';
	$PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

    }elsif (($opt_x == 2) || ($opt_x eq 'Mobydev')) {

	$URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://moby-dev.inab.org/cgi-bin/MOBY-Central.pl';
	$URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://moby-dev.inab.org/MOBY/Central';
	$PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

    }elsif (($opt_x == 3) || ($opt_x eq 'Inab')) {

	$URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://www.inab.org/cgi-bin/MOBY-Central.pl';
	$URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://www.inab.org/MOBY/Central';
	$PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

    }elsif (($opt_x == 4) || ($opt_x eq 'BioMoby')) {

	$URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://mobycentral.icapture.ubc.ca/cgi-bin/MOBY05/mobycentral.pl';
	$URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://mobycentral.icapture.ubc.ca/MOBY/Central';
	$PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

    }else {
	print STDERR help;
	exit 0;
    }

}else {
    print STDERR help;
    exit 0;
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
    exit 0;
}

my $service;

# Execute runGeneID Web service

###################################
#
# Disabled !!
#
###################################

# print STDERR "testing runGeneID...\n\n";

# my $service = MobyServiceInstantiation ($C, "runGeneID", $AUTH);
# if (defined $service) {
#    my $result = $service->execute(
#				   XMLinputlist => [
#						    ['sequences', $nucleotide_sequence_xml]
#						    ]
#				   ) ;
#
#    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);
#
#    print $results_fh "$result\n";
#
#    my @diff_results = qx/diff $control_data_dir\/$runGeneID_control_file $results_file/;
#
#    # The execution date report is obviously different so don't check this !!
#
#    if ((@diff_results > 0) && (! ($diff_results[1] =~ /date/))) {
#
#	# Don't check either the temporary file name which is also different !
#
#	if (! ($diff_results[5] =~ /geneid\d+monstre1.imim.es.fasta/)) {
#
#	    print STDERR "runGeneID failed!\n";
#	    print STDERR "diff_results: @diff_results\n";
#
#	    close $results_fh;
#	    unlink $results_file;
#
#	}
#	else {
#
#	    print STDERR "runGeneID okay...\n";
#
#	    close $results_fh;
#	    unlink $results_file;
#	}
#    }
#    else {
#
#	print STDERR "runGeneID okay...\n";
#
#	close $results_fh;
#	unlink $results_file;
#    }
#}

# Execute runGeneIDGFF Web service

print "\ntesting runGeneIDGFF...\n\n";

$service = MobyServiceInstantiation ($C, "runGeneIDGFF", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['sequences', $nucleotide_sequence_xml]
						    ]
				   );

    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);

    print $results_fh "$result\n";

    my @diff_results = qx/diff $control_data_dir\/$runGeneIDGFF_control_file $results_file/;

    if ((@diff_results > 0) && (! ($diff_results[1] =~ /date/))) {
	print STDERR "runGeneIDGFF service failed!\n";
	print STDERR "diff_results: @diff_results\n";

	close $results_fh;
	unlink $results_file;

    }
    else {

	print "runGeneIDGFF okay...\n";

	close $results_fh;
	unlink $results_file;
    }
}

# Execute translateGeneIDGFFPredictions Web service

print  "\ntesting translateGeneIDGFFPredictions...\n\n";

$service = MobyServiceInstantiation ($C, "translateGeneIDGFFPredictions", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['sequence', $nucleotide_sequence_xml, 'geneid_predictions', $GeneIDGFF_xml]
						    ]
				   );

    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);

    print $results_fh "$result\n";

    my @diff_results = qx/diff $control_data_dir\/$translateGeneIDGFFPredictions_control_file $results_file/;

    if (@diff_results > 0) {
	print STDERR "translateGeneIDGFFPredictions service failed!\n";
	print STDERR "diff_results: @diff_results\n";

	close $results_fh;
	unlink $results_file;

    }
    else {

	print  "translateGeneIDGFFPredictions okay...\n";

	close $results_fh;
	unlink $results_file;
    }
}

# Execute getUpstreamSeqFromEnsembl Web service

print  "\ntesting getUpstreamSeqFromEnsembl...\n\n";

$service = MobyServiceInstantiation ($C, "getUpstreamSeqFromEnsembl", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['genes', $geneIds_lst_xml]
						    ]
				   );

    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);

    print $results_fh "$result\n";

    my @diff_results = qx/diff $control_data_dir\/$getUpstreamSeqFromEnsembl_control_file $results_file/;

    if (@diff_results > 0) {
	print STDERR "getUpstreamSeqFromEnsembl service failed!\n";
	print STDERR "diff_results: @diff_results\n";

	close $results_fh;
	unlink $results_file;

    }
    else {

	print  "getUpstreamSeqFromEnsembl okay...\n";

	close $results_fh;
	unlink $results_file;
    }
}

# Execute fromGenericSequenceToFASTA Web service

print  "\ntesting fromGenericSequenceToFASTA...\n\n";

$service = MobyServiceInstantiation ($C, "fromGenericSequenceToFASTA", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['sequences', $nucleotide_sequence_xml]
						    ]
				   );

    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);

    print $results_fh "$result\n";

    my @diff_results = qx/diff $control_data_dir\/$fromGenericSequenceToFASTA_control_file $results_file/;

    if ((@diff_results > 0) && (! ($diff_results[1] =~ /date/))) {
	print STDERR "fromGenericSequenceToFASTA service failed!\n";
	print STDERR "diff_results: @diff_results\n";

	close $results_fh;
	unlink $results_file;

    }
    else {

	print  "fromGenericSequenceToFASTA okay...\n";

	close $results_fh;
	unlink $results_file;
    }
}

# Execute fromGenericSequenceCollectionToFASTA Web service

print  "\ntesting fromGenericSequenceCollectionToFASTA...\n\n";

$service = MobyServiceInstantiation ($C, "fromGenericSequenceCollectionToFASTA", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['sequences', [$nucleotide_sequence_xml]]
						    ]
				   );

    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);

    print $results_fh "$result\n";

    my @diff_results = qx/diff $control_data_dir\/$fromGenericSequenceCollectionToFASTA_control_file $results_file/;

    if ((@diff_results > 0) && (! ($diff_results[1] =~ /date/))) {
	print STDERR "fromGenericSequenceCollectionToFASTA service failed!\n";
	print STDERR "diff_results: @diff_results\n";

	close $results_fh;
	unlink $results_file;

    }
    else {

	print  "fromGenericSequenceCollectionToFASTA okay...\n";

	close $results_fh;
	unlink $results_file;
    }
}

# Execute runMatScanGFF Web service

print  "\ntesting runMatScanGFF...\n\n";

$service = MobyServiceInstantiation ($C, "runMatScanGFF", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['upstream_sequences', $ENSG00000197785_upstream_sequence_xml]
						    ]
				   );

    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);

    print $results_fh "$result\n";

    my @diff_results = qx/diff $control_data_dir\/$runMatScanGFF_control_file $results_file/;

    if ((@diff_results > 0) && (! ($diff_results[1] =~ /date/))) {
	print STDERR "runMatScanGFF service failed!\n";
	print STDERR "diff_results: @diff_results\n";

	close $results_fh;
	unlink $results_file;

    }
    else {

	print  "runMatScanGFF okay...\n";

	close $results_fh;
	unlink $results_file;
    }
}

# Execute runMatScanGFFCollection Web service

print  "\ntesting runMatScanGFFCollection...\n\n";

$service = MobyServiceInstantiation ($C, "runMatScanGFFCollection", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['upstream_sequences', [$ENSG00000197785_upstream_sequence_xml]]
						    ]
				   );

    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);

    print $results_fh "$result\n";

    my @diff_results = qx/diff $control_data_dir\/$runMatScanGFFCollection_control_file $results_file/;

    if ((@diff_results > 0) && (! ($diff_results[1] =~ /date/))) {
	print STDERR "runMatScanGFFCollection service failed!\n";
	print STDERR "diff_results: @diff_results\n";

	close $results_fh;
	unlink $results_file;

    }
    else {

	print  "runMatScanGFFCollection okay...\n";

	close $results_fh;
	unlink $results_file;
    }
}

# Execute runMetaAlignment Web service

print  "\ntesting runMetaAlignment...\n\n";

$service = MobyServiceInstantiation ($C, "runMetaAlignment", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['map1', "$ENSG00000197785_matscan_xml", 'map2', "$ENSPTRG00000000031_matscan_xml"]
						    ]
				   );

    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);

    print $results_fh "$result\n";

    my @diff_results = qx/diff $control_data_dir\/$runMetaAlignment_control_file $results_file/;

    if ((@diff_results > 0) && (! ($diff_results[1] =~ /date/))) {
	print STDERR "runMetaAlignment service failed!\n";
	print STDERR "diff_results: @diff_results\n";

	close $results_fh;
	unlink $results_file;

    }
    else {

	print  "runMetaAlignment okay...\n";

	close $results_fh;
	unlink $results_file;
    }

}

# Execute runMetaAlignmentGFF Web service

print  "\ntesting runMetaAlignmentGFF...\n\n";

$service = MobyServiceInstantiation ($C, "runMetaAlignmentGFF", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['map1', "$ENSG00000197785_matscan_xml", 'map2', "$ENSPTRG00000000031_matscan_xml"]
						    ]
				   );

    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);

    print $results_fh "$result\n";

    my @diff_results = qx/diff $control_data_dir\/$runMetaAlignmentGFF_control_file $results_file/;

    if ((@diff_results > 0) && (! ($diff_results[1] =~ /date/))) {
	print STDERR "runMetaAlignmentGFF service failed!\n";
	print STDERR "diff_results: @diff_results\n";

	close $results_fh;
	unlink $results_file;

    }
    else {

	print  "runMetaAlignmentGFF okay...\n";

	close $results_fh;
	unlink $results_file;
    }

}

# Execute fromFASTAToDNASequenceCollection Web service

print  "\ntesting fromFASTAToDNASequenceCollection...\n\n";

$service = MobyServiceInstantiation ($C, "fromFASTAToDNASequenceCollection", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['sequences', $ENSG00000197785_fasta_xml]
						    ]
				   );

    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);

    print $results_fh "$result\n";

    my @diff_results = qx/diff $control_data_dir\/$fromFASTAToDNASequenceCollection_control_file $results_file/;

    if ((@diff_results > 0) && (! ($diff_results[1] =~ /date/))) {
	print STDERR "fromFASTAToDNASequenceCollection service failed!\n";
	print STDERR "diff_results: @diff_results\n";

	close $results_fh;
	unlink $results_file;

    }
    else {

	print  "fromFASTAToDNASequenceCollection okay...\n";

	close $results_fh;
	unlink $results_file;
    }
}

# Execute fromMetaAlignmentsToTextScoreMatrix Web service

print  "\ntesting fromMetaAlignmentsToTextScoreMatrix...\n\n";

# Construct the meta output moby objects array from the concatenated string

my @meta_output_objects = split ("\n\n", $Dmel_MultiMeta_xml);

$service = MobyServiceInstantiation ($C, "fromMetaAlignmentsToTextScoreMatrix", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['meta_predictions', \@meta_output_objects]
						    ]
				   );

    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);

    print $results_fh "$result\n";

    my @diff_results = qx/diff $control_data_dir\/$fromMetaAlignmentsToTextScoreMatrix_control_file $results_file/;

    if ((@diff_results > 0) && (! ($diff_results[1] =~ /date/))) {
	print STDERR "fromMetaAlignmentsToTextScoreMatrix service failed!\n";
	print STDERR "diff_results: @diff_results\n";

	close $results_fh;
	unlink $results_file;

    }
    else {

	print  "fromMetaAlignmentsToTextScoreMatrix okay...\n";

	close $results_fh;
	unlink $results_file;
    }
}

# Execute runSGP2GFF Web service

print  "\ntesting runSGP2GFF...\n\n";

$service = MobyServiceInstantiation ($C, "runSGP2GFF", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['sequence', $nucleotide_sequence_xml, 'tblastx_report', $tblastx_output_xml]
						    ]
				   );

    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);

    print $results_fh "$result\n";

    my @diff_results = qx/diff $control_data_dir\/$runSGP2GFF_control_file $results_file/;

    if ((@diff_results > 0) && (! ($diff_results[1] =~ /date/))) {
	print STDERR "runSGP2GFF service failed!\n";
	print STDERR "diff_results: @diff_results\n";

	close $results_fh;
	unlink $results_file;

    }
    else {

	print  "runSGP2GFF okay...\n";

	close $results_fh;
	unlink $results_file;
    }

}

# Execute runRepeatMasker Web service

print  "\ntesting runRepeatMasker...\n\n";

$service = MobyServiceInstantiation ($C, "runRepeatMasker", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['sequence', $ENSG00000197785_upstream_sequence_xml]
						    ]
				   );

    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);

    print $results_fh "$result\n";

    my @diff_results = qx/diff $control_data_dir\/$runRepeatMasker_control_file $results_file/;

    if ((@diff_results > 0) && (! ($diff_results[1] =~ /date/))) {
	print STDERR "runRepeatMasker service failed!\n";
	print STDERR "diff_results: @diff_results\n";

	close $results_fh;
	unlink $results_file;

    }
    else {

	print  "runRepeatMasker okay...\n";

	close $results_fh;
	unlink $results_file;
    }

}

print  "\ntesting runDust...\n\n";

$service = MobyServiceInstantiation ($C, "runDust", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['sequence', $ENSG00000197785_upstream_sequence_xml]
						    ]
				   );

    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);

    print $results_fh "$result\n";

    my @diff_results = qx/diff $control_data_dir\/$runDust_control_file $results_file/;

    if ((@diff_results > 0) && (! ($diff_results[1] =~ /date/))) {
	print STDERR "runDust service failed!\n";
	print STDERR "diff_results: @diff_results\n";

	close $results_fh;
	unlink $results_file;

    }
    else {

	print  "runDust okay...\n";

	close $results_fh;
	unlink $results_file;
    }

}

print  "\ntesting runSOTAClustering...\n\n";

$service = MobyServiceInstantiation ($C, "runSOTAClustering", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['sequence', $Lep_Gasdermin_score_matrix_xml]
						    ]
				   );

    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);

    print $results_fh "$result\n";

    my @diff_results = qx/diff $control_data_dir\/$runSOTAClustering_control_file $results_file/;

    if ((@diff_results > 0) && (! ($diff_results[1] =~ /date/))) {
	print STDERR "runSOTAClustering service failed!\n";
	print STDERR "diff_results: @diff_results\n";

	close $results_fh;
	unlink $results_file;

    }
    else {

	print  "runSOTAClustering okay...\n";

	close $results_fh;
	unlink $results_file;
    }

}



if (($opt_x == 1) || ($opt_x eq 'Chirimoyo')) {

    # Execute runGOstat Web service

    print  "\ntesting runGOstat...\n\n";

    $service = MobyServiceInstantiation ($C, "runGOstat", $AUTH);
    if (defined $service) {
	my $result = $service->execute(
				       XMLinputlist => [
							['regulated_genes', $gostat_regulated_xml, 'reference_genes', $gostat_allArray_xml]
							]
				       );

	my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);

	print $results_fh "$result\n";

	my @diff_results = qx/diff $control_data_dir\/$runGOstat_control_file $results_file/;

	if ((@diff_results > 0) && (! ($diff_results[1] =~ /date/))) {
	    print STDERR "runGOstat service failed!\n";
	    print STDERR "diff_results: @diff_results\n";

	    close $results_fh;
	    unlink $results_file;

	}
	else {

	    print  "runGOstat okay...\n";

	    close $results_fh;
	    unlink $results_file;
	}

    }

}

my $t2 = Benchmark->new ();
print   "\nTotal : ", timestr (timediff ($t2, $t1)), "\n";

#
# End
#

sub MobyServiceInstantiation {

    my ($central, $serviceName, $serviceAuthority) = @_;

    ##################################################################
    #
    # Moby Service instantiation
    #
    ##################################################################

    my ($service_instances, $reg) = $central->findService (
							   serviceName  => $serviceName,
							   authURI      => $serviceAuthority,
							   );

    if ($_debug) {
	print STDERR "Service instances references: " . @$service_instances . "\n";
    }

    my $service_instance = $service_instances->[0];

    if (@$service_instances > 1) {
	print STDERR "more than one service instance for service, $serviceName, at \"$serviceAuthority\"\n";
    }

    my $wsdl = $C->retrieveService($service_instance);

    if ($_debug) {
	print STDERR "wsdl: $wsdl\n";
    }

    if (!$wsdl || ($wsdl !~ /\<definitions/)){
	# print STDERR "test \t\t[FAIL]\tWSDL was not retrieved\n\n";
    }

    my $service = MOBY::Client::Service->new(service => $wsdl);

    if (! (defined $service)) {
	print STDERR "Error, cound not instanciate MOBY::Client::Service object for service, $serviceName, at \"$serviceAuthority\"\n";
    }

    return $service;

}
