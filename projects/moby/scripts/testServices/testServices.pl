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

# BioMoby and SOAP libraries

use MOBY::Client::Central;
use MOBY::Client::Service;
# use SOAP::Lite + 'trace';

use File::Temp qw/tempfile/;

sub help {
return <<"END_HELP";
Description: Execute GeneID Moby services available from genome.imim.es
Usage:

    testServices.pl [-h] -x {Moby Central}
	-h help
	-x MOBY Central: Chirimoyo, Xistral, Inab or BioMoby
		<1> or Chirimoyo
		<2> or Xistral
		<3> or Inab
		<4> or BioMoby
	-s Service Name
	
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

my $input_data_dir            = "./inputData";
my $control_data_dir          = "./outputControl";

my $nucleotide_sequence_xml_file = "Hsap_BTK.msk.xml";
my $tblastx_output_xml_file   = "Hsap_BTK.tbx.xml";
my $GeneIDGFF_xml_file        = "Hsap_BTK.msk.GeneIDGFF.xml";
my $geneIds_lst_xml_file      = "geneIds.lst.xml";
my $gostat_regulated_xml_file = "mut1_downreg.fbgn.xml";
my $gostat_allArray_xml_file  = "allArray.fbgn.xml";

# Check that the files exist !!!

if ((not -f "$input_data_dir/$nucleotide_sequence_xml_file") || (not -f "$input_data_dir/$tblastx_output_xml_file") || (not -f "$input_data_dir/$geneIds_lst_xml_file") || (not -f "$input_data_dir/$GeneIDGFF_xml_file") || (not -f "$input_data_dir/$gostat_regulated_xml_file") || (not -f "$input_data_dir/$gostat_allArray_xml_file")) {
    print STDERR "Error, can't find one of the input files in directory, $input_data_dir!\n";
    exit 0;
}

my $tblastx_output_xml   = qx/cat $input_data_dir\/$tblastx_output_xml_file/;
my $nucleotide_sequence_xml = qx/cat $input_data_dir\/$nucleotide_sequence_xml_file/;
my $GeneIDGFF_xml        = qx/cat $input_data_dir\/$GeneIDGFF_xml_file/;
my $geneIds_lst_xml      = qx/cat $input_data_dir\/$geneIds_lst_xml_file/;
my $gostat_regulated_xml = qx/cat $input_data_dir\/$gostat_regulated_xml_file/;
my $gostat_allArray_xml  = qx/cat $input_data_dir\/$gostat_allArray_xml_file/;

my $runGeneID_control_file            = "Hsap_BTK.msk.runGeneID.control";
my $runGeneIDGFF_control_file         = "Hsap_BTK.msk.runGeneIDGFF.control";
my $runSGP2GFF_control_file           = "Hsap_BTK.msk.runSGP2GFF.control";
my $translateGeneIDGFFPredictions_control_file = "Hsap_BTK.msk.GeneIDGFF.translateGeneIDGFFPredictions.control";
my $getUpstreamSeqfromEnsembl_control_file = "geneIds.lst.getUpstreamSeqfromEnsembl.control";
my $runGOstat_control_file = "mut1_downreg.fbgn.runGOstat.control";

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
	
    }elsif (($opt_x == 2) || ($opt_x eq 'Xistral')) {
	
	$URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://xistral/cgi-bin/MOBY-Central.pl';
	$URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://xistral/MOBY/Central';
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

# Execute runGeneID Web service

print STDERR "testing runGeneID...\n\n";

my $service = MobyServiceInstantiation ($C, "runGeneID", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['sequences', $nucleotide_sequence_xml]
						    ]
				   ) ;
    
    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);
    
    print $results_fh "$result\n";
    
    my @diff_results = qx/diff $control_data_dir\/$runGeneID_control_file $results_file/;
    
    # The execution date report is obviously different so don't check this !!

    if ((@diff_results > 0) && (! ($diff_results[1] =~ /date/))) {

	# Don't check either the temporary file name which is also different !

	if (! ($diff_results[5] =~ /geneid\d+monstre1.imim.es.fasta/)) {
	    
	    print STDERR "runGeneID failed!\n";
	    print STDERR "diff_results: @diff_results\n";
	    
	    close $results_fh;
	    unlink $results_file;
	    
	}
	else {
	
	    print STDERR "runGeneID okay...\n";
	    
	    close $results_fh;
	    unlink $results_file;
	}
    }
    else {
	
	print STDERR "runGeneID okay...\n";
	
	close $results_fh;
	unlink $results_file;
    }
}

# Execute runGeneIDGFF Web service

print STDERR "\ntesting runGeneIDGFF...\n\n";

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
	
	print STDERR "runGeneIDGFF okay...\n";
	
	close $results_fh;
	unlink $results_file;
    }
}

# Execute translateGeneIDGFFPredictions Web service

print STDERR "\ntesting translateGeneIDGFFPredictions...\n\n";

$service = MobyServiceInstantiation ($C, "translateGeneIDGFFPredictions", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['sequences', $nucleotide_sequence_xml, 'geneid_predictions', $GeneIDGFF_xml]
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
	
	print STDERR "translateGeneIDGFFPredictions okay...\n";
	
	close $results_fh;
	unlink $results_file;
    }
}

# Execute getUpstreamSeqfromEnsembl Web service

print STDERR "\ntesting getUpstreamSeqfromEnsembl...\n\n";

$service = MobyServiceInstantiation ($C, "getUpstreamSeqfromEnsembl", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['genes', $geneIds_lst_xml]
						    ]
				   );
    
    my ($results_fh, $results_file) = tempfile ("/tmp/MOBY_RESULTS.XXXXX", UNLINK => 0);

    print $results_fh "$result\n";
    
    my @diff_results = qx/diff $control_data_dir\/$getUpstreamSeqfromEnsembl_control_file $results_file/;
    
    if (@diff_results > 0) {
	print STDERR "getUpstreamSeqfromEnsembl service failed!\n";
	print STDERR "diff_results: @diff_results\n";
	
	close $results_fh;
	unlink $results_file;
	
    }
    else {
	
	print STDERR "getUpstreamSeqfromEnsembl okay...\n";
	
	close $results_fh;
	unlink $results_file;
    }
}

# Execute runSGP2GFF Web service

print STDERR "\ntesting runSGP2GFF...\n\n";

$service = MobyServiceInstantiation ($C, "runSGP2GFF", $AUTH);
if (defined $service) {
    my $result = $service->execute(
				   XMLinputlist => [
						    ['sequences', $nucleotide_sequence_xml, 'tblastx', $tblastx_output_xml]
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
	
	print STDERR "runSGP2GFF okay...\n";
	
	close $results_fh;
	unlink $results_file;
    }
    
}

# Execute runGOstat Web service

print STDERR "\ntesting runGOstat...\n\n";

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
	
	print STDERR "runGOstat okay...\n";
	
	close $results_fh;
	unlink $results_file;
    }
    
}

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
        print STDERR "test \t\t[FAIL]\tWSDL was not retrieved\n\n";
    }
    
    my $service = MOBY::Client::Service->new(service => $wsdl);

    if (! (defined $service)) {
	print STDERR "Error, cound not instanciate MOBY::Client::Service object for service, $serviceName, at \"$serviceAuthority\"\n";
    }

    return $service;

}
