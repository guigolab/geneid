#!/usr/local/bin/perl -w

##################################################################
#
# MatScan Moby Service Client
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

# Bioperl Libraries

use Bio::SeqIO;

# Benchmark Module
use Benchmark;
##################################################################

sub help {
return <<"END_HELP";
Description: Execute MatScan Moby services available from genome.imim.es
Usage:

    MatScan_Moby_Service.pl [-h] -x {Moby Central} -s {Service Name} -f {sequence FASTA file}
	-h help
	-x MOBY Central: Chirimoyo, Mobydev, Inab or BioMoby
		<1> or Chirimoyo
		<2> or Mobydev
		<3> or Inab
		<4> or BioMoby
	-s Service Name
	-i Sequence(s) input file, in FASTA format - optional
	
Examples using some combinations:
	perl MatScan_Moby_Service.pl -x 1 -s runMatScanGFF -f /home/ug/arnau/data/promoterExtraction/ENSG00000197785.upstream_region.5000.fa -m /home/ug/arnau/data/MeMe/meme2matrix.text-formatted.out

END_HELP

}

BEGIN {
	
	# Determines the options with values from program
	use vars qw/$opt_h $opt_x $opt_s $opt_f $opt_m/;
	   
	# these are switches taking an argument (a value)
	my $switches = 'hxsfm';
	   
	# Get the switches
	getopt($switches);
	
	# If the user does not write nothing, skip to help
	if (defined($opt_h) || !defined($opt_x) || !defined($opt_s)){
		print help;
		exit 0;
	}
	
}

my $t1 = Benchmark->new ();

my $_debug = 0;

##################################################################
#
# Setup sequence input data and moby parameters
#
##################################################################

my $serviceName = $opt_s;

print STDERR "service name, $serviceName\n";

my $sequence_articleName = "upstream_sequences";
my $matrix_articleName   = "matrix";
$::authURI = 'genome.imim.es';

my $serviceType = "Simple";
if ($serviceName =~ /collection/i) {
    $serviceType = "Collection";
}

my $in_file_1    = $opt_f || "/home/ug/arnau/data/promoterExtraction/ENSG00000197785.upstream_region.5000.fa";
my $in_file_2    = "/home/ug/arnau/data/promoterExtraction/ENSG00000160087.upstream_region.5000.fa";
my $matrix_file  = $opt_m || "/home/ug/arnau/data/MeMe/meme2matrix.text-formatted.out";
my $datasource   = "EMBL";

# Parameters

my $threshold   = "0.8";
my $strands     = "Both";
my $matrix_mode = "log-likelihood";
my $matrix;
if ($serviceName ne "runMatScanGFFCollectionVsInputMatrix") {
    print STDERR "initialising matrix...\n";
    $matrix = "Transfac";
}

##################################################################
#
# Setup Moby configuration parameters
#
##################################################################

##############################################
#
# ASSIGN THE MOBY URI AND MOBY SERVER
#
##############################################

my $URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';
my $URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://chirimoyo.ac.uma.es/MOBY/Central';
my $PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

if (defined($opt_x)) {

	# Delete spaces
	$opt_x =~ s/\s//g;

	# Assign the MOBY Server and MOBY URI
	if (($opt_x == 1) || ($opt_x eq 'Chirimoyo')) {
	    
	    # export MOBY_URI
	    # export MOBY_SERVER
	    $URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';
	    $URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://chirimoyo.ac.uma.es/MOBY/Central';
	    $PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';
	
	}elsif (($opt_x == 2) || ($opt_x eq 'Mobydev')) {
	
		# export MOBY_URI
		# export MOBY_SERVER
		$URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://moby-dev.inab.org/cgi-bin/MOBY-Central.pl';
		$URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://moby-dev.inab.org/MOBY/Central';
		$PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

	}elsif (($opt_x == 3) || ($opt_x eq 'Inab')) {

		# export MOBY_URI
		# export MOBY_SERVER
		$URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://www.inab.org/cgi-bin/MOBY-Central.pl';
		$URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://www.inab.org/MOBY/Central';
		$PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';
	
	}elsif (($opt_x == 4) || ($opt_x eq 'BioMoby')) {

		# export MOBY_URI
		# export MOBY_SERVER
		$URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://mobycentral.icapture.ubc.ca/cgi-bin/MOBY05/mobycentral.pl';
		$URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://mobycentral.icapture.ubc.ca/MOBY/Central';
		$PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

	}else {
		print help;
		exit 0;
	}

}else {
	print help;
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

##################################################################
#
# Moby Service instantiation
#
##################################################################

if ($_debug) {
    print STDERR "finding service, $serviceName...\n";
}

my ($service_instances, $reg) = $C->findService (
						 serviceName  => "$serviceName",
						 authURI      => $::authURI
						 );

if ($_debug) {
    print STDERR "Service instances references: " . @$service_instances . "\n";
}

if (@$service_instances == 0) {
    print STDERR "Error, can't find any service called $serviceName!\n";
    exit 0;
}

my $service_instance = $service_instances->[0];

my $wsdl = $C->retrieveService($service_instance);

if ($_debug) {
    print STDERR "wsdl: $wsdl\n";
}

if (!$wsdl || ($wsdl !~ /\<definitions/)){
    print "test \t\t[FAIL]\tWSDL was not retrieved\n\n";
}

my $Service = MOBY::Client::Service->new(service => $wsdl);

##################################################################
#
# Instanciate a Bioperl Sequence Factory to parse the sequence entries from the FASTA file
#
##################################################################

my $inputs = [];
my $files = [$in_file_1, $in_file_2];
my $i = 0;
while ($i < 2) {
    my $seqin = Bio::SeqIO->new (
				 -file   => $files->[$i],
				 -format => 'fasta',
				 );
    
    while (my $seqobj = $seqin->next_seq) {
	my $nucleotides  = $seqobj->seq;
	my $seq_id       = $seqobj->display_id;
	my $lnucleotides = length($nucleotides);
	
	##################################################################
	#
	# Set up the service input and the secondary articles in XML format 
	#
	##################################################################
	
	#
	# Sequence Input
	#
	
	my $input = <<PRT;
<DNASequence namespace="$datasource" id="$seq_id">
<Integer namespace="" id="" articleName="Length">$lnucleotides</Integer>
<String namespace="" id=""  articleName="SequenceString">$nucleotides</String>
</DNASequence>
PRT
	    
	push (@$inputs, $input);
    }

    $i++;
}

# Input Matrix

my $input_matrix_xml;
if ($serviceName eq "runMatScanGFFCollectionVsInputMatrix") {
    my $input_matrix = qx/cat $matrix_file/;
    $input_matrix_xml = <<PRT;
<text-formatted namespace="MEME" id="id">
<![CDATA[
$input_matrix
]]>
</text-formatted>
PRT
}

#
# Parameters (secondary articles)
#

# Threshold parameter

my $threshold_xml = <<PRT;
<Value>$threshold</Value>
PRT

# Matrix parameter

my $matrix_xml;
if (defined $matrix) {
   $matrix_xml = <<PRT;
<Value>$matrix</Value>
PRT
}

# Matrix parameter

my $matrix_mode_xml = <<PRT;
<Value>$matrix_mode</Value>
PRT

# Strands parameters

my $strands_xml = <<PRT;
<Value>$strands</Value>
PRT

##################################################################
#
# Service execution
#
##################################################################

my $result;

if ($serviceType eq "Collection") {

    print STDERR "Collection mode\n";

    if ($serviceName eq "runMatScanGFFCollection") {
	
	print STDERR "serviceName, $serviceName\n";
	print STDERR "matrix parameter, $matrix_xml\n";
	
	$result = $Service->execute(XMLinputlist => [
						     ["$sequence_articleName", $inputs, 'threshold', $threshold_xml, 'matrix', $matrix_xml, 'matrix mode', $matrix_mode_xml, 'strands', $strands_xml]
						     ]);
    }
    else {

	print STDERR "serviceName, $serviceName\n";
	
	$result = $Service->execute(XMLinputlist => [
						     ["$sequence_articleName", $inputs, "$matrix_articleName", $input_matrix_xml, 'threshold', $threshold_xml, 'matrix mode', $matrix_mode_xml, 'strands', $strands_xml]
						     ]);
    }
}
else {
    
    print STDERR "Simple mode\n";
    print STDERR "matrix parameter, $matrix_xml\n";
    
    my $input = $inputs->[0];
    $result = $Service->execute(XMLinputlist => [
						 ["$sequence_articleName", $input, 'threshold', $threshold_xml, 'matrix' => $matrix_xml, 'matrix mode' => $matrix_mode_xml, 'strands' => $strands_xml]
						 ]);
}

##################################################################
#
# Result processing
#
##################################################################

print STDERR "result\n";
print $result;
print STDERR "\n";

my $t2 = Benchmark->new ();
print STDERR "\nTotal : ", timestr (timediff ($t2, $t1)), "\n";
