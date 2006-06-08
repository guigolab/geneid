#!/usr/local/bin/perl -w

##################################################################
#
# GeneID Moby Service Client
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

# MIME encoding/decoding
use MIME::Base64;

# Bioperl Libraries

use Bio::SeqIO;

# Benchmark Module
use Benchmark;
##################################################################

sub help {
return <<"END_HELP";
Description: Execute filterSequencesByLength Moby services available from genome.imim.es to filter out small sequences
Usage:

    Phrap_Moby_Service.pl [-h] -x {Moby Central} -s {Service Name} -f {sequences FASTA file}
	-h help
	-x MOBY Central: Chirimoyo, Mobydev, Inab or BioMoby
		<1> or Chirimoyo
		<2> or Mobydev
		<3> or Inab
		<4> or BioMoby
	-s Service Name
	-i Sequence(s) input file, in FASTA format - optional
	
Examples using some combinations:
	perl FilteringByLength_Moby_Service.pl -x 2 -s filterSequencesByLength -f /home/ug/arnau/data/Assembly_Galicia/seqs.fa

END_HELP

}

BEGIN {
	
	# Determines the options with values from program
	use vars qw/$opt_h $opt_x $opt_s $opt_f/;
	   
	# these are switches taking an argument (a value)
	my $switches = 'hxsf';
	   
	# Get the switches
	getopt($switches);
	
	# If the user does not write nothing, skip to help
	if (defined($opt_h) || !defined($opt_x)){
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

my $serviceName = "filterSequencesByLength";
defined $opt_s and $serviceName = $opt_s;
my $sequences_articleName       = "sequences";
my $quality_articleName         = "base_quality_data";
my $input_sequences_object_type = "FASTA_NA_multi";
my $input_quality_object_type   = "FASTA_Base_Quality_multi";
$::authURI = 'genome.imim.es';

my $in_file    = $opt_f || "/home/ug/arnau/data/Assembly_Galicia/seqs.fa";
my $in_quality_file = $in_file . ".qual";
my $datasource = "INB";
my $seqId = "assembly_00001";

if (! -f $in_file) {
    die "can't find , input file, $in_file!\n";
}

if ($serviceName eq "runPhrapWithQualityData") {
    if (! -f $in_quality_file) {
	die "can't find , input quality file, $in_quality_file!\n";
    }
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
    print STDERR "test \t\t[FAIL]\tWSDL was not retrieved\n\n";
}

my $Service = MOBY::Client::Service->new(service => $wsdl);

##################################################################
#
# Instanciate a Bioperl Sequence Factory to parse the sequence entries from the FASTA file
#
##################################################################

#
# Sequence Input
#

my $input_sequences_data = qx/cat $in_file/;

my $input_sequences_data_b64;
# Compression
# my $system_result = qx/gzip $in_file/;
# $input_sequences_data_b64 = qx/cat $in_file".gz"/;

# Binary encoding
# $input_sequences_data_b64 = encode_base64 ($input_sequences_data);

my $input_sequences_xml = <<PRT;
<$input_sequences_object_type namespace="$datasource" id="$seqId">
  <String namespace="" id=""  articleName="content"><![CDATA[$input_sequences_data]]></String>
</$input_sequences_object_type>
PRT

undef $input_sequences_data;
undef $input_sequences_data_b64;

#
# Parameters (secondary articles)
#

my $trim_masked_regions = "On";
my $length_cutoff       = 200;

my $trim_masked_regions_xml = "<Value>$trim_masked_regions</Value>";
my $length_cutoff_xml       = "<Value>$length_cutoff</Value>";

my $results;
if ($serviceName eq "filterSequencesByLength") {
    $results = $Service->execute(XMLinputlist => [
						 ["$sequences_articleName", $input_sequences_xml, 'trim_masked_regions', $trim_masked_regions_xml, 'length_cutoff', $length_cutoff_xml]
						]);
}
else {
    my $input_quality_data = qx/cat $in_quality_file/;

    my $input_quality_xml = <<PRT;
<$input_quality_object_type namespace="$datasource" id="$seqId">
<String namespace="" id=""  articleName="content"><![CDATA[$input_quality_data]]></String>
</$input_quality_object_type>
PRT

    undef $input_quality_data;

    $results = $Service->execute(XMLinputlist => [
						["$sequences_articleName", $input_sequences_xml, "$quality_articleName", $input_quality_xml,, 'trim_masked_regions', $trim_masked_regions, 'length_cutoff', $length_cutoff]
					       ]);
}

undef $input_sequences_xml;

##################################################################
#
# Result processing
#
##################################################################

print STDERR "results:\n";
print "$results\n";

my $t2 = Benchmark->new ();
print STDERR "\nTotal : ", timestr (timediff ($t2, $t1)), "\n";
