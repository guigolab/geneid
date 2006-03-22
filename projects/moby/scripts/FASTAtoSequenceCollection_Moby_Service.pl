#!/usr/local/bin/perl -w

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
use SOAP::Lite;
# use SOAP::Lite + 'trace';

# Benchmark Module
use Benchmark;
##################################################################

sub help {
return <<"END_HELP";
Description: Converts FASTA sequences to a collection of DNASequence or AminoAcidSequence objects Moby services available from genome.imim.es
Usage:

    FASTAtoSequenceCollection_Moby_Service.pl [-h] -x {Moby Central} -f {sequence FASTA file} -t {sequence type}
	-h help
	-x MOBY Central: Chirimoyo, Mobydev, Inab or BioMoby
		<1> or Chirimoyo
		<2> or Mobydev
		<3> or Inab
		<4> or BioMoby
	-f Sequence(s) input file, in FASTA format
	-t Sequence type [DNASequence, AminoAcidSequence]
	
Examples using some combinations:
	perl FASTAtoSequenceCollection_Moby_service.pl -x 2 -f /home/ug/arnau/data/AC005155.fa -t DNASequence

END_HELP

}

BEGIN {
	
	# Determines the options with values from program
	use vars qw/$opt_h $opt_x $opt_f $opt_t/;
	   
	# these are switches taking an argument (a value)
	my $switches = 'hxft';
	   
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

my $articleName = "sequences";
$::authURI = 'genome.imim.es';
my $sequence_type = $opt_t || "DNASequence";

my $in_file    = $opt_f || "/home/ug/arnau/data/AB065701.fa";
if (not (-f $in_file)) {
    print STDERR "Error, can't find file, $in_file\n";
    exit 0;
}
my $datasource = "EMBL";

my $serviceName;
if ($sequence_type eq "DNASequence") {
    $serviceName = "fromFASTAToDNASequenceCollection";
}
elsif ($sequence_type eq "AminoAcidSequence") {
    $serviceName = "fromFASTAToAminoAcidSequenceCollection";
}
else {
    print STDERR "sequence type, $sequence_type, unknown!\n";
    exit 0;
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

my $fasta_sequences = qx/cat $in_file/;

# <String namespace="" id="" articleName="content">
# </String>

my $fasta_xml = <<PRT;
<FASTA_NA_multi namespace="$datasource" id="">
<![CDATA[
$fasta_sequences
]]>
</FASTA_NA_multi>
PRT

if ($_debug) {
    print STDERR "input xml,\n $fasta_xml\n";
}

##################################################################
#
# Service execution
#
##################################################################

my $result = $Service->execute (XMLinputlist => [
						 ["$articleName", $fasta_xml]
						]);

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
