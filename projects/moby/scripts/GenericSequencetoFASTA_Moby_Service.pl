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
# use SOAP::Lite + 'trace';

# Bioperl Libraries

use Bio::SeqIO;

# Benchmark Module
use Benchmark;
##################################################################

sub help {
return <<"END_HELP";
Description: Execute fromGenericSequencetoFASTA Moby services available from genome.imim.es
Usage:

    GenericSequencetoFASTA_Moby_Service.pl [-h] -x {Moby Central} -s {Service Name} -f {sequence FASTA file} -t {sequence format type}
	-h help
	-x MOBY Central: Chirimoyo, Xistral, Inab or BioMoby
		<1> or Chirimoyo
		<2> or Xistral
		<3> or Inab
		<4> or BioMoby
	-s Service Name
	-f Sequence(s) input file, in FASTA format or XML
	-t sequence format type (xml or fasta)
	
Examples using some combinations:
	perl GenericSequencetoFASTA_Moby_service.pl -x 1 -s fromGenericSequencetoFASTA -f /home/ug/arnau/data/AC005155.fa -t fasta

END_HELP

}

BEGIN {
	
	# Determines the options with values from program
	use vars qw/$opt_h $opt_x $opt_s $opt_f $opt_t/;
	   
	# these are switches taking an argument (a value)
	my $switches = 'hxsft';
	   
	# Get the switches
	getopt($switches);
	
	# If the user does not write nothing, skip to help
	if (defined($opt_h) || !defined($opt_x) || !defined($opt_s) || !defined($opt_t)){
	    print STDERR help;
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
my $articleName = "sequences";
$::authURI = 'genome.imim.es';

my $in_file    = $opt_f || "/home/ug/arnau/data/AC005155.fa";
my $datasource = "EMBL";
my $format_type = $opt_t;
if (not ((lc ($format_type) eq "fasta") || (lc ($format_type) eq "xml"))) {
    print STDERR "don't know about this sequence format type, $format_type\n";
    print STDERR "should be xml or fasta\n";
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
	
	}elsif (($opt_x == 2) || ($opt_x eq 'Xistral')) {
	
		# export MOBY_URI
		# export MOBY_SERVER
		$URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://xistral/cgi-bin/MOBY-Central.pl';
		$URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://xistral/MOBY/Central';
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

my $result;

if ($format_type eq "xml") {

    my $xml_input = qx/cat $in_file/;
    
    ##################################################################
    #
    # Service execution
    #
    ##################################################################
    
    if ($serviceName =~ /collection/i) {
	
	# Collection
	
	# Make a collection
	my $xml_inputs = [$xml_input,$xml_input];
	
	$result = $Service->execute(XMLinputlist => [
						     ["$articleName", $xml_inputs]
						     ]);
    }
    else {
	
	# Simple
	
	$result = $Service->execute(XMLinputlist => [
						     ["$articleName", $xml_input]
						     ]);
    }
}
else {
    my $seqin = Bio::SeqIO->new (
				 -file   => $in_file,
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
	
	my $xml_input = <<PRT;
<DNASequence namespace="$datasource" id="$seq_id">
    <Integer namespace="" id="" articleName="Length">$lnucleotides</Integer>
    <String namespace="" id=""  articleName="SequenceString">$nucleotides</String>
</DNASequence>
PRT
	    
        ##################################################################
        #
        # Service execution
	#
        ##################################################################

        if ($serviceName =~ /collection/i) {
	    
	    # Collection
	    
	    # Make a collection
	    
	    my $namespace_1 = "$datasource" . "_1";
	    my $namespace_2 = "$datasource" . "_2";

	    my $namespace   = "$datasource";
	    
	    my $xml_input_1 = <<PRT;
<DNASequence namespace="$namespace" id="$seq_id">
    <Integer namespace="" id="" articleName="Length">$lnucleotides</Integer>
    <String namespace="" id=""  articleName="SequenceString">$nucleotides</String>
</DNASequence>
PRT
            my $xml_input_2 = <<PRT;
<DNASequence namespace="$namespace" id="$seq_id">
    <Integer namespace="" id="" articleName="Length">$lnucleotides</Integer>
    <String namespace="" id=""  articleName="SequenceString">$nucleotides</String>
</DNASequence>
PRT

	    my $xml_inputs = [$xml_input_1,$xml_input_2];
	    
	    $result = $Service->execute(XMLinputlist => [
							 ["$articleName", $xml_inputs]
							 ]);
	}
	else {
	    
	    # Simple
	    
	    $result = $Service->execute(XMLinputlist => [
							 ["$articleName", $xml_input]
							 ]);
	}
    }
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
