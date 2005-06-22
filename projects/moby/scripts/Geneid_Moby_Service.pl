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
use SOAP::Lite;
# use SOAP::Lite + 'trace';

# Bioperl Libraries

use Bio::SeqIO;

# Benchmark Module
use Benchmark;
##################################################################

sub help {
return <<"END_HELP";
Description: Execute GeneID Moby services available from genome.imim.es
Usage:

    retrieveService.pl [-h] -x {Moby Central} -s {Service Name} -i {sequence FASTA file}
	-h help
	-x MOBY Central: Chirimoyo, Xistral, Inab or BioMoby
		<1> or Chirimoyo
		<2> or Xistral
		<3> or Inab
		<4> or BioMoby
	-s Service Name
	-i Sequence(s) input file, in FASTA format - optional
	
Examples using some combinations:
	perl retrieveService.pl -x 1 -s runGeneIDGFF -f /home/ug/arnau/data/AC005155.fa

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
my $articleName = "sequences";
$::authURI = 'genome.imim.es';

my $in_file    = $opt_f || "/home/ug/arnau/data/AC005155.fa";
my $datasource = "EMBL";

# Parameters

my $profile    = "Human";
my $strands    = "Both";
my @exons      = ();
my @signals    = ("Start codons", "Stop codons");


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
		$URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://mobycentral.cbr.nrc.ca/cgi-bin/MOBY05/mobycentral.pl';
		$URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://mobycentral.cbr.nrc.ca/cgi-bin/MOBY05/Central';
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

my $C = MOBY::Client::Central->new(
				   Registries => {mobycentral => {URL => $URL,URI => $URI}}
				   );

if ($_debug) {
    print "TESTING MOBY CLIENT with\n\tURL: $URL\n\tURI: $URI\n\tProxy: $PROXY\n\n";
}

##################################################################
#
# Moby Service instantiation
#
##################################################################

my ($service_instances, $reg) = $C->findService (
						 serviceName  => "$serviceName",
						 authURI      => $::authURI
						 );

if ($_debug) {
    print STDERR "Service instances references: " . @$service_instances . "\n";
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

my $seqin = Bio::SeqIO->new (
			     -file   => $in_file,
			     -format => 'fasta',
			     );

# Execute GeneID Web service on each individual sequence
# Another way would be to set up a collection of sequences and Run GeneID web service for the whole lot
# See Geneid_Moby_Service.v2.pl for this

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

    #
    # Parameters (secondary articles)
    #

    # Profile parameter

    my $profile_xml = <<PRT;
<Value>$profile</Value>
PRT

    # Strands parameter

    my $strand_xml = <<PRT;
<Value>$strands</Value>
PRT

    # Signals parameters

    my $signals_xml = <<PRT;
<Value>None</Value>
PRT

    my $exons_xml = <<PRT;
<Value>None</Value>
PRT

    ##################################################################
    #
    # Service execution
    #
    ##################################################################

# with some reported signals options

#    my $result = $Service->execute(XMLinputlist => [
#						    ["$articleName", $input, 'profile', $profile_xml, 'strands', $strand_xml, 'signals' => $signal_start_codons_xml, 'signals' => $signal_stop_codons_xml]
#						   ]);

# with some reported exons options

    my $result = $Service->execute(XMLinputlist => [
						    ["$articleName", $input, 'profile', $profile_xml, 'strands', $strand_xml, 'signals' => $signals_xml, 'exons' => $exons_xml]
						   ]);

# my $result = $Service->execute(XMLinputlist => [
# 						["$articleName", $input, 'profile', $profile_xml, 'strands', $strand_xml]
# 						]);
    
    ##################################################################
    #
    # Result processing
    #
    ##################################################################

    print "result\n", $result, "\n";

}

my $t2 = Benchmark->new ();
print STDERR "\nTotal : ", timestr (timediff ($t2, $t1)), "\n";
