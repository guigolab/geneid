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

# Benchmark Module
use Benchmark;
##################################################################

sub help {
return <<"END_HELP";
Description: Execute FASTA to a collection of DNASequence objects Moby services available from genome.imim.es
Usage:

    GenericSequencetoFASTA_Moby_Service.pl [-h] -x {Moby Central} -s {Service Name} -f {sequence FASTA file}
	-h help
	-x MOBY Central: Chirimoyo, Xistral, Inab or BioMoby
		<1> or Chirimoyo
		<2> or Xistral
		<3> or Inab
		<4> or BioMoby
	-s Service Name - not required right now as only one!
	-f Sequence(s) input file, in FASTA format - optional
	
Examples using some combinations:
	perl FASTAtoDNASequenceCollection_Moby_service.pl -x 1 -f /home/ug/arnau/data/AC005155.fa -c config.txt

END_HELP

}

BEGIN {
	
	# Determines the options with values from program
	use vars qw/$opt_h $opt_x $opt_f $opt_c/;
	   
	# these are switches taking an argument (a value)
	my $switches = 'hxfc';
	   
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

# input file

my $in_file = $opt_f || "~/data/promoterExtraction/Homo_sapiens.fa";
if (not (-f $in_file)) {
    print STDERR "Error, can't find input file, $in_file\n";
    exit 1;
}

# parameter file

# generate a hash of hash
# $parameters{$serviceName} = $serviceName_parameters  

my %parameters;
my $serviceName = undef;

my $config_file = $opt_c || "~/config/config.txt";
if (not (-f $config_file)) {
    print STDERR "Error, can't find config file, $config_file\n";
    exit 1;
}

open CONFIG, "<$config_file" or die "can't open config file, $config_file\n";
while (<CONFIG>) {
  my $line = $_;
  if ($line =~ /^>(.+)/) {
    $serviceName = $1;
    chomp $serviceName;
    if (! defined $serviceName) {
      print STDERR "pb when parsing config file - can't get the service name\n";
      print STDERR "line , $line\n";
    }
    else {
      print STDERR "parsed following article name, $serviceName\n";
    }
  }
  elsif ($line =~ /(.+)=(.+)/) {
  
    print STDERR "parsing a parameter...\n";
  
    if (! defined $serviceName) {
      print STDERR "pb - no service name defined!\n";
    }
  
    my $parameter_name  = $1;
    my $parameter_value = $2;
    my $parameters_tmp;
    if (defined $parameters{$serviceName}) {
      $parameters_tmp = $parameters{$serviceName};
    }
    $parameters_tmp->{$parameter_name} = $parameter_value;
  }
}
close CONFIG;

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

# Get input information from config file
# works this way ????
my $namespace = $parameters{input}->{namespace};
print STDERR "namespace, $namespace\n";
my $input_type = $parameters{input}->{input_type};

if (! ($input_type eq "FASTA" || $input_type eq "identifiers")) {
	print STDERR "don't know about this input type, $input_type!\n";
	print STDERR "must be either 'FASTA' or 'identifiers'\n";
	exit 1;
}

##################################################################
#
# Moby Service instantiation
#
##################################################################

# 1/ fromFASTAtoDNASequenceCollection

if ($input_type eq "FASTA") {
	# fromFASTAtoDNASequenceCollection
	$serviceName = "fromFASTAtoDNASequenceCollection";
}
else {
	# getUpstreamSeqfromEnsembl
	$serviceName = "fromFASTAtoDNASequenceCollection";
}

# Get the parameters from the configuration file

my $authURI           = "";
my $articleName       = "";
my $species           = "";
my $upstream_length   = "";
my $downstream_length = "";
my $intergenic_only   = "";

if ($_debug) {
    print STDERR "finding service, $serviceName...\n";
}

my ($service_instances, $reg) = $C->findService (
						 serviceName  => "$serviceName",
						 authURI      => $authURI
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

my $fasta_xml = <<PRT;
<FASTA namespace="$namespace" id="">
<![CDATA[
$fasta_sequences
]]>
</FASTA>
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

if ($_debug) {
	print STDERR "result\n";
	print STDERR $result;
	print STDERR "\n";
}

# runMatScanGFFCollection

$serviceName = "runMatScanGFFCollection";
# ...

# runMultiMetaAlignment

$serviceName = "runMultiMetaAlignment";
# ...

# generateScoreMatrix

$serviceName = "generateScoreMatrix";
# ...

# inbHierarchicalCluster

$serviceName = "inbHierarchicalCluster";
# ...

# inbTreeView

$serviceName = "inbTreeView";
# ...

my $t2 = Benchmark->new ();
print STDERR "\nTotal : ", timestr (timediff ($t2, $t1)), "\n";
