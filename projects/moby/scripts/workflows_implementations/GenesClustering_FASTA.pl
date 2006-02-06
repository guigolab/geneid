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

# XML DOM parsing library
use XML::LibXML;

# MIME encoding/decoding
use MIME::Base64;

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

my $serviceName = undef;

my $config_file = $opt_c || "~/cvs/GRIB/projects/moby/scripts/workflows_implementations/config.txt";
if (not (-f $config_file)) {
    print STDERR "Error, can't find config file, $config_file\n";
    exit 1;
}

my %parameters = setConfigurationData ($config_file);

##################################################################
#
# Setup Moby configuration parameters
#
##################################################################

# Get input information from config file
# works this way ????

# my $input_parameter = $parameters{input};
my $namespace = $parameters{input}->{namespace};
print STDERR "namespace, $namespace\n";
my $input_type = $parameters{input}->{input_type};

if (! ($input_type eq "FASTA" || $input_type eq "identifiers")) {
	print STDERR "don't know about this input type, $input_type!\n";
	print STDERR "must be either 'FASTA' or 'identifiers'\n";
	exit 1;
}

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

##################################################################
#
# Setup the moby objects
#
##################################################################

my $input_xml;

if ($input_type eq "FASTA") {

    # fromFASTAtoDNASequenceCollection

    my $input = qx/cat $in_file/;
    
    $input_xml = <<PRT;
<FASTA namespace="$namespace" id="">
<String articleName="Content">
<![CDATA[
$input
]]>
</String>
</FASTA>
PRT

}
else {

    # getUpstreamSeqfromEnsembl

    my $input = qx/cat $in_file/;
    
    $input_xml = <<PRT;
<text-formatted namespace="$namespace" id="">
<String articleName="Content">
<![CDATA[
$input
]]>
</String>
</text-formatted>
PRT

}

if ($_debug) {
    print STDERR "input xml,\n $input_xml\n";
}

# 1/ fromFASTAtoDNASequenceCollection

if ($input_type eq "FASTA") {
	# fromFASTAtoDNASequenceCollection
	$serviceName = "fromFASTAtoDNASequenceCollection";
}
else {
	# getUpstreamSeqfromEnsembl
	$serviceName = "getUpstreamSeqfromEnsembl";
}

# Get the parameters from the configuration file

my $authURI           = $parameters{$serviceName}->{authURI} || die "no URI\n";
my $articleName       = $parameters{$serviceName}->{articleName} || die "no article name\n";
my $species_xml           = "";
my $upstream_length_xml   = "";
my $downstream_length_xml = "";
my $intergenic_only_xml   = "";

if ($input_type eq "identifiers") {
    my $species           = $parameters{$serviceName}->{species} || die "no species\n";
    my $upstream_length   = $parameters{$serviceName}->{upstream_length}   || die "no upstream length\n";
    my $downstream_length = $parameters{$serviceName}->{downstream_length} || die "no downstream length\n";
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

my $result;
if ($input_type eq "FASTA") {
    $result = $Service->execute (XMLinputlist => [
						  ["$articleName", $input_xml],
						 ]);
}
else {
    $result = $Service->execute (XMLinputlist => [
						  ["$articleName", $input_xml, "species", $species_xml, "upstream_length", $upstream_length_xml, "downstream_length", $downstream_length_xml, "intergenic_only", $intergenic_only_xml]
						 ]);
}

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

$input_xml = parseResults ($result, "DNASequence");

if ($_debug) {
	print STDERR "input xml for next service:\n";
	print STDERR join (', ', @$input_xml);
	print STDERR ".\n";
}

# runMatScanGFFCollection

$serviceName        = "runMatScanGFFCollection";
$authURI            = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
$articleName        = $parameters{$serviceName}->{articleName} || die "article name for $serviceName\n";
my $threshold       = $parameters{$serviceName}->{threshold}   || die "no threshold for $articleName\n";
my $matrix          = $parameters{$serviceName}->{matrix}      || die "no matrix for $articleName\n";
my $matrix_mode     = $parameters{$serviceName}->{matrix_mode} || die "no matrix mode for $articleName\n";
my $threshold_xml   = "<Value>$threshold</Value>";
my $matrix_xml      = "<Value>$matrix</Value>";
my $matrix_mode_xml = "<Value>$matrix_mode</Value>";

$Service = getService ($C, $serviceName, $authURI);

$result = $Service->execute (XMLinputlist => [
					      ["$articleName", $input_xml, "threshold", $threshold_xml, "matrix", $matrix_xml, "matrix mode", $matrix_mode_xml]
					     ]);

if ($_debug) {
	print STDERR "result\n";
	print STDERR $result;
	print STDERR "\n";
}

$input_xml = parseResults ($result, "GFF");

if ($_debug) {
	print STDERR "input xml for next service:\n";
	print STDERR join (', ', @$input_xml);
	print STDERR ".\n";
}

# runMultiMetaAlignment

$serviceName = "runMultiMetaAlignment";
$authURI     = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
$articleName = $parameters{$serviceName}->{articleName} || die "article name for $serviceName\n";

$Service = getService ($C, $serviceName, $authURI);

$result = $Service->execute (XMLinputlist => [
					      ["$articleName", $input_xml]
					     ]);

if ($_debug) {
	print STDERR "result\n";
	print STDERR $result;
	print STDERR "\n";
}

$input_xml = parseResults ($result, "text-formatted");

if ($_debug) {
	print STDERR "input xml for next service:\n";
	print STDERR join (', ', @$input_xml);
	print STDERR ".\n";
}

# generateScoreMatrix

$serviceName = "generateScoreMatrix";
$authURI     = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
$articleName = $parameters{$serviceName}->{articleName} || die "article name for $serviceName\n";

$Service = getService ($C, $serviceName, $authURI);

$result = $Service->execute (XMLinputlist => [
					      ["$articleName", $input_xml]
					     ]);

if (1 || $_debug) {
	print STDERR "$serviceName result\n";
	print STDERR $result;
	print STDERR "\n";
}

my $input_xml_tmp = parseResults ($result, "MicroArrayData_Text");

if ($_debug) {
	print STDERR "input xml for next service:\n";
	print STDERR join (', ', @$input_xml_tmp);
	print STDERR ".\n";
}

# convert the input xml into a scalar
$input_xml =  $input_xml_tmp->[0];

# Request Inab registry for the two following ones
$URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://www.inab.org/cgi-bin/MOBY-Central.pl';
$URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://www.inab.org/MOBY/Central';
$PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

# inbHierarchicalCluster

$serviceName   = "inbHierarchicalCluster";
$authURI       = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
$articleName   = $parameters{$serviceName}->{articleName} || die "article name for $serviceName\n";
my $method     = $parameters{$serviceName}->{method}      || die "no method for $articleName\n";
my $method_xml = "<Value>$method</Value>";

$Service = getService ($C, $serviceName, $authURI);

$result = $Service->execute (XMLinputlist => [
					      ["$articleName", $input_xml, "method", $method]
					     ]);

if (1 || $_debug) {
	print STDERR "result\n";
	print STDERR $result;
	print STDERR "\n";
}

$input_xml_tmp = parseResults ($result, "Clustering");

if ($_debug) {
	print STDERR "input xml for next service:\n";
	print STDERR join (', ', @$input_xml_tmp);
	print STDERR ".\n";
}

# convert the input xml into a scalar
$input_xml =  $input_xml_tmp->[0];

# inbTreeView

$serviceName   = "inbTreeView";
$authURI       = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
$articleName   = $parameters{$serviceName}->{articleName} || die "article name for $serviceName\n";

$Service = getService ($C, $serviceName, $authURI);

$result = $Service->execute (XMLinputlist => [
					      ["$articleName", $input_xml]
					     ]);

if ($_debug) {
	print STDERR "result\n";
	print STDERR $result;
	print STDERR "\n";
}

$input_xml_tmp = parseResults ($result, "b64_Encoded_PNG");

if ($_debug) {
	print STDERR "input xml for next service:\n";
	print STDERR join (', ', @$input_xml_tmp);
	print STDERR ".\n";
}

# convert the input xml into a scalar
$input_xml =  $input_xml_tmp->[0];
my $picture_b64 = parseTextContent ($input_xml);

# Convert in to a picture and store into a file

my $picture = decode_base64($picture_b64);

print STDERR "returning picture data\n";
print $picture;
print STDERR ".\n";

my $t2 = Benchmark->new ();
print STDERR "\nTotal : ", timestr (timediff ($t2, $t1)), "\n";

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
	    
	    print STDERR "parsed name, value, $parameter_name, $parameter_value.\n";
	    
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
    
    return $Service;
}


# parse the MobyData bloc to get the collection or simple bloc
# return an array of objects
sub parseResults {
    my ($XML, $object_type) = @_;
    my $inputs_xml = [];
    
    print STDERR "parsing XML results...\n";

    my $parser = XML::LibXML->new();
    my $doc    = $parser->parse_string( $XML );
    $XML       = $doc->getDocumentElement();
    my $elements = $XML->getElementsByTagName( "moby:$object_type" ) 
	|| $XML->getElementsByTagName( "object_type" );
    
    my $size = $elements->size();
    
    if ($size == 0) {
	print STDERR "Error, can't parse the moby output from the moby XML...\n";
    }
    
    my $i = 0;
    while (my $element = $elements->get_node($i)) {
	my $input_xml = $element->toString();
	
	print STDERR "input_xml, $input_xml\n";
	
	push (@$inputs_xml, $input_xml);
	
	$i++;
    }

    return $inputs_xml;
}

sub parseTextContent {
    my ($XML) = @_;
    my $data;

    my $parser = XML::LibXML->new();
    my $doc    = $parser->parse_string( $XML );
    $XML       = $doc->getDocumentElement();

    $data = $XML->toString();

    return $data;
}
