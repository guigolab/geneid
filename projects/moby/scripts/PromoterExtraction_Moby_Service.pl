#!/usr/local/bin/perl -w

##################################################################
#
# Upstream Sequence retrieval from Ensembl Moby Service Client
#
##################################################################

use strict;
use Data::Dumper;

# BioMoby and SOAP libraries

use MOBY::Client::Central;
use MOBY::Client::Service;
# use SOAP::Lite + 'trace';

# Benchmark Module
use Benchmark;
##################################################################

my $t1 = Benchmark->new ();

my $_debug = 0;

my $serviceName   = "getUpstreamSeqfromEnsembl";
my $articleName = "genes";

# URI
$::authURI = 'genome.imim.es';
# $::authURI = 'www.pcm.uam.es';

print STDERR "\nrequesting, $::authURI...\n\n";

##################################################################
#
# Setup Moby configuration parameters
#
##################################################################

my $URL = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';
my $URI = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://chirimoyo.ac.uma.es/MOBY/Central';
my $PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

# CNB - Production

$URL = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://www.inab.org/cgi-bin/MOBY-Central.pl';
$URI = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://www.inab.org/MOBY/Central';
$PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

# icapture

$URI = $ENV{MOBY_URI}='http://mobycentral.icapture.ubc.ca/MOBY/Central';
$URL = $ENV{MOBY_SERVER}='http://mobycentral.icapture.ubc.ca/cgi-bin/MOBY05/mobycentral.pl';
$PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

##################################################################
#
# Moby Server Instanciation
#
##################################################################

my $C = MOBY::Client::Central->new(
				   Registries => {mobycentral => {URL => $URL,URI => $URI}}
				   );
    
print STDERR "TESTING MOBY CLIENT with\n\tURL: $URL\n\tURI: $URI\n\tProxy: $PROXY\n\n";

##################################################################
#
# Moby Service instantiation
#
##################################################################

my ($service_instances, $reg) = $C->findService (
						 serviceName => $serviceName,
						 authURI     => $::authURI,
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
# Setup sequence input data
#
##################################################################

my $in_file  = shift @ARGV || "/home/ug/arnau/data/promoterExtraction/Homo_sapiens.lst";
my $datasource = "Ensembl";

if ((! -f $in_file)) {
    print STDERR "Error, can't find one or both of these files, $in_file!\n";
    exit 0;
}

my $genes = qx/cat $in_file/;

my $genes_xml = <<PRT;
<text-formatted namespace="$datasource" id="id">
<String id='' articleName='content'>
<![CDATA[
$genes
]]>
</String>
</text-formatted>
PRT

#
# Parameters
#

my $organism          = "homo sapiens";
# $organism             = "Drosophila melanogaster";
# $organism             = "Mus musculus";
my $dbrelease         = "35";
my $upstream_length   = 1000;
my $downstream_length = 0;
my $intergenic_only   = "true";
my $orthologous_mode  = "false";

my $organism_xml = <<PRT;
<Value>$organism</Value>
PRT

my $dbrelease_xml = <<PRT;
<Value>$dbrelease</Value>
PRT

my $upstream_length_xml = <<PRT;
<Value>$upstream_length</Value>
PRT

my $downstream_length_xml = <<PRT;
<Value>$downstream_length</Value>
PRT

my $intergenic_only_xml = <<PRT;
<Value>$intergenic_only</Value>
PRT

my $orthologous_mode_xml = <<PRT;
<Value>$orthologous_mode</Value>
PRT

##################################################################
#
# Service execution
#
##################################################################

# my $result = $Service->execute(XMLinputlist => [
#     [ $articleName, $genes_xml, 
#       'organism', $organism_xml, 'upstream length', $upstream_length_xml, 'downstream length', $downstream_length_xml, 'intergenic only', $intergenic_only_xml, 'orthologous mode', $orthologous_mode_xml]
#    ]);

my $result = $Service->execute(XMLinputlist => [
                                                [ $articleName, $genes_xml]
                                               ]
                              );

##################################################################
#
# Result processing
#
##################################################################

if (defined $result) {
    print STDERR "result\n";
    print $result;
    print STDERR "\n";
}
else {
    print STDERR "The service didn't return any results!\n";
}

my $t2 = Benchmark->new ();
print  STDERR "\nTotal : ", timestr (timediff ($t2, $t1)), "\n";
