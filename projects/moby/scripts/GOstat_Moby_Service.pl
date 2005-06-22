#!/usr/local/bin/perl -w

##################################################################
#
# GOstat Moby Service Client
#
##################################################################

use strict;
use Data::Dumper;

# BioMoby and SOAP libraries

use MOBY::Client::Central;
use MOBY::Client::Service;
use SOAP::Lite;
# use SOAP::Lite + 'trace';

# Benchmark Module
use Benchmark;
##################################################################

my $t1 = Benchmark->new ();

my $_debug = 0;

my $serviceName   = "runGOstat";
my $articleName_1 = "regulated genes";
my $articleName_2 = "array genes";

# URI
$::authURI = 'genome.imim.es';

##################################################################
#
# Setup Moby configuration parameters
#
##################################################################

my $URL = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';
my $URI = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://chirimoyo.ac.uma.es/MOBY/Central';
my $PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

# CNB - Production

# $URL = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://www.inab.org/cgi-bin/MOBY-Central.pl';
# $URI = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://www.inab.org/MOBY/Central';
# $PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

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

my $in_file_1  = shift @ARGV || "/home/ug/arnau/projects/gostat/data/mut1_downreg.fbgn";
my $in_file_2  = shift @ARGV || "/home/ug/arnau/projects/gostat/data/allArray.fbgn";
my $datasource = "FB";

my @regulated_genes_xml = ();
my @array_genes_xml     = ();

open (FILE, "<$in_file_1") or die "can't open file with regulated genes, $in_file_1!\n";

while (<FILE>) {
    my $seq_id = $_;
    chomp ($seq_id);

    if ($_debug) {
	print STDERR "Sequence identifier, $seq_id\n";
    }
    
    #
    # Sequence Input in XML format
    #
    
    my $regulated_genes_xml = <<PRT;
<Object namespace="$datasource" id="$seq_id">
</Object>
PRT

push (@regulated_genes_xml, $regulated_genes_xml);

} # Next sequence identifier

close FILE;

open (FILE, "<$in_file_2") or die "can't open file with all array genes, $in_file_2!\n";;

while (<FILE>) {
    my $seq_id = $_;
    chomp ($seq_id);

    if ($_debug) {
	print STDERR "Sequence identifier, $seq_id\n";
    }

    #
    # Sequence Input in XML format
    #
    
    my $array_genes_xml = <<PRT;
<Object namespace="$datasource" id="$seq_id">
</Object>
PRT

push (@array_genes_xml, $array_genes_xml);

} # Next sequence identifier

close FILE;
	
#
# Parameters
#

# ...

##################################################################
#
# Service execution
#
##################################################################

my $result = $Service->execute(XMLinputlist => [
						[$articleName_1, \@regulated_genes_xml, $articleName_2, \@array_genes_xml]
						]);

##################################################################
#
# Result processing
#
##################################################################

if (defined $result) {
    # NB : it doesn't it worked, a way to know is by using SOAP::Lite +trace
    # Any better way ??????
    
    print "result\n", $result, "\n";
}
else {
    print STDERR "The service didn't return any results!\n";
}

my $t2 = Benchmark->new ();
print  STDERR "\nTotal : ", timestr (timediff ($t2, $t1)), "\n";
