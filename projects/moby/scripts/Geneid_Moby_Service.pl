#!/usr/local/bin/perl -w

##################################################################
#
# GeneID Moby Service Client
#
##################################################################

use strict;
use Data::Dumper;

# BioMoby and SOAP libraries

use MOBY::Client::Central;
use MOBY::Client::Service;
use SOAP::Lite;
# use SOAP::Lite + 'trace';

# Bioperl Libraries

use Bio::SeqIO;

my $_debug = 0;

##################################################################
#
# Setup sequence input data
#
##################################################################

my $nucleotide = shift @ARGV;
my $datasource = shift @ARGV || "";

my $in_file    = "/home/ug/arnau/data/AC005155.fa";
$datasource    = "EMBL";
my $profile    = "Human";
my $strands    = "Reverse";

my $seqin = Bio::SeqIO->new (
			     -file => $in_file,
			     -format => 'fasta',
			     );

# Execute GeneID Web service on each individual sequence

# Another way would be to set up a collection of sequences and Run GeneID web service for the whole lot

while (my $seqobj = $seqin->next_seq) {
    my $default_nucleotide = $seqobj->seq;
    
    my $seq_id = $seqobj->display_id;
    
    $nucleotide = $default_nucleotide if !$nucleotide;
    my $lnucleotide = length($nucleotide);
    
    ##################################################################
    #
    # Setup Moby configuration parameters
    #
    ##################################################################

    my $URL = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';
    my $URI = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://chirimoyo.ac.uma.es/MOBY/Central';
    my $PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';
    
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
						     serviceName => "runGeneIDGFF",
						     authURI      => "genome.imim.es",
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
    # Service execution
    #
    ##################################################################

    #
    # Sequence Input
    #

    my $input = <<PRT;
<DNASequence namespace="$datasource" id="$seq_id">
  <Integer namespace="" id="" articleName="Length">$lnucleotide</Integer>
  <String namespace="" id=""  articleName="SequenceString">$nucleotide</String>
</DNASequence>
PRT
	
    #
    # Parameters
    #

    # Profile parameter

    my $profile_xml = <<PRT;
<Value>$profile</Value>
PRT

    # Strands parameter

    my $strand_xml = <<PRT;
<Value>$strands</Value>
PRT

    my $result = $Service->execute(XMLinputlist => [
						    ['sequence', $input, 'profile', $profile_xml, 'strands', $strand_xml]
						   ]);
    
    ##################################################################
    #
    # Result processing
    #
    ##################################################################

    print "result\n", $result, "\n";

}
