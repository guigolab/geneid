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

my $tblastx_output_file = "/home/ug/arnau/projects/sgp2/sgp2/samples/Hsap_BTK.tbx";
my $tblastx_output = qx/cat $tblastx_output_file/;

my $in_file    = "/home/ug/arnau/projects/sgp2/sgp2/samples/Hsap_BTK.msk.fa";
$datasource    = "EMBL";

my $seqin = Bio::SeqIO->new (
			     -file   => $in_file,
			     -format => 'fasta',
			     );

# Execute SGP2 Web service on each individual sequence

# Another way would be to set up a collection of sequences and Run SGP2 web service for the whole lot

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

    # Chirimoyo (Development)

    my $URL = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';
    my $URI = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://chirimoyo.ac.uma.es/MOBY/Central';
    my $PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';
    
    # Inab (Production)

    $URL = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://www.inab.org/cgi-bin/MOBY-Central.pl';
    $URI = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://www.inab.org/MOBY/Central';
    $PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

    ##################################################################
    #
    # Moby Server Instanciation
    #
    ##################################################################

    print STDERR "TESTING MOBY CLIENT with\n\tURL: $URL\n\tURI: $URI\n\tProxy: $PROXY\n\n";

    my $C = MOBY::Client::Central->new(
				       Registries => {mobycentral => {URL => $URL,URI => $URI}}
				      );
    
    ##################################################################
    #
    # Moby Service instantiation
    #
    ##################################################################
    
    my ($service_instances, $reg) = $C->findService (
						     serviceName  => "runSGP2GFF",
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

    my $sequence_xml = <<PRT;
<DNASequence namespace="$datasource" id="$seq_id">
  <Integer namespace="" id="" articleName="Length">$lnucleotide</Integer>
  <String namespace="" id=""  articleName="SequenceString">$nucleotide</String>
</DNASequence>
PRT
	

    #
    # TBLASTX report Input
    #

    my $tblastx_xml = <<PRT;
<Blast-Text namespace="$datasource" id="$seq_id">
<![CDATA[
$tblastx_output
]]>
</Blast-Text>
PRT

    my $result = $Service->execute(
			       XMLinputlist => [
						['sequences', $sequence_xml, 'tblastx', $tblastx_xml]
					       ]
			       ) ;
    
    ##################################################################
    #
    # Result processing
    #
    ##################################################################

    print STDERR "result\n";
    print $result;
    print STDERR "\n";

}
