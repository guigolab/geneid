#!/usr/local/bin/perl -w

##################################################################
#
# GOstat Moby Service Client
#
##################################################################

use strict;

# BioMoby and SOAP libraries

use MOBY::Client::Central;
use MOBY::Client::Service;

# use SOAP::Lite + 'trace';

# Benchmark Module
use Benchmark;
##################################################################

# Bioperl
use Bio::SeqIO;

my $t1 = Benchmark->new ();

my $_debug = 0;

my $serviceName   = "translateGeneIDGFFPredictions";
my $articleName_1 = "sequences";
my $articleName_2 = "geneid_predictions";

my $translation_table = "Standard (1)";

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

$URL = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://www.inab.org/cgi-bin/MOBY-Central.pl';
$URI = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://www.inab.org/MOBY/Central';
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

if (not defined $Service) {
    print STDERR "Error, service not instanciated!\n";
    exit 0;
}

##################################################################
#
# Setup sequence input data
#
##################################################################

my $in_file_1  = shift @ARGV || "/home/ug/arnau/data/AC005155.fa";
my $in_file_2  = shift @ARGV || "/home/ug/arnau/data/AC005155.geneid.gff.out";
my $datasource = "EMBL";

my $sequence_xml;
my @sequences_xml = ();

my $seqin = Bio::SeqIO->new (
			     -file   => $in_file_1,
			     -format => 'fasta',
			     );

while (my $seqobj = $seqin->next_seq) {
    my $nucleotide = $seqobj->seq;
    my $seq_id = $seqobj->display_id;
    my $lnucleotide = length($nucleotide);

    #
    # Sequence Input in XML format
    #

    $sequence_xml = <<PRT;
<DNASequence namespace="$datasource" id="$seq_id">
  <Integer namespace="" id="" articleName="Length">$lnucleotide</Integer>
  <String namespace="" id=""  articleName="SequenceString">$nucleotide</String>
</DNASequence>
PRT

    push (@sequences_xml, $sequence_xml);

} # Next sequence

# Hardcoded, should match the query sequence !!
my $seq_id = "AC005155";

my @predictions_xml = ();

my @GFFs = qx/cat $in_file_2/;

my $GFF = "";
foreach my $GFF_line (@GFFs) {
    $GFF .= "$GFF_line";
}
undef @GFFs;

chomp $GFF;

my $prediction_xml = <<PRT;
<GFF namespace="$datasource" id="$seq_id">
<![CDATA[
$GFF
]]>
</GFF>
PRT

$prediction_xml = <<PRT;
<GFF namespace="$datasource" id="$seq_id">
<String namespace="" id="" articleName="Content">
<![CDATA[
$GFF
]]>
</String>
</GFF>
PRT

push (@predictions_xml, $prediction_xml);

#
# Parameters
#

my $translation_table_parameter_xml = <<PRT;
<Value>$translation_table</Value>
PRT

##################################################################
#
# Service execution
#
##################################################################

    print STDERR "Executing Moby request...\n";

my $result = $Service->execute(XMLinputlist => [
						[$articleName_1, $sequence_xml, $articleName_2, $prediction_xml, "translation table", $translation_table]
						]);

##################################################################
#
# Result processing
#
##################################################################

if (defined $result) {
    # NB : it doesn't it worked, a way to know is by using SOAP::Lite +trace
    # Any better way ??????
    
    print STDERR "result\n";
    print "$result\n";
    print STDERR "\n";
}
else {
    print STDERR "The service didn't return any results!\n";
}

my $t2 = Benchmark->new ();
print  STDERR "\nTotal : ", timestr (timediff ($t2, $t1)), "\n";
