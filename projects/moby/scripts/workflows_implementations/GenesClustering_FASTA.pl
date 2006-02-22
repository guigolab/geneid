#!/usr/local/bin/perl -w

##################################################################
#
# GenericSequence to FASTA conversion Moby Service Client
#
##################################################################

use strict;

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
Description: Execute a gene clustering workflow, based on patterns found in the upstream regions of a given set of genes. This workflow takes a set of gene upstream sequences in FASTA format and return in STDOUT a clustering tree picture in PNG format.
Usage:

    GenesClustering_FASTA.pl [-h] -x {Moby Central} -f {sequence FASTA file} -t {MatScan threshold} -d {MatScan database} -m {Hierarchical clustering method} -o {Output directory}
	-h help
	-x MOBY Central: Inab, BioMoby, Mobydev (optional - Default is Inab registry)
		<1> or Inab
		<2> or BioMoby
		<3> or Mobydev
	-f Sequence(s) input file, in FASTA format
        -t MatScan probability threshold (Default is 0.85)
        -d MatScan Motifs database [Jaspar, Transfac] (Default is Transfac)
        -m HierarchicalCluster method, e.g nearest neighbour joining or furthest neighbour joining [nearest, furthest] (Default is nearest)
        -o Output directory name, if not specified, the output is turned off, the script will just return a tree clustering picture in STDOUT.
	-c workflow configuration file (default is \$HOME/.workflow.config)

Examples using some combinations:
	perl GenesClustering_FASTA.pl -x 2 -f /home/ug/arnau/data/ENSRNOG00000007726.orthoMode.withRat.1000.fa -c \$HOME/.workflow.config -t 0.80 -d jaspar -m nearest -o output

END_HELP

}

BEGIN {
	
	# Determines the options with values from program
	use vars qw/$opt_h $opt_x $opt_f $opt_c $opt_t $opt_d $opt_m $opt_o/;
	   
	# these are switches taking an argument (a value)
	my $switches = 'hxfctdmo';
	   
	# Get the switches
	getopt($switches);
	
	# If the user does not write nothing, skip to help
	if (defined($opt_h) || !defined ($opt_f)){
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

# input file

my $in_file = $opt_f || $ENV{HOME} . "/data/ENSRNOG00000007726.orthoMode.withRat.1000.fa";
if (not (-f $in_file)) {
    print STDERR "Error, can't find input file, $in_file\n";
    exit 1;
}

# Command-line parameters

my $method    = $opt_m || "nearest";
if (lc ($method) eq "nearest") {
  $method = "Nearest neighbor (single linkage)";
}
elsif (lc ($method) eq "furthest") {
  $method = "Furthest neighbor (complete linkage)";
}
else {
  print STDERR "don't know about the value for method parameter, $method!";
  exit 1;
}
my $threshold = $opt_t || 0.85;
my $database  = $opt_d || "transfac";

# parameter file

# generate a hash of hash
# $parameters{$serviceName} = $serviceName_parameters  

my $serviceName = undef;

my $config_file = $opt_c || $ENV{HOME} . "/.workflow.config";
if (not (-f $config_file)) {
    print STDERR "Error, can't find config file, $config_file\n";
    exit 1;
}

my %parameters = setConfigurationData ($config_file);

# Output
my $output_dir = $opt_o || undef;
if (defined $output_dir) {
  if (-d $output_dir) {
    qx/rm -rf $output_dir/;
  }
  qx/mkdir $output_dir/;
}

##################################################################
#
# Setup Moby configuration parameters
#
##################################################################

# Get input information from config file
# works this way ????

# my $input_parameter = $parameters{input};
my $namespace = $parameters{input}->{namespace};
if ($_debug) {
  print STDERR "namespace, $namespace\n";
}
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

# Default is Inab production registry

my $URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://www.inab.org/cgi-bin/MOBY-Central.pl';
my $URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://www.inab.org/MOBY/Central';
my $PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

if (defined($opt_x)) {
    
	# Delete spaces
        $opt_x =~ s/\s//g;
    
	# Assign the MOBY Server and MOBY URI
	if (($opt_x == 1) || (lc ($opt_x) eq 'inab')) {

	    $URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://www.inab.org/cgi-bin/MOBY-Central.pl';
	    $URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://www.inab.org/MOBY/Central';
	    $PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

	}elsif (($opt_x == 2) || (lc ($opt_x) eq 'biomoby')) {

	    $URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://mobycentral.icapture.ubc.ca/cgi-bin/MOBY05/mobycentral.pl';
	    $URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://mobycentral.icapture.ubc.ca/MOBY/Central';
	    $PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';
		
	}elsif (($opt_x == 3) || (lc ($opt_x) eq 'mobydev')) {

	    $URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://moby-dev.inab.org/cgi-bin/MOBY-Central.pl';
	    $URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://moby-dev.inab.org/MOBY/Central';
	    $PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

	}else {
	    print STDERR "Don't anything about this registry, $opt_x!\n";
	    print help;
	    exit 0;
	}

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

print STDERR "First step, binding sites predictions...\n";
print STDERR "Executing MatScan...\n";

$serviceName        = "runMatScanGFFCollection";
$authURI            = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
$articleName        = $parameters{$serviceName}->{articleName} || die "article name for $serviceName\n";
# Set it up from the command-line now
# my $threshold       = $parameters{$serviceName}->{threshold}   || die "no threshold for $articleName\n";
# my $database        = $parameters{$serviceName}->{database}    || die "no matrix database for $articleName\n";
my $matrix_mode     = $parameters{$serviceName}->{matrix_mode} || die "no matrix mode for $articleName\n";
my $threshold_xml   = "<Value>$threshold</Value>";
my $matrix_xml      = "<Value>$database</Value>";
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
if (defined $output_dir) {
  saveResults ($result, "GFF", "MatScan", $output_dir);
}

if ($_debug) {
	print STDERR "input xml for next service:\n";
	print STDERR join (', ', @$input_xml);
	print STDERR ".\n";
}

print STDERR "First step done\n\n";

# runMultiMetaAlignment

print STDERR "Second step, making the pairwise alignments of the binding site maps...\n";
print STDERR "Executing meta-alignment...\n";

$serviceName = "runMultiMetaAlignment";
$authURI     = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
$articleName = $parameters{$serviceName}->{articleName} || die "article name for $serviceName\n";

$Service = getService ($C, $serviceName, $authURI);

$result = $Service->execute (XMLinputlist => [
					      ["$articleName", $input_xml]
					     ]);

if ($_debug) {
	print STDERR "$serviceName result\n";
	print STDERR $result;
	print STDERR "\n";
}

$input_xml = parseResults ($result, "text-formatted");
if (defined $output_dir) {
  saveResults ($result, "text-formatted", "Meta", $output_dir);
}

if ($_debug) {
	print STDERR "input xml for next service:\n";
	print STDERR join (', ', @$input_xml);
	print STDERR ".\n";
}

print STDERR "Second step done!\n\n";

# generateScoreMatrix

print STDERR "Third step, generating a score matrix by parsing meta-alignment data...\n";

$serviceName = "generateScoreMatrix";
$authURI     = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
$articleName = $parameters{$serviceName}->{articleName} || die "article name for $serviceName\n";

$Service = getService ($C, $serviceName, $authURI);

$result = $Service->execute (XMLinputlist => [
					      ["$articleName", $input_xml]
					     ]);

if ($_debug) {
	print STDERR "$serviceName result\n";
	print STDERR $result;
	print STDERR "\n";
}

my $input_xml_aref = parseResults ($result, "MicroArrayData_Text");
if (defined $output_dir) {
  saveResults ($result, "MicroArrayData_Text", "score_matrix", $output_dir);
}

if ($_debug) {
	print STDERR "input xml for next service:\n";
	print STDERR join (', ', @$input_xml_aref);
	print STDERR ".\n";
}

# convert the input xml into a scalar
$input_xml =  $input_xml_aref->[0];

print STDERR "Third step done!\n\n";

# Request Inab registry for the two following ones
$URL   = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://www.inab.org/cgi-bin/MOBY-Central.pl';
$URI   = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://www.inab.org/MOBY/Central';
$PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

$C = MOBY::Client::Central->new(
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


# inbHierarchicalCluster

print STDERR "Fourth step, gene clustering using a neighbour joining clustering algorithm...\n";

if ($_debug) {
  print STDERR "\nExecuting inbHierarchicalCluster...\n\n";
}

$serviceName   = "inbHierarchicalCluster";
$authURI       = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
$articleName   = $parameters{$serviceName}->{articleName} || "";
# set it up from the command-line now
# my $method     = $parameters{$serviceName}->{method}      || die "no method for $articleName\n";
my $method_xml = "<Value>$method</Value>";

$Service = getService ($C, $serviceName, $authURI);

$result = $Service->execute (XMLinputlist => [
					      ["$articleName", "$input_xml\n", "method", $method_xml]
					     ]);

if ($_debug) {
	print STDERR "\n$serviceName result\n";
	print STDERR $result;
	print STDERR "\n";
}

$input_xml_aref = parseResults ($result, "Clustering");
if (defined $output_dir) {
  saveResults ($result, "Newick_Text", "newick", $output_dir);
}

if ($_debug) {
	print STDERR "input xml for next service:\n";
	print STDERR join (', ', @$input_xml_aref);
	print STDERR ".\n";
}

# convert the input xml into a scalar
$input_xml =  $input_xml_aref->[0];

print STDERR "Clustering done!\n\n";

# inbTreeView

print STDERR "Generating a picture of the clustering tree...\n";

$serviceName   = "inbTreeView";
$authURI       = $parameters{$serviceName}->{authURI}     || die "no URI for $serviceName\n";
$articleName   = $parameters{$serviceName}->{articleName} || "";

$Service = getService ($C, $serviceName, $authURI);

$result = $Service->execute (XMLinputlist => [
					      ["$articleName", $input_xml]
					     ]);

if ($_debug) {
	print STDERR "$serviceName result\n";
	print STDERR $result;
	print STDERR "\n";
}

if (defined $output_dir) {
  saveResults ($result, "b64_Encoded_PNG", "clustering_tree", $output_dir);
}
my $picture_b64_aref = parseTextContent ($result, "b64_Encoded_PNG");
my $picture_b64 = $picture_b64_aref->[0];

# Convert in to a picture and store into a file

my $picture = decode_base64($picture_b64);

print STDERR "Picture done\n\n";
print STDERR "Workflow has terminated.\n";
print $picture;

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
	        if ($_debug) {
  		  print STDERR "parsed following article name, $serviceName\n";
                }
	    }
	}
	elsif ($line =~ /(.+)=(.+)/) {
	    
	    if ($_debug) {
  	      print STDERR "parsing a parameter...\n";
            }
	    
	    if (! defined $serviceName) {
		print STDERR "pb - no service name defined!\n";
	    }
	    
	    my $parameter_name  = $1;
	    my $parameter_value = $2;
	    
	    if ($_debug) {
              print STDERR "parsed name, value, $parameter_name, $parameter_value.\n";
            }
	    
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
	if ($_debug) {
	    print STDERR "test \t\t[FAIL]\tWSDL was not retrieved\n\n";
	}
    }
    
    my $Service = MOBY::Client::Service->new(service => $wsdl);
    
    return $Service;
}


# parse the MobyData bloc to get the collection or simple bloc
# return an array of objects
sub parseResults {
    my ($XML, $object_type) = @_;
    my $inputs_xml = [];
    
    my $parser = XML::LibXML->new();
    my $doc    = $parser->parse_string( $XML );
    $XML       = $doc->getDocumentElement();
    my $elements = $XML->getElementsByTagName( "moby:$object_type" );
    
    my $size = $elements->size();
    
    if ($size == 0) {
	$elements = $XML->getElementsByTagName( "object_type" );
	if ($size == 0) {
	    print STDERR "Error, can't parse the moby output from the moby XML...\n";
	}
    }
    
    my $i = 0;
    while (my $element = $elements->get_node($i)) {
	my $input_xml = $element->toString();
	push (@$inputs_xml, $input_xml);
	
	$i++;
    }

    return $inputs_xml;
}

sub parseTextContent {
    my ($XML, $object_type) = @_;
    my $inputs_text = [];
    
    my $parser = XML::LibXML->new();
    my $doc    = $parser->parse_string( $XML );
    $XML       = $doc->getDocumentElement();
    my $elements = $XML->getElementsByTagName ( "moby:$object_type" );
    my $size   = $elements->size();
    
    if ($size == 0) {
      $elements = $XML->getElementsByTagName ("$object_type");
      if ($size == 0) {
        print STDERR "Error, can't parse the moby output from the moby XML...\n";
      }
    }
    
    my $i = 0;
    while (my $element = $elements->get_node ($i)) {
      my $input_text = $XML->textContent();
      push (@$inputs_text, $input_text);
      $i++;
    }

    return $inputs_text;
}

sub saveResults {
  my ($XML, $object_type, $softwareName, $output_dir) = @_;
  
  if ($_debug) {
    print STDERR "saving results for tool, $softwareName...\n";
  }
  
  my %inputs_text;
    
  my $parser = XML::LibXML->new();
  my $doc    = $parser->parse_string( $XML );
  $XML       = $doc->getDocumentElement();
  my $elements = $XML->getElementsByTagName ( "moby:$object_type" );
  my $size   = $elements->size();
    
  if ($size == 0) {
    $elements = $XML->getElementsByTagName ("$object_type");
    $size = $elements->size();
    if ($size == 0) {
      print STDERR "Error, can't parse the moby output from the moby XML...\n";
    }
  }

  if ($_debug) {
    print STDERR "found $size elements\n";
  }

  # parsing
  
  my $i = 0;
  while ($i < $size) {
    my $element = $elements->[$i];
    my $id = $element->getAttribute ("id");
    
    # if more than one simple object, create an identifier index
    if ((! defined $id) || $id eq "") {
      if ($size > 1) {
        $id = $i + 1;
      }
      else {
        $id = "";
      }
    }
    my $input_text = $element->textContent();
    if ($object_type eq "b64_Encoded_PNG") {
      $input_text = decode_base64($input_text);
    }
    
    $inputs_text{$id} = $input_text;
    $i++;
  }

  my $suffix = getSuffix ($object_type);

  # Save on the local disk the data
  my @ids = keys (%inputs_text);
  if (@ids == 1) {
    # simple => just one file
    my $filename = $softwareName . "." . $suffix;
    my $id = $ids[0];
    if ($id ne "") {
      $filename = $id . "." . $filename;
    }
    
    my $data = $inputs_text{$id};
    open FILE, ">$output_dir/$filename" or die "can' open output file, $output_dir/$filename";
    print FILE "$data";
    close FILE;
  }
  else {
    # collection of simples => a sub directory
    my $subdir = "$softwareName";
    
    if (! -d "$output_dir/$subdir") {
      qx/mkdir $output_dir\/$subdir/;
    }
    
    foreach my $id (@ids) {
      my $filename = $id . "." . $softwareName . "." . $suffix;
      my $data = $inputs_text{$id};
      open FILE, ">$output_dir/$subdir/$filename" or die "can' open output file, $output_dir/$subdir/$filename";
      print FILE "$data";
      close FILE;
    }
  }
}

sub getSuffix {
  my ($object_type) = @_;
  my $suffix = "txt";
  
  if (lc ($object_type) eq "gff") {
    return "gff";
  }
  elsif (lc ($object_type) eq "text-formatted" || lc ($object_type) eq "microarraydata_text" || lc ($object_type) eq "newick_text") {
    return "txt";
  }
  elsif (lc ($object_type) eq "b64_encoded_png") {
    return "png";
  }

  print STDERR "Suffix for object_type, $object_type unknown, return txt!\n"; 
  
  return $suffix;
}
