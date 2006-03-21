package INB::GRIB::Utils::CommonUtilsSubs;

use strict;
use warnings;
use Carp;
use Data::Dumper;

use File::Temp qw/tempfile/;
use FileHandle;

# XML DOM parsing library
use XML::LibXML;

# Bioperl
use Bio::Seq;
use Bio::SeqIO;

# Biomoby
use MOBY::CommonSubs qw(:all);
# This module contains the Node type constants
use MOBY::MobyXMLConstants;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use INB::UPC::Blast ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw() ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

# Aqui pondremos las funciones que exportaremos a otros módulos. 
# 
# @EXPORT = qw( &func1 &func2);
# 
our @EXPORT = qw(
  &getTextContentfromXML
  &createSequenceObjectsFromFASTA
  &parseSingleGFFIntoCollectionGFF
  &parseMobySequenceObjectFromDOM
  &parseMobySequenceObjectFromDOMintoBioperlObject
  &convertSequencesIntoFASTA
  &validateDataType
  &getNamespace
  &setMobyResponse
  &MOBY_EMPTY_RESPONSE
  &MOBY_EMPTY_SIMPLE_RESPONSE
  &MOBY_EMPTY_COLLECTION_RESPONSE
  &convert_tabularPositionWeightMatrix_into_MobyMatrix
  &convert_tabularScoreMatrix_into_MobyMatrix
  &convert_MobyMatrix_into_tabularMatrix
);

our $VERSION = '1.0';

# Common subs

sub MOBY_EMPTY_RESPONSE {
    my $self = shift;
    my ($queryID, $output_article_name) = @_;
    return "<moby:mobyData moby:queryID='$queryID'><moby:Simple moby:articleName='$output_article_name'/></moby:mobyData>";
}

sub MOBY_EMPTY_SIMPLE_RESPONSE {
    my $self = shift;
    my ($queryID, $output_article_name) = @_;
    return "<moby:mobyData moby:queryID='$queryID'><moby:Simple moby:articleName='$output_article_name'/></moby:mobyData>";
}

sub MOBY_EMPTY_COLLECTION_RESPONSE {
    my $self = shift;
    my ($queryID, $output_article_name) = @_;
    return "<moby:mobyData moby:queryID='$queryID'><moby:Collection moby:articleName='$output_article_name'/></moby:mobyData>";
}

# Works for both raw text content and CDATA bloc

sub getTextContentFromXML {
    my $self = shift;
    my ( $XML, $tagName ) = @_;
    
    unless ( ref( $XML ) =~ /XML\:\:LibXML/ ) {
	my $parser = XML::LibXML->new();
	my $doc    = $parser->parse_string( $XML );
	$XML       = $doc->getDocumentElement();
    }	

    return '' unless (   
			 my $elements = $XML->getElementsByTagName( "$tagName" ) 
			                || $XML->getElementsByTagName( "moby:$tagName" )
		     );
    my $element = $elements->get_node(0);
    my $content = $element->textContent;

    # Do some cleaning...
    # Get rid of empty lines (those with only spaces or tab...)
    # Get rid of spaces at beginning of lines
    
    $content =~ s/^\s+//;
    $content =~ s/\n\s+/\n/;
    $content =~ s/\s+$//;

    $content =~ s/^\n//g;
    
    # Get rid of empty lines at the end
    while (chomp $content) {
	# just chomp !!
    }
    
    return $content;
}


sub createSequenceObjectsFromFASTA {
    my $self = shift;
    my ($fasta_sequences, $object_type, $namespace) = @_;
    
    my $moby_sequence_objects    = [];
    my $bioperl_sequence_objects = _parse_fasta_sequences ($fasta_sequences);
    
    # Parsing the Bioperl objects returned by _parse_fasta_sequences method
    
    foreach my $seqobj (@$bioperl_sequence_objects) {
	
	# Check the alphabet of the bioperl sequence,
	# if it doesn't match the type of the Moby object, then exit !
	my $alphabet = $seqobj->alphabet();
	
	if ($object_type eq "GenericSequence") {
	    # Don't bother
	}
	elsif ($object_type =~ /AminoAcidSequence|CommentedAASequence/) {
	    if ($alphabet ne "protein") {
		print STDERR "Error FASTA sequences are not protein sequences!\n";
		exit 0;
	    }
	}
	elsif ($object_type eq "NucleotideSequence") {
	    if (not ($alphabet =~ /dna|rna/)) {
		print STDERR "Error FASTA sequences are not nucleotide sequences!\n";
		exit 0;
	    }
	}
	elsif ($object_type =~ /DNASequence/){
	    if (not ($alphabet eq "dna")) {
		print STDERR "Error FASTA sequences are not DNA sequences!\n";
		exit 0;
	    }
	}
	elsif (($object_type =~ /RNASequence/)) {
	    if (not ($alphabet eq "rna")) {
		print STDERR "Error FASTA sequences are not RNA sequences!\n";
		exit 0;
	    }
	}
	else {
	    print STDERR "problem with the fasta sequence alphabet, $alphabet - not matching the returned moby object type, $object_type!\n";
	    exit 0;
	}
	
	my $seq_id     = $seqobj->display_id;
	my $seq_desc   = $seqobj->desc || "";
	my $sequence   = $seqobj->seq;
	my $seq_length = length ($sequence);
	
	my $moby_object;

	if ($object_type =~ /commented/i) {
	    $moby_object = <<PRT;
<moby:$object_type namespace='$namespace' id='$seq_id'>
  <Integer namespace="" id="" articleName="Length">$seq_length</Integer>
  <String namespace="" id=""  articleName="SequenceString">$sequence</String>
  <String namespace="" id=""  articleName="Description">$seq_desc</String>
</moby:$object_type>
PRT
        }
	else {

	    # No description

	    $moby_object = <<PRT;
<moby:$object_type namespace='$namespace' id='$seq_id'>
  <Integer namespace="" id="" articleName="Length">$seq_length</Integer>
  <String namespace="" id=""  articleName="SequenceString">$sequence</String>
</moby:$object_type>
PRT
	}

        push (@$moby_sequence_objects, $moby_object);
    }

    return $moby_sequence_objects;
}

sub parseSingleGFFIntoCollectionGFF {
    my $self = shift;
    my ($report, $output_format, $namespace) = @_;
    my $output_objects = [];
    
    my @lines = split ('\n', $report);
    
    # print STDERR "got " . @lines . "\n";
    # my $line_tmp = $lines[0];
    # print STDERR "first line, $line_tmp.\n";
    # $line_tmp = $lines[1];
    # print STDERR "second line, $line_tmp.\n";
    
    my $sequenceIdentifier;    
    my $report_tmp = "";
    
    while (my $line = shift @lines) {
	my $sequenceIdentifier_tmp;

	# Get the sequence identifier
	if ($line =~ /^([^\t]+)\t.+/) {
	    $sequenceIdentifier_tmp = $1;
	}
	else {
	    print STDERR "in parseSingleGFFIntoCollectionGFF, can't parse sequence identifier from GFF,\n$line\n";
	    exit 0;
	}
	
	if (not defined $sequenceIdentifier) {
	    $sequenceIdentifier = $sequenceIdentifier_tmp;
	    $report_tmp .= $line . "\n";
	}
	elsif ($sequenceIdentifier_tmp eq $sequenceIdentifier) {
	    $report_tmp .= $line . "\n";
	}
	else {

	    # New sequence report
	    # Build the GFF Moby object

	    my $input = <<PRT;
<moby:$output_format namespace='$namespace' id='$sequenceIdentifier'>
<String namespace='' id='' articleName='content'>
<![CDATA[
$report_tmp
]]>
</String>
</moby:$output_format>
PRT

            push (@$output_objects, $input);

	    # Reinitialisation
	    $sequenceIdentifier = $sequenceIdentifier_tmp;
	    $report_tmp         = $line . "\n";
	    
        }
    }

    # Add the last report !

    my $input = <<PRT;
<moby:$output_format namespace='$namespace' id='$sequenceIdentifier'>
<String namespace='' id='' articleName='content'>
<![CDATA[
$report_tmp
]]>
</String>
</moby:$output_format>
PRT

    push (@$output_objects, $input);
    
    return $output_objects;
    
}

# @param sequences in FASTA format as a string
# @return a set of bioperl objects

sub _parse_fasta_sequences {
    my ($fasta_sequences) = @_;
    
    my @fasta_sequence_lines = split (/\n/, $fasta_sequences);
    
    my $seqobjs = [];
    
    my $seq_id_previous;
    my $seq_desc_previous = "";
    my $seq_str = "";
    
    my $seq_id_current;
    my $seq_desc_current = "";

    foreach my $line (@fasta_sequence_lines) {
	
	if ($line =~ /^>/) {
	    
	    # It's a definition line

	    if ($line =~ /\s/) {
		$line =~ /^>([^\s]+)\s+(.*)/;
		
		$seq_id_current   = $1;
		$seq_desc_current = $2 || "";
	    }
	    else {
		$line =~ /^>(.+)/;
		
		$seq_id_current = $1;
	    }
	    
	    if (!defined $seq_id_current) {
		print STDERR "Problem parsing sequence identifier FASTA line,\n$line\n";
	    }

	    if (defined ($seq_id_previous)) {
		# Parsing new sequence entry
		# Store the previous sequence entry
		
		my $seqobj = Bio::Seq->new (
					    -display_id => $seq_id_previous,
					    -desc       => $seq_desc_previous,
					    -seq        => $seq_str,
					   );
		push (@$seqobjs, $seqobj);
	    }
	    
	    # (Re)initialisation
	    
	    $seq_id_previous   = $seq_id_current;
	    $seq_desc_previous = $seq_desc_current;
	    $seq_str  = "";
	}
	elsif (length ($line) > 0) {
	    
	    # It's a sequence line
	    
	    $seq_str .= $line;
	}
	else {
	    # empty line
	}
    }
    
    if (defined $seq_id_current) {
	# Store the last sequence entry
	
	my $seqobj = Bio::Seq->new (
				    -display_id => $seq_id_current,
				    -desc       => $seq_desc_current,
				    -seq        => $seq_str,
				   );

	push (@$seqobjs, $seqobj);
    }
    
    return $seqobjs;
}

# @Input $sequences_hash, a reference to a hashtable of sequences
# The hash keys are the sequence identifiers, and the hash values are the string sequences.
# @Input $DOM, The XML Moby Document containing the sequence data

# @Return %sequences_hash, the input Hashtable (as a hash not a reference). The parsed sequences will be appended to it.

sub parseMobySequenceObjectFromDOM {
  my $self = shift;
  my ($DOM, $sequences_hash) = @_;

  my $sequenceIdentifier;
  my $sequence_str;
  
  my @articles = ($DOM);
  ($sequenceIdentifier) = getSimpleArticleIDs (\@articles);
  if (not defined $sequenceIdentifier) {
    print STDERR "Error - no sequence identifier!!!\n";
    exit 0;
  }                                    
  
  # Los contenidos los devuelve como una lista, dado que
  # el objeto de la ontologia podria tener una relacion
  # "has" n-aria. Bien, en nuestro caso solo habia un peptido.

  # The Sequence as a string  

  ($sequence_str) = getNodeContentWithArticle($DOM, "String", "SequenceString");
  # Lo que hacemos aqui es limpiar un sting de caracteres raros
  # (espacios, \n, ...) pq nadie asegura que no los hayan.
  $sequence_str =~ s/\W+//sg; # trim trailing whitespace

  # Add the sequence into a hash table
   
  $sequences_hash->{$sequenceIdentifier} = $sequence_str;
  
  return %$sequences_hash;
}

# @Input $DOM, The XML Moby Document containing the sequence data

# @Return $seqobj, a Bioperl sequence object

sub parseMobySequenceObjectFromDOMintoBioperlObject {
  my $self = shift;
  my ($DOM) = @_;

  my $sequenceIdentifier;
  my $sequence_str;
  my $description;
  
  my @articles = ($DOM);
  ($sequenceIdentifier) = getSimpleArticleIDs (\@articles);
  if (not defined $sequenceIdentifier) {
    print STDERR "Error - no sequence identifier!!!\n";
    exit 0;
  }                                    
  
  # Los contenidos los devuelve como una lista, dado que
  # el objeto de la ontologia podria tener una relacion
  # "has" n-aria. Bien, en nuestro caso solo habia un peptido.

  # The Sequence as a string  

  ($sequence_str) = getNodeContentWithArticle($DOM, "String", "SequenceString");
  # Lo que hacemos aqui es limpiar un sting de caracteres raros
  # (espacios, \n, ...) pq nadie asegura que no los hayan.
  $sequence_str =~ s/\W+//sg; # trim trailing whitespace

  # The description
  ($description) = getNodeContentWithArticle($DOM, "String", "Description");
  if (not defined $description) {
      $description = "";
  }

  # Instanciate the bioperl object
   
  my $seqobj = Bio::Seq->new (
			      -display_id => $sequenceIdentifier,
			      -seq        => $sequence_str,
			      -desc       => $description,
			      );
  
  return $seqobj;
}

sub convertSequencesIntoFASTA {
    my ($seqobjs)  = @_;

    my ($fh, $tempfile) = tempfile("/tmp/FASTA.XXXXXX", UNLINK => 0);
    my $out = Bio::SeqIO->new(
			      -fh => $fh,
			      -format => 'fasta'
			      );

    foreach my $seqobj (@$seqobjs) {
	$out->write_seq ($seqobj);
    }

    close $fh;
    $fh = new FileHandle "$tempfile", "r";
    my @fasta_sequences = $fh->getlines;
    unlink $tempfile;

    my $fasta_sequences = join ('', @fasta_sequences);

    return $fasta_sequences;
}

# Check that the input datatype matched the specified one

# Be careful with this validation as the input datatype could well have the 'moby:' prefix

sub validateDataType {
    my $self = shift;
    my ($DOM, $specifiedType) = @_;

    my $inputDataType = undef;
    my $rightType     = undef;

    eval {
    
    # input 
    # * DOM containing articles we want to validate the type
    # * the specified expected type

    my @object_nodes = ();
    
    if ($DOM->nodeName =~ /collection/i) {
	my @simple_nodes = getCollectedSimples ($DOM);
	foreach my $simple_node (@simple_nodes) {
	    # Get the object node
	    my ($node) = $simple_node->getElementsByTagName ('*');
	    push (@object_nodes, $node);
	}
    }
    else {
	# it is a simple - get directly the object node from the DOM
	my ($node) = $DOM->getElementsByTagName ('*');
	push (@object_nodes, $node);
    }
    
    foreach my $node (@object_nodes) {

	# must all match the specified type
	
	my $nodeType = $node->nodeType();
	# should be already this type !!
	if ($node->nodeType() == ELEMENT_NODE) {
	    
	    # That should be the object - Validate code is here!!
	    
	    $inputDataType = $node->nodeName;
	    
	    if ($specifiedType eq "GenericSequence") {
		if ($inputDataType =~ /GenericSequence|AminoAcidSequence|AASequence|NucleotideSequence|DNASequence|RNASequence/) {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }
	    
	    if ($specifiedType eq "NucleotideSequence") {
		if ($inputDataType =~ /NucleotideSequence|DNASequence|RNASequence/) {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }
	    
	    if ($specifiedType =~ /DNASequence/) {
		if ($inputDataType =~ /DNASequence/) {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }
	    
	    if ($specifiedType =~ /RNASequence/) {
		if ($inputDataType =~ /RNASequence/) {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }
	    
	    if ($specifiedType eq "AminoAcidSequence") {
		if ($inputDataType =~ /AminoAcidSequence|AASequence/) {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }

	    if ($specifiedType =~ /AASequence/) {
		if ($inputDataType =~ /AASequence/) {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }

	    if ($specifiedType eq "FASTA") {
		if ($inputDataType =~ /^FASTA/) {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }
	    
	    if ($specifiedType eq "FASTA_NA") {
		if ($inputDataType eq "FASTA_NA") {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }
	    
	    if ($specifiedType eq "FASTA_NA_multi") {
		if ($inputDataType =~ /^FASTA_NA/) {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }
	    
	    if ($specifiedType eq "FASTA_AA") {
		if ($inputDataType eq "FASTA_AA") {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }
	    
	    if ($specifiedType eq "FASTA_AA_multi") {
		if ($inputDataType =~ /^FASTA_AA/) {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }

	    if (($specifiedType eq "BLAST-Text") || ($specifiedType eq "NCBI_BLAST_Text") || ($specifiedType eq "WU_BLAST_Text")) {
		if ($inputDataType =~ /BLAST-Text|NCBI_BLAST_Text|WU_BLAST_Text/) {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }
	    
	    if ($specifiedType eq "GFF2") {
		if ($inputDataType =~ /GFF2$/) {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }

	    if ($specifiedType eq "GFF3") {
		if ($inputDataType =~ /GFF3$/) {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }

	    if ($specifiedType eq "GFF") {
		if ($inputDataType =~ /GFF\d*$/) {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }
	    
	    if ($specifiedType eq "MatrixFloat") {
		if ($inputDataType =~ /MatrixFloat$/) {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }
	    
	    if ($specifiedType eq "MatrixInteger") {
		if ($inputDataType =~ /MatrixInteger$/) {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }

	    if ($specifiedType eq "Matrix") {
		if ($inputDataType =~ /Matrix$/) {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }

	    if ($specifiedType eq "Meta_Alignment_Text") {
		if ($inputDataType =~ /Meta_Alignment_Text$/) {
		    $rightType = 1;
		}
		else {
		    # Wrong input type
		    return (0, $inputDataType);
		}
	    }
	    
	    # ...
	    
	}
	else {
	    print STDERR "not an ELEMENT node...\n";
	}
    }
    
    };
    if ($@) {
	my $note = "Internal System Error. Can not parse input XML for validating type, $specifiedType!\n";
	print STDERR "$note\n";
	return (undef, undef);
    }

    return ($rightType, $inputDataType);
}

# Return an array of all namespaces associated with a simple or a collection
# If the collection contains simples with different namespace, then return all of them.

sub getNamespace {
    my $self = shift;
    my ($DOM) = @_;
    
    # input 
    # * DOM containing articles we want to validate the type
    
    my @namespaces = ();
    my @object_nodes = ();
    
    if ($DOM->nodeName =~ /collection/i) {
	my @simple_nodes = getCollectedSimples ($DOM);
	foreach my $simple_node (@simple_nodes) {
	    # Get the object node
	    my ($node) = $simple_node->getElementsByTagName ('*');
	    push (@object_nodes, $node);
	}
    }
    else {
	# it is a simple - get directly the object node from the DOM
	my ($node) = $DOM->getElementsByTagName ('*');
	push (@object_nodes, $node);
    }
    
    foreach my $node (@object_nodes) {
	# should be already this type !!
	if ($node->nodeType() == ELEMENT_NODE) {
	    my $namespace = $node->getAttributeNode ("namespace")->getValue();
	    if (! is_in (\@namespaces, $namespace)) {
		push (@namespaces, $namespace);
	    }
	}
    }
    
    return @namespaces;
}

sub is_in {
    my ($aref, $input_element) = @_;
    foreach my $element (@$aref) {
	if ($element eq $input_element) {
	    return 1;
	}
    }
    return 0;
}

# Set a Moby response that will return to a user
# It takes:
# * the output embedded in a mobyData block
# * a set of Exceptions
# * the moby_logger instance so we can log the executino status of the requested service
# * the service name
# it returns a mobyContent object

sub setMobyResponse {
    my $self = shift;
    my ($MOBY_RESPONSE, $moby_exceptions, $moby_logger, $serviceName) = @_;
    
    if (@$moby_exceptions > 0) {
	# build the moby exception response
	my $moby_exception_response = "";
	my %severities;
	foreach my $moby_exception (@$moby_exceptions) {
	    my $severity = $moby_exception->getExceptionType;
	    $severities{$severity} = $moby_exception;
	    $moby_exception_response .= $moby_exception->retrieveExceptionResponse() . "\n";
	}
	
	# logging report
	# Check 'error' first then 'warning' or 'information'
	if (defined $severities{error}) {
	    my $exception = $severities{error};
	    $moby_logger->info ("$serviceName failed");
	    $moby_logger->info ("Exception code, " . $exception->getExceptionCode);
	    $moby_logger->info ("Exception message, " . $exception->getExceptionMessage);
	}
	elsif (defined $severities{warning} || defined $severities{information}) {
	    my $exception = $severities{error};
	    $moby_logger->info ("$serviceName terminated successfully with warning or information notes");
	    $moby_logger->info ("Exception code, " . $exception->getExceptionCode);
	    $moby_logger->info ("Exception message, " . $exception->getExceptionMessage);
	}
	
	return responseHeader(
			      -authority => "genome.imim.es",
			      -note      => "$moby_exception_response"
			      )
	    . $MOBY_RESPONSE . responseFooter;
    }
    else {
	$moby_logger->info ("$serviceName terminated successfully");
	$moby_logger->info ("Exception code, 700");
	
	my $note = "Service execution succeeded";
	return responseHeader (
			       -authority => "genome.imim.es",
			       -note      => "<Notes>$note</Notes>"
			       )
	    . $MOBY_RESPONSE . responseFooter;
    }
}

# PWM syntax - tab-delimited
# Transfac or jaspar compliant format

# motif0
# 1    0.206132  0.759466  0.027833  0.006569
# 2    0.027765  0.591999  0.188669  0.191568
# 3    0.292155  0.060594  0.048565  0.598686
# 4    0.028134  0.779645  0.008042  0.184179
# 5    0.069995  0.019596  0.905345  0.005063
# 6    0.013693  0.973960  0.005145  0.007203
# 7    0.086467  0.207722  0.372804  0.333006
# 8    0.002733  0.015888  0.245232  0.736148
# 9    0.011034  0.873492  0.080148  0.035326
# 10    0.002281  0.983921  0.007792  0.006006
# 11    0.056512  0.865194  0.056446  0.021847
# 12    0.206273  0.375163  0.292359  0.126204
# //

# matrix_object_name  = ['Matrix', 'MatrixFloat', 'MatrixInteger', 'Distance_Matrix']
# matrix_element_type = ['String', 'Integer', 'Float']

sub convert_tabularPositionWeightMatrix_into_MobyMatrix {
    my $self = shift;
    my ($tab_matrix, $matrix_object_name, $matrix_element_type) = @_;

    # intermediary matrix object as an array of array
    my @matrix;
    # output moby matrix object
    my $moby_matrix_object;
    
    #############################################
    #
    # parsing the matrix in tab-delimited format
    # 
    #############################################
    
    my @lines = split ('\n', $tab_matrix);
    
    if (@lines == 0) {
	print STDERR "can't parse any lines...\n";
	exit 1;
    }

    # The first line is the identifier

    my $first_line = shift @lines;
    chomp $first_line;
    $first_line =~ /^([^\s]+)/;
    my $motif_identifier = $1;
    
    # print STDERR "motif identifier, $motif_identifier.\n";
    
    # the last line is a separator '//'

    pop @lines;

    # Now we are only left with a set of arrays of values, 
    # first column of each array is a label

    my $i = 0;
    foreach my $line (@lines) {
	
	chomp $line;
	if (($line =~ /\t$/) || ($line =~ /\s$/)) {
	    chop $line;
	}
	
	# print STDERR "parsing line,$line...\n";
	
	if (! ($line =~ /\t/)) {
	    print STDERR "Warning, it could be a problem, matrix doesn't have a tab-delimited syntax (\t).\n";
	}
	
	my @values =  split ('\s+', $line);
	
	if (@values == 0) {
	    print STDERR "can't parse any values...\n";
	    exit 1;
	}
	
	# Get rid of the label (first column)
	shift @values;
	
	my $j = 0;
	foreach my $value (@values) {
	    # print STDERR "parsing value,$value...\n";
	    chomp $value;

	    $matrix[$i][$j] = $value;
	    $j++;
	}

	$i++;
    }
    
    #####################################
    #
    # Writing out the moby matrix object
    #
    #####################################
    
    $moby_matrix_object = "<$matrix_object_name namespace='' id='$motif_identifier'>\n";
    $moby_matrix_object .= "<Integer articleName='Key'/>\n";
    
    for my $row (0..$#matrix) {
	$moby_matrix_object .= "<Array" . $matrix_element_type . " articleName='Array'>\n";
        $moby_matrix_object .= "<Integer articleName='Key'>$row</Integer>\n";
	
	for my $col (0..$#{$matrix[$row]}) {
	    my $value = $matrix[$row][$col];
	    
            $moby_matrix_object .= "<Element" . $matrix_element_type . " articleName='Element'>";
	    $moby_matrix_object .= "\<Integer articleName='Key'>$col</Integer>";
	    $moby_matrix_object .= "<" . $matrix_element_type . " articleName='Value'>" . $value . "</" . $matrix_element_type . ">";
	    $moby_matrix_object .= "</Element" . $matrix_element_type . ">\n";
	}
	
	$moby_matrix_object .= "</Array" . $matrix_element_type . ">\n";
	
    }
    
    $moby_matrix_object .= "</$matrix_object_name>\n";

    # print STDERR "debugging moby matrix object, $moby_matrix_object\n";

    return $moby_matrix_object;
}

# Score distance matrix format
# The first line is a commented line with the list of identifiers
# These identifiers are also in the first column
# So we will parse the identifiers from the first column and will not take into account the commented line at the top of the file !

# e.g.

#       CG9855  CG31911 CG1916  CG40178 CG5390  CG6936  CG5304  CG11331 CG5731
# CG9855  -       28.38   17.67   35.56   11.06   14.51   21.82   22.16   -5.02
# CG31911 28.38   -       41.74   32.83   12.47   34.92   38.27   34.93   -6.61
# CG1916  17.67   41.74   -       42.22   -2.50   17.58   37.41   35.19   -6.71
# CG40178 35.56   32.83   42.22   -       20.95   22.08   21.05   32.10   -7.65
# CG5390  11.06   12.47   -2.50   20.95   -       -2.31   -3.56   -0.16   5.25
# CG6936  14.51   34.92   17.58   22.08   -2.31   -       44.56   29.38   -7.30
# CG5304  21.82   38.27   37.41   21.05   -3.56   44.56   -       41.23   -7.93
# CG11331 22.16   34.93   35.19   32.10   -0.16   29.38   41.23   -       -6.09
# CG5731  -5.02   -6.61   -6.71   -7.65   5.25    -7.30   -7.93   -6.09   -

# matrix_object_name  = ['Matrix', 'MatrixFloat', 'MatrixInteger', 'Distance_Matrix']
# matrix_element_type = ['String', 'Integer', 'Float']

sub convert_tabularScoreMatrix_into_MobyMatrix {
    my $self = shift;
    my ($tab_matrix, $matrix_object_name, $matrix_element_type) = @_;
    
    my $debug = 0;
    
    # intermediary matrix object as an array of array
    my @matrix;
    my @labels = ();
    # output moby matrix object
    my $moby_matrix_object;
    
    #############################################
    #
    # parsing the matrix in tab-delimited format
    # 
    #############################################
    
    my @lines = split ('\n', $tab_matrix);
    
    if (@lines == 0) {
	print STDERR "can't parse any lines...\n";
	exit 1;
    }

    my $i = 0;
    foreach my $line (@lines) {
	
	chomp $line;
	
	if ($line =~ /^#/) {
	    # don't parse commented line
	    next;
	}
	
	if (($line =~ /\t$/) || ($line =~ /\s$/)) {
	    chop $line;
	}

	# print STDERR "parsing line,$line...\n";
	
	if (! ($line =~ /\t/)) {
	    print STDERR "Warning, it could be a problem, matrix doesn't have a tab-delimited syntax (\t).\n";
	}
	
	my @values =  split ('\s+', $line);
	
	if (@values == 0) {
	    print STDERR "can't parse any values...\n";
	    exit 1;
	}
	
	# Get the label (first column)
	my $identifier = shift @values;
	push (@labels, $identifier);
	
	my $j = 0;
	foreach my $value (@values) {
	
	    # print STDERR "parsing value,$value...\n";
	    
	    chomp $value;
	    
	    $matrix[$i][$j] = $value;
	    $j++;
	}
	
	$i++;
    }
    
    #####################################
    #
    # Writing out the moby matrix object
    #
    #####################################
    
    $moby_matrix_object = "<$matrix_object_name namespace='' id=''>\n";
    $moby_matrix_object .= "<Integer articleName='Key'/>\n";
    
    # The labels in an ArrayString
    
    $moby_matrix_object .= "<ArrayString articleName='Label'>\n";
    $moby_matrix_object .= "<Integer articleName='Key'>0</Integer>\n";
    
    for $i  (0..$#labels) {
	my $identifier = $labels[$i];

	$moby_matrix_object .= "<ElementString articleName='Element'>";
	$moby_matrix_object .= "<Integer articleName='Key'>$i</Integer>";
	$moby_matrix_object .= "<String articleName='Value'>$identifier</String>";
	$moby_matrix_object .= "</ElementString>\n";
    }
    $moby_matrix_object .= "</ArrayString>\n";

    for my $row (0..$#matrix) {
	$moby_matrix_object .= "<Array" . $matrix_element_type . " articleName='Array'>\n";
        $moby_matrix_object .= "<Integer articleName='Key'>$row</Integer>\n";
	
	for my $col (0..$#{$matrix[$row]}) {
	    my $value = $matrix[$row][$col];
	    
	    if ($debug) {
		print STDERR "value, $value\n";
	    }
	    
	    if ($value eq "-") {
		$value = undef;
	    }
	    
            $moby_matrix_object .= "<Element" . $matrix_element_type . " articleName='Element'>";
	    $moby_matrix_object .= "<Integer articleName='Key'>$col</Integer>";
	    if (defined $value) {
		$moby_matrix_object .= "<" . $matrix_element_type . " articleName='Value'>" . $value . "</" . $matrix_element_type . ">";
	    }
	    else {
		$moby_matrix_object .= "<" . $matrix_element_type . " articleName='Value'/>";
	    }
	    $moby_matrix_object .= "</Element" . $matrix_element_type . ">\n";
	}
	
	$moby_matrix_object .= "</Array" . $matrix_element_type . ">\n";
	
    }
    
    $moby_matrix_object .= "</$matrix_object_name>\n";
    
    # print STDERR "debugging moby matrix object, $moby_matrix_object\n";
    
    return $moby_matrix_object;
}


# Check if the matrix mode matches the type of the Matrix

sub convert_MobyMatrix_into_tabularMatrix {
    my $self = shift;
    my ($DOM, $matrix_type, $debug) = @_;
    
    # output
    my $matrix_text = "";
    my $moby_exceptions = [];
    
    my $matrixIdentifier;
    my $_delimitor = "   ";
    
    my @articles = ($DOM);
    ($matrixIdentifier) = getSimpleArticleIDs (\@articles);
    
    if (not defined $matrixIdentifier) {
	print STDERR "Error - no matrix identifier!!!\n";
	# return exception
	# ...
	return ("", $moby_exceptions);
    }
    
    $matrix_text = "$matrixIdentifier\n";
    
    if ($debug) {
	print STDERR "motif identifer, $matrixIdentifier\n";
    }
    
    # Parse the matrix object
    
    my @matrix = _parseMatrixMobyObject ($DOM, $matrix_type, $debug);
    
    if (@matrix == 0) {
	print STDERR "Empty matrix!!! - Error parsing the Matrix object\n";
	# return exception
	# ...
	return ("", $moby_exceptions);
    }
    
    # Build the matrix in text format from the matrix array reference - matrix[$row][$col]
    
    for my $row (0..$#matrix) {
	my $row_index = $row + 1;
	$matrix_text .= "$row_index" . "$_delimitor";
	for my $col (0..$#{$matrix[$row]}) {
	    $matrix_text .= "$_delimitor" if ($col);
	    $matrix_text .= $matrix[$row][$col];
	}
	$matrix_text .= "\n";
    }
    
    $matrix_text .= "//\n";

    if ($debug) {
	print STDERR "matrix text, $matrix_text\n";
    }
    
    return ($matrix_text, $moby_exceptions);
}

sub _parseMatrixMobyObject {
    my ($matrix_DOM, $matrix_type, $debug) = @_;
    my @matrix;
    
    unless ( ref( $matrix_DOM ) =~ /XML\:\:LibXML/ ) {
	my $parser  = XML::LibXML->new();
	my $doc     = $parser->parse_string( $matrix_DOM );
	$matrix_DOM = $doc->getDocumentElement();
    }
    
    # Be careful with the labels as they are also in ArrayString, it is confusing !!!

    my $array_elements = $matrix_DOM->getElementsByTagName ("Array" . $matrix_type);
    my $size = $array_elements->size();
    if ($size == 0) {
	$array_elements = $matrix_DOM->getElementsByTagName ("moby:Array" . $matrix_type);
	$size = $array_elements->size();
	if ($size == 0) {
	    print STDERR "Error, can't parse any array element from the Matrix moby XML...\n";
	    return ();
	}
    }
    
    my $i = 0;
    while ($i < $size) {
	my $array_element = $array_elements->get_node ($i);

	# Get the row index
	
	my $row_elements = $array_element->getElementsByTagName ("Integer");
	my $row_size = $row_elements->size();
	if ($row_size == 0) {
	    $row_elements = $array_element->getElementsByTagName ("moby:Integer");
	    $row_size = $row_elements->size();
	    if ($row_size == 0) {
		print STDERR "Error, can't parse any row index element from the Matrix moby XML...\n";
		return ();
	    }
	}
	my $row = $row_elements->[0]->textContent();
	
	if ($debug) {
	    print STDERR "row, $row\n";
	}
	
	# Get the elements in the array
	
	my $array_element_elements = $array_element->getElementsByTagName ("Element" . $matrix_type);
	my $element_size  = $array_element_elements->size();
	if ($element_size == 0) {
	    $array_element_elements = $array_element->getElementsByTagName ("moby:Element" . $matrix_type);
	    $element_size = $array_element_elements->size();
	    if ($element_size == 0) {
		print STDERR "Error, can't parse any array_element element from the Matrix moby XML...\n";
		return ();
	    }
	}
	
	if ($debug) {
	    print STDERR "element size, $element_size\n";
	}
	
	my $j = 0;
	while ($j < $element_size) {
	    my $array_element_element = $array_element_elements->[$j];
	
	    if ($debug) {
		print STDERR "array_element element dumping, " . $array_element_element->toString() . "\n";
	    }
	    
	    # The key
	    
	    my $key_elements = $array_element_element->getElementsByTagName ("Integer");
	    my $sub_element_size  = $key_elements->size();
	    if ($sub_element_size == 0) {
		$key_elements = $array_element_element->getElementsByTagName ("moby:Integer");
		$sub_element_size = $key_elements->size();
		if ($sub_element_size == 0) {
		    print STDERR "Error, can't parse any array_element element from the Matrix moby XML...\n";
		    return ();
		}
	    }
	    
	    if ($debug) {
		print STDERR "sub element size, $sub_element_size\n";
	    }
	    
	    my $key_element = $key_elements->[0];
	    my $col = $key_element->textContent();

	    # The value
	    
	    my $value_elements = $array_element_element->getElementsByTagName ($matrix_type);
	    $sub_element_size  = $key_elements->size();
	    if ($sub_element_size == 0) {
		$value_elements = $array_element_element->getElementsByTagName ("moby:" . $matrix_type);
		$sub_element_size = $value_elements->size();
		if ($sub_element_size == 0) {
		    print STDERR "Error, can't parse any array_element element from the Matrix moby XML...\n";
		    return ();
		}
	    }
	    
	    if ($debug) {
		print STDERR "sub element size, $sub_element_size\n";
	    }
	    
	    my $value_element = $value_elements->[0];
	    my $value = $value_element->textContent();
	    
	    if ($debug) {
		print STDERR "column,$col.\n";
		print STDERR "value,$value.\n";
	    }
	    
	    $matrix[$row][$col] = $value;

	    $j++;
	}
	
	$i++;
    }
    
    return @matrix;
}

1;

__END__
