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
);

our $VERSION = '1.0';

# Common subs

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
    return $element->textContent;
}


sub createSequenceObjectsFromFASTA {
    my $self = shift;
    my ($fasta_sequences, $object_type, $namespace) = @_;

    my $moby_sequence_objects = [];
    my $sequence_objects      = _parse_fasta_sequences ($fasta_sequences);

    # Parsing the Bioperl objects returned by _parse_fasta_sequences method

    foreach my $seqobj (@$sequence_objects) {

	my $seq_id     = $seqobj->display_id;
	my $sequence   = $seqobj->seq;
	my $seq_length = length ($sequence);
	
	my $moby_object = <<PRT;
<moby:$object_type namespace='$namespace' id='$seq_id'>
  <Integer namespace="" id="" articleName="Length">$seq_length</Integer>
  <String namespace="" id=""  articleName="SequenceString">$sequence</String>
</moby:$object_type>
PRT

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
<![CDATA[
$report_tmp
]]>
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
<![CDATA[
$report_tmp
]]>
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
	else {
	    # It's a sequence line

	    $seq_str .= $line;
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

1;

__END__
