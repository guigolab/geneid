package INB::GRIB::Utils::CommonUtilsSubs;

use strict;
use warnings;
use Carp;
use Data::Dumper;

use XML::LibXML;

use Bio::Seq;

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
    my ($fasta_sequences, $object_type) = @_;

    my $moby_sequence_objects = [];
    my $sequence_objects      = _parse_fasta_sequences ($fasta_sequences);

    # Parsing the Bioperl objects returned by _parse_fasta_sequences method

    foreach my $seqobj (@$sequence_objects) {

	my $seq_id     = $seqobj->display_id;
	my $sequence   = $seqobj->seq;
	my $seq_length = length ($sequence);
	
	my $moby_object = <<PRT;
<moby:$object_type namespace='' id='$seq_id'>
  <Integer namespace="" id="" articleName="Length">$seq_length</Integer>
  <String namespace="" id=""  articleName="SequenceString">$sequence</String>
</moby:$object_type>
PRT

        push (@$moby_sequence_objects, $moby_object);
    }

    return $moby_sequence_objects;
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

1;

__END__
