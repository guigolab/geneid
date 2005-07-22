package INB::GRIB::Utils::CommonUtilsSubs;

use strict;
use warnings;
use Carp;
use XML::LibXML;

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

1;

__END__
