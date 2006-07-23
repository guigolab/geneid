# $Id: GFF2PSServices.pm,v 1.2 2006-07-23 20:55:59 gmaster Exp $
#
# This file is an instance of a template written
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#
# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::GFF2PSServices  - Implementation Package of MetaAlignment services.

=head1 SYNOPSIS

Con este package podremos parsear las llamadas xml de BioMoby para
poder llamar al servicio runGeneid. Una vez que tengamos la salida de llamando
a las funciones de Factory.pm, podremos encapsularla a objectos BioMoby.

  # Esta llamada/s nos devuelve una variable que contiene el texto con la
  # salida del programa Geneid encapsulado en objetos Moby.

=head1 DESCRIPTION

Este package sirve para parsear las llamadas BioMoby de entrada y salida.
De esta forma hace de puente entre el cliente que llama el servicio y la
aplicacion geneid.

Tipicamente la libreria que necesitamos en este módulo es la MOBY::CommonSubs,
que podremos ver tecleando en el terminal:
> perldoc MOBY::CommonSubs

No obstante, en esta libreria puede que no encontremos todo lo que necesitamos.
Si es así tendiamos que utilizar un parser de XML y una libreria que lo
interprete. Estas librerias podrian ser:

use XML::DOM;
use XML::Parser

y hacer cosas como:

my $parser = new XML::DOM::Parser;
my $doc = $parser->parse($message);
my $nodes = $doc->getElementsByTagName ("moby:Simple");
my $xml   = $nodes->item(0)->getFirstChild->toString;

Atencion: recomiendo mirarse esta libreria, solo cuando CommonSubs no nos de
nuestra solucion. La razón de ello es que, si cambia el estandard bioMoby, en
principio las llamadas a CommonSubs no tendrían que cambiar, pero si no las
hacemos servir, puede que debamos modificar el código.


=head1 AUTHOR

Arnaud Kerhornou, akerhornou@imim.es

=head1 COPYRIGHT

Copyright (c) 2006, Arnaud Kerhornou and INB - Nodo 1 - GRIB/CRG.
 All Rights Reserved.

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.


=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package INB::GRIB::Services::GFF2PSServices;

use strict;
use warnings;
use Carp;

use INB::GRIB::Services::Factory;
use MOBY::CommonSubs qw(:all);

# Moby Exceptions
use INB::Exceptions::MobyException;

# Logging
use Log::Log4perl qw(get_logger :levels);
use Data::Dumper;

# b64 image Encoding
use MIME::Base64;

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

# Aqui pondremos las funciones que reciben el mensaje BioMoby. Ejemplo
#
# @EXPORT = qw( &func1 &func2);
#
our @EXPORT = qw(
  &runGFF2JPEG
);

our $VERSION = '1.0';

my $_debug = 0;

# Preloaded methods go here.

###############################################################################

sub _do_query_GFF2JPEG {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby
    my $queryInput_DOM = shift @_;
    # The moby output format
    my $moby_output_format = shift @_;

    # Output definition
    my $moby_exceptions = [];
    my $MOBY_RESPONSE   = ""; # set empty response
    my $output_article_name = "image";
    
    # Variables that will be passed to GFF2JPEG_call
    
    my %parameters;
    my $maps_gff = [];
    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);
    
    my $title = "Transcription factor binding site maps";
    
    $parameters{output_format} = $moby_output_format;
    $parameters{title} = $title;
    
    # Tratamos a cada uno de los articulos
    foreach my $article (@articles) {

	# El articulo es una tupla que contiene el nombre de este
	# y su texto xml.

	my ($articleName, $DOM) = @{$article}; # get the named article

	if ($_debug) {
	    print STDERR "processing article, $articleName...\n";
	}

	# Si le hemos puesto nombre a los articulos del servicio,
	# podemos recoger a traves de estos nombres el valor.
	# Sino sabemos que es el input articulo porque es un simple/collection articulo

	# It's not very nice but taverna doesn't set up easily article name for input data so we let the users not setting up the article name of the input (which should be 'sequences')
	# In case of GeneID, it doesn't really matter as there is only one input anyway

	if (isSimpleArticle($DOM)) {
	    my $note = "Received a simple input article instead of a collection";
	    print STDERR "$note\n";
	    my $code = "201";
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      refElement => "maps",
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    push (@$moby_exceptions, $moby_exception);
	    
	    # Return an empty moby data object, as well as an exception telling what nothing got returned
	    
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_SIMPLE_RESPONSE ($queryID, $output_article_name);
	    return ($MOBY_RESPONSE, $moby_exceptions);
	}

	if ((isCollectionArticle ($DOM)) || (defined ($articleName) && ($articleName eq "maps"))) {
	    
	    if ($_debug) {
		print STDERR "node ref, " . ref ($DOM) . "\n";
		print STDERR "DOM: " . $DOM->toString () . "\n";
	    }
	    
	    # Validate the type of the simples in the collection - should all be GFF objects
	    my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "GFF");
	    if (!$rightType) {
		my $note = "Expecting a set of GFF objects, and receiving a $inputDataType object";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => 'maps',
									  code       => $code,
									  type       => 'error',
									  queryID    => $queryID,
									  message    => "$note",
									  );
		push (@$moby_exceptions, $moby_exception);
		
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_SIMPLE_RESPONSE ($queryID, $output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
	    }
	    
	    my @maps_article_DOMs = getCollectedSimples ($DOM);
	    foreach my $map_DOM (@maps_article_DOMs) {
		
		my $map = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($map_DOM, "GFF");
		
		if ($_debug) {
		    print STDERR "map, $map\n";
		}
		
		push (@$maps_gff, $map);
		
	    }
	}
	
    } # Next article
    
    if ($_debug) {
	print STDERR "parsed " . @$maps_gff . " maps.\n";
    }
    
    if (@$maps_gff < 1) {
	my $note = "Can't parse any input maps!\n";
	print STDERR "$note\n";
	my $code = "201";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "maps",
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	# Return an empty moby data object, as well as an exception telling why nothing got returned
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_SIMPLE_RESPONSE ($queryID, $output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    else {
	
	if ($_debug) {
	    print STDERR "calling GFF2JPEG_call...\n";
	}
	
	my ($gff2jpeg_report, $moby_exceptions_tmp) = GFF2JPEG_call (maps => $maps_gff, debug => $_debug, queryID => $queryID, parameters => \%parameters);
	push (@$moby_exceptions, @$moby_exceptions_tmp);
	
	if ($_debug) {
	    print STDERR "GFF2JPEG_call call done.\n";
	}
	
	# Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX
	# nos queda encapsularla en un Objeto bioMoby. Esta operacio
	# la podriamos realizar en una funcion a parte si fuese compleja.
	
	if (defined $gff2jpeg_report) {
	    
	    # Encode the JPEG image
	    
	    my $gff2jpeg_encoded_report = encode_base64 ($gff2jpeg_report);
	    
	    my $namespace = "";
	    
	    # Build the Moby object
	    
	    # Image_Encoded (Inab)
	    
	    my $output_object = <<PRT;
<moby:$moby_output_format namespace='' id=''>
<String namespace='' id='' articleName='rawdata'>
<![CDATA[
$gff2jpeg_encoded_report
]]>
</String>
<String namespace='' id='' articleName='mimeType'>
<![CDATA[
image/jpeg
]]>
</String>
</moby:$moby_output_format>
PRT

            # b64_encoded_jpeg (icapture)
	    
	    $output_object = <<PRT;
<moby:$moby_output_format namespace='' id=''>
<String namespace='' id='' articleName='content'>
<![CDATA[
$gff2jpeg_encoded_report
]]>
</String>
</moby:$moby_output_format>
PRT

            $MOBY_RESPONSE = simpleResponse($output_object, $output_article_name, $queryID);
	    
        }
	else {
	    my $note = "gff2jpeg call failed\n";
	    print STDERR "$note\n";
	    my $code = "701";
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    push (@$moby_exceptions, $moby_exception);
	    
	    # Return an empty moby data object, as well as an exception telling what nothing got returned
	    
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_SIMPLE_RESPONSE ($queryID, $output_article_name);
	    return ($MOBY_RESPONSE, $moby_exceptions);
	}
    }
    
    # Bien!!! ya tenemos el objeto de salida del servicio , solo nos queda
    # volver a encapsularlo en un objeto biomoby de respuesta. Pero
    # en este caso disponemos de una funcion que lo realiza. Si tuvieramos
    # una respuesta compleja (de verdad, esta era simple ;) llamariamos
    # a collection response.
    # IMPORTANTE: el identificador de la respuesta ($queryID) debe ser
    # el mismo que el de la query.

    return ($MOBY_RESPONSE, $moby_exceptions);
}

###########################################################################################

sub runGFF2JPEG {
    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_; # get the incoming MOBY query XML
    
    my $_output_format = "Image_Encoded";
    my $moby_logger = get_logger ("MobyServices");
    my $serviceName = "runGFF2JPEG";
    
    if ($_debug) {
	print STDERR "processing Moby runGFF2JPEG query...\n";
    }
    
    # Hasta el momento, no existen objetos Perl de BioMoby paralelos
    # a la ontologia, y debemos contentarnos con trabajar directamente
    # con objetos DOM. Por consiguiente lo primero es recolectar la
    # lista de peticiones (queries) que tiene la peticion.
    #
    # En una misma llamada podemos tener mas de una peticion, y de
    # cada servicio depende la forma de trabajar con ellas. En este
    # caso las trataremos una a una, pero podriamos hacer Threads para
    # tratarlas en paralelo, podemos ver si se pueden aprovechar resultados
    # etc..
    my @queries = getInputs($message);  # returns XML::DOM nodes
    #
    # Inicializamos la Respuesta a string vacio. Recordar que la respuesta
    # es una coleccion de respuestas a cada una de las consultas.
    my $MOBY_RESPONSE   = "";             # set empty response
    my $moby_exceptions = [];
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){
	
	# En este punto es importante recordar que el objeto $query
	# es un XML::DOM::Node, y que si queremos trabajar con
	# el mensaje de texto debemos llamar a: $query->toString()
	
	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}
	
	my ($query_response, $moby_exceptions_tmp) = _do_query_GFF2JPEG ($queryInput, $_output_format);
	push (@$moby_exceptions, @$moby_exceptions_tmp);
	
	# $query_response es un string que contiene el codigo xml de
	# la respuesta.  Puesto que es un codigo bien formado, podemos
	# encadenar sin problemas una respuesta con otra.
	$MOBY_RESPONSE .= $query_response;
    }
    # Una vez tenemos la coleccion de respuestas, debemos encapsularlas
    # todas ellas con una cabecera y un final. Esto lo podemos hacer
    # con las llamadas de la libreria Common de BioMoby.
    my $response = INB::GRIB::Utils::CommonUtilsSubs->setMobyResponse ($MOBY_RESPONSE, $moby_exceptions, $moby_logger, $serviceName);
    
    return $response;
}

1;

__END__
