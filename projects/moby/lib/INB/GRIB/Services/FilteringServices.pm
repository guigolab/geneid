# $Id: FilteringServices.pm,v 1.1 2006-06-08 14:28:15 gmaster Exp $
#
# This file is an instance of a template written
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#
# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::FilteringServices  - Package for parser the Moby message to call the Sequence Filtering program .

=head1 SYNOPSIS

Con este package podremos parsear las llamadas xml de BioMoby para
poder llamar al servicio runFilteringGFF. Una vez que tengamos la salida de llamando
a las funciones de Factory.pm, podremos encapsularla a objectos BioMoby.

  #
  # En este módulo parsearemos los mensajes de BioMoby para realizar la
  # llamada al programa geneid.
  # Este modulo requiere las librerias de BioMoby. Para ver su
  # funcionamiento en modo linea de comandos (antes de activar el servicio)
  # enviaremos directamente un mensaje <xml> en BioMoby. Por ejemplo:
  $in = <<EOF
 <?xml version='1.0' encoding='UTF-8'?>
 <moby:MOBY xmlns:moby='http://www.biomoby.org/moby-s'>
 <moby:mobyContent>
   <moby:mobyData queryID='1'>
     <moby:Simple moby:articleName='seq'>
     <moby:DNASequence namespace='Global_Keyword' id=''>
       <moby:Integer namespace="" id="" articleName="Length">126</moby:Integer>
       <moby:String namespace="" id=""  articleName="SequenceString">
    ACTGCATGCTAAAGGTACATGACCGATCGGACTGTGACTGAACCTGCATTGA
    </moby:String>
     </moby:DNASequence>
     </moby:Simple>
     <moby:Parameter moby:articleName='arg2'><Value>swissprot</Value>
     </moby:Parameter>
   </moby:mobyData>
 </moby:mobyContent>
 </moby:MOBY>
 EOF
  my $result = MatScan("call", $in);


  # Esta llamada/s nos devuelve una variable que contiene el texto con la
  # salida del programa Geneid encapsulado en objetos Moby.

=head1 DESCRIPTION

Este package sirve para parsear las llamadas BioMoby de entrada y salida.
De esta forma hace de puente entre el cliente que llama el servicio y la
aplicacion MatScan.

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

Copyright (c) 2006, Arnaud Kerhornou and INB - Nodo INB 1 - GRIB/CRG.
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

package INB::GRIB::Services::FilteringServices;

use strict;
use warnings;
use Carp;

use INB::GRIB::Services::Factory;
use INB::GRIB::Utils::CommonUtilsSubs;
use MOBY::CommonSubs qw(:all);

# Logging
use Log::Log4perl qw(get_logger :levels);
use Data::Dumper;

# MIME encoding/decoding
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
  &filterSequencesByLength
  &filterSequencesAndQualityDataByLength
);

our $VERSION = '1.0';

my $_debug = 0;

# Preloaded methods go here.

###############################################################################

=head2 _do_query_sequence_filtering

 Title   : _do_query_sequence_filtering
	 :
	 : private function (NOT EXPORTED)
	 :
=cut

sub _do_query_sequence_filtering {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby
    my $queryInput_DOM   = shift @_;
    my $sequences_input_format     = shift @_;
    my $sequences_output_format    = shift @_;
    my $base_quality_output_format = shift @_;
    my $base_quality_input_format  = shift @_;
    # ServiceName tells us whether we receive quality data or just DNA sequences 
    my $serviceName = shift @_;
    
    if ($_debug) {
      print STDERR "Sequence filtering...\n";
    }
    
    my $MOBY_RESPONSE   = "";
    my $moby_exceptions = [];
    my $sequences_output_article_name    = "filtered_sequences";
    my $quality_data_output_article_name = "filtered_base_quality_data";
    my $namespace = "";
    
    # Aqui escribimos las variables que necesitamos para la funcion.
    my $trim_masked_regions;
    my $length_cutoff;
    
    # Variables that will be passed to Phrap_call
    my $fasta_quality_data_str;
    my $fasta_sequences_str;
    my %parameters;
    
    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);
    
    # Get the parameters
    
    ($trim_masked_regions) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "trim_masked_regions");
    if (not defined $trim_masked_regions) {
	$trim_masked_regions = 0;
    }
    elsif (! ($trim_masked_regions =~ /on|off/i)) {
	my $note = "trim_masked_regions parameter, '$trim_masked_regions', not valid, should be ['On', 'Off']";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "trim_masked_regions",
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	# Return an empty moby data object, as well as an exception telling what nothing got returned
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    # Map it into a boolean
    $trim_masked_regions =~ /^on$/i  and $trim_masked_regions = 1;
    $trim_masked_regions =~ /^off$/i and $trim_masked_regions = 0;
    
    ($length_cutoff) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "length_cutoff");
    if (not defined $length_cutoff) {
	$length_cutoff = 75;
    }
    elsif (! ($length_cutoff =~ /\d+/)) {
	my $note = "length_cutoff parameter, '$length_cutoff', not numerical";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "length_cutoff",
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	# Return an empty moby data object, as well as an exception telling what nothing got returned
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    # Add the parsed parameters in a hash table
    
    $parameters{trim_masked_regions} = $trim_masked_regions;
    $parameters{length_cutoff}       = $length_cutoff;
    
    if ($_debug) {
	print STDERR "trim_masked_regions, $trim_masked_regions\n";
	print STDERR "length_cutoff, $length_cutoff\n";
    }
    
    # Free some memory !
    undef $queryInput_DOM;
    
    # Tratamos a cada uno de los articulos
    foreach my $article (@articles) {
	
	# El articulo es una tupla que contiene el nombre de este
	# y su texto xml.
	
	my ($articleName, $DOM) = @{$article}; # get the named article
	
	# Si le hemos puesto nombre a los articulos del servicio,
	# podemos recoger a traves de estos nombres el valor.
	# Sino sabemos que es el input articulo porque es un simple/collection articulo
	
	if ($articleName =~ /sequences/i) {
	    
	    if (isCollectionArticle ($DOM)) {
		
		if ($_debug) {
		    print STDERR "$articleName tag is a collection article...\n";
		}
		
		# not allowed
		
		my $note = "Received a collection input article instead of a simple";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => 'sequences',
									  code       => $code,
									  type       => 'error',
									  queryID    => $queryID,
									  message    => "$note",
									  );
		push (@$moby_exceptions, $moby_exception);
		
		# Return an empty moby data object, as well as an exception telling what nothing got returned
		
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
	    }
	    
	    if ($_debug) {
		print STDERR "processing the sequences article...\n";
	    }
	    
	    # Validate the type of the simple article
	    my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "FASTA_NA_multi");
	    if (!$rightType) {
		my $note = "Expecting a FASTA_NA_multi object, and receiving a $inputDataType object";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => 'sequences',
									  code       => $code,
									  type       => 'error',
									  queryID    => $queryID,
									  message    => "$note",
									  );
		push (@$moby_exceptions, $moby_exception);
		
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
	    }
	    
	    if ($_debug) {
		print STDERR "Parsing the fasta sequences string...\n";
	    }
	    
	    $fasta_sequences_str  = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "String");
	    # Testing compression
	    # $fasta_sequences = decode_base64 ($fasta_sequences);
	    	    
	} # End parsing sequences article tag
	
	# Parse the quality data if any - check the serviceName, if not "filterSequencesAndQualityDataByLength" - tell that it is the wrong service !!
	
	if ($articleName =~ /base_quality_data/i) {
	    
	    # Chech also that the service is not "filterSequencesByLength", because not allowed with that one !!
	    
	    if ($serviceName eq "filterSequencesAndQualityDataByLength") {
		
		if ($_debug) {
		    print STDERR "processing the quality data article\n";
		}
		
		if (isCollectionArticle ($DOM)) {
		    
		    if ($_debug) {
			print STDERR "$articleName tag is a collection article - not allowed...\n";
		    }
		    
		    # not allowed
		    
		    my $note = "Received a collection input article instead of a simple";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => 'base_quality_data',
									      code       => $code,
									      type       => 'error',
									      queryID    => $queryID,
									      message    => "$note",
									      );
		    push (@$moby_exceptions, $moby_exception);
		    
		    # Return an empty moby data object, as well as an exception telling what nothing got returned
		    
		    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name);
		    return ($MOBY_RESPONSE, $moby_exceptions);
		}
		
		if ($_debug) {
		    print STDERR "$articleName is a simple article...\n";
		    print STDERR "Simple DOM: " . $DOM->toString() . "\n";
		}
		
		# Validate the type of the object
		my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "FASTA_Base_Quality_multi");
		if (!$rightType) {
		    my $note = "Expecting a FASTA_Base_Quality_multi object, and receiving a $inputDataType object";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => 'base_quality_data',
									      code       => $code,
									      type       => 'error',
									      queryID    => $queryID,
									      message    => "$note",
									      );
		    push (@$moby_exceptions, $moby_exception);
		    
		    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name);
		    return ($MOBY_RESPONSE, $moby_exceptions);
		}
		
		# Parsing
		# String instead of FASTA_Base_Quality_multi
		
		$fasta_quality_data_str = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "String");
		
	    }
	    else {
		my $note = "Found an article with quality data, but the service, $serviceName, doesn't require them, so they won't be taken into account. Use filterSequencesAndQualityDataByLength instead.";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => 'base_quality_data',
									  code       => $code,
									  type       => 'warning',
									  queryID    => $queryID,
									  message    => "$note",
									  );
		push (@$moby_exceptions, $moby_exception);
		
	    }
	} # End parsing base quality data
	
    } # Next article
    
    # Check that we have parsed properly the sequences
    
    # Free some memory!
    undef @articles;
    
    if (length ($fasta_sequences_str) < 1) {
	my $note = "can't parse any sequences...\n";
	print STDERR "$note\n";
	my $code = "201";
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name);
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "sequences",
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    if ($serviceName eq "filterSequencesAndQualityDataByLength") {
	if (length ($fasta_quality_data_str) < 1) {
	    my $note = "can't parse any quality data...\n";
	    print STDERR "$note\n";
	    my $code = "201";
	    
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name);
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      refElement => "base_quality_data",
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    push (@$moby_exceptions, $moby_exception);
	    
	    return ($MOBY_RESPONSE, $moby_exceptions);
	}
	
	#################################################################
	
	# Check also they have the same number of sequences
	
	# to do !!!
	
	my $nb_sequences = -1;
	my $nb_quality_sequences = -1;
	if ($nb_sequences != $nb_quality_sequences) {
	    my $note = "The number of DNA sequences, $nb_sequences, doesn't match the number of quality data sequences, $nb_quality_sequences...\n";
	    print STDERR "$note\n";
	    my $code = "201";
	    
	    $MOBY_RESPONSE     = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name);
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code    => $code,
								      type    => 'error',
								      queryID => $queryID,
								      message => "$note",
								      );
	    push (@$moby_exceptions, $moby_exception);
	    
	    return ($MOBY_RESPONSE, $moby_exceptions);
	}
	
	###################################################################
    }
	
    # Una vez recogido todos los parametros necesarios, llamamos a
    # la funcion que nos devuelve el report.
    
    if ($_debug) {
	print STDERR "parsing of the Moby input data done, preparing now the sequence filtering script execution...\n";
    }
    
    my ($filtered_seqs_fasta, $filtered_quality_data, $moby_exceptions_tmp) = SequenceFilteringByLength_call (sequences  => $fasta_sequences_str, quality_data => $fasta_quality_data_str, parameters => \%parameters, queryID => $queryID, debug => $_debug);
    push (@$moby_exceptions, @$moby_exceptions_tmp);
    
    undef $fasta_sequences_str;
    undef $fasta_quality_data_str;
    
    if (defined $filtered_seqs_fasta) {
	my $sequences_moby_seqobj = "<$sequences_output_format namespace='$namespace' id='Default'>\n<String id='' namespace='' articleName='content'><![CDATA[$filtered_seqs_fasta]]></String>\n</$sequences_output_format>\n";
	if ($serviceName eq "filterSequencesAndQualityDataByLength") {
	    if (defined $filtered_quality_data) {
		my $quality_data_moby_seqobj = "<$base_quality_output_format namespace='$namespace' id='Default'>\n<String id='' namespace='' articleName='content'><![CDATA[$filtered_quality_data]]></String>\n</$base_quality_output_format>\n";
		
		$MOBY_RESPONSE  = INB::GRIB::Utils::CommonUtilsSubs->MOBY_DOUBLE_SIMPLE_RESPONSE ($sequences_moby_seqobj, $sequences_output_article_name, $quality_data_moby_seqobj, $quality_data_output_article_name, $queryID);
	    }
	    else {
		if ($_debug) {
		    print STDERR "filtered output sequences fine, but not the filtered quality data!\n";
		}
		
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $quality_data_output_article_name);
	    }
	}
	else {
	    $MOBY_RESPONSE  = simpleResponse ($sequences_moby_seqobj, $sequences_output_article_name, $queryID);
	}
    }
    else {
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name);
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

##############################################

# Service interfaces

sub filterSequencesByLength {
    
    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_; # get the incoming MOBY query XML
    
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
    
    my $_sequences_input_format  = "FASTA_NA_multi";
    my $_sequences_output_format = "FASTA_NA_multi";
    my $moby_logger       = get_logger ("MobyServices");
    my $serviceName       = "filterSequencesByLength";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){
	
	# En este punto es importante recordar que el objeto $query
	# es un XML::DOM::Node, y que si queremos trabajar con
	# el mensaje de texto debemos llamar a: $query->toString()
	
	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}

	my ($query_response, $moby_exceptions_tmp) = _do_query_sequence_filtering ($queryInput, $_sequences_input_format, $_sequences_output_format, undef, undef, $serviceName);
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


sub filterSequencesAndQualityDataByLength {

    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_;        # get the incoming MOBY query XML

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
    
    my $_sequences_input_format     = "FASTA_NA_multi";
    my $_sequences_output_format    = "FASTA_NA_multi";
    my $_base_quality_input_format  = "FASTA_Base_Quality_multi";
    my $_base_quality_output_format = "FASTA_Base_Quality_multi";
    my $moby_logger       = get_logger ("MobyServices");
    my $serviceName       = "filterSequencesAndQualityDataByLength";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){
	
	# En este punto es importante recordar que el objeto $query
	# es un XML::DOM::Node, y que si queremos trabajar con
	# el mensaje de texto debemos llamar a: $query->toString()
	
	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}
	
	my ($query_response, $moby_exceptions_tmp) = _do_query_sequence_filtering ($queryInput, $_sequences_input_format, $_sequences_output_format, $_base_quality_input_format, $_base_quality_output_format, $serviceName);
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
