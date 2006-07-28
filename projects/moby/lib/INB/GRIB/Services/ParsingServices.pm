# $Id: ParsingServices.pm,v 1.8 2006-07-28 13:47:03 gmaster Exp $
#
# This file is an instance of a template written 
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#
# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::ParsingServices  - Package for Parsing services

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

Copyright (c) 2005, Arnaud Kerhornou and INB - Nodo 1 - GRIB/IMIM.
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

package INB::GRIB::Services::ParsingServices;

use strict;
use warnings;
use Carp;

use INB::GRIB::Services::Factory;
use INB::GRIB::Utils::CommonUtilsSubs;
# Official CommonSubs
# that one uses DOM XML parsing which is not optimal for huge dataset
use MOBY::CommonSubs qw(:all);
# Seb Carrere CommonSubs which uses XSLT which is better, because faster !!!
use MOBY::MOBYXSLT;

# Logging
use Log::Log4perl qw(get_logger :levels);
use Data::Dumper;

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
  &fromMetaAlignmentstoScoreMatrix
  &fromMetaAlignmentstoTextScoreMatrix
  &parseMotifMatricesFromMEME
);

our $VERSION = '1.0';

my $_debug = 0;

# Preloaded methods go here.

###############################################################################

=head2 _do_query_generateScoreMatrixFromMetaAlignment

 Title   : _do_query_generateScoreMatrixFromMetaAlignment
         : 
         : private function (NOT EXPORTED)
         : 
 Usage   : my $query_response = _do_query_generateScoreMatrix($query);
         : 
         : donde:
         :   $query es un XML::DOM::Node que contiene el arbol xml que
         :   engloba:  
         :     <moby:mobyData queryID='1'>...</moby:mobyData>
         : 
 Returns : Devuelve un string que contiene el resultado de la ejecución
         : para una sola query.
         : Un ejemplo sería: 
         : 
         : <moby:mobyData moby:queryID='1'>
         :   <moby:Simple moby:articleName='report'>
         :     <moby:text-plain namespace='Global_Keyword' id=''>
	 :    ....
         :     </moby:text-plain>
         :   </moby:Simple>
         : </moby:mobyData>

=cut

sub _do_query_generateScoreMatrixFromMetaAlignment {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby 
    my $queryInput_DOM      = shift @_;
    my $output_format       = shift @_;
    my $output_moby_type    = shift @_;
    my $output_article_name = shift @_;
    
    my $MOBY_RESPONSE   = "";     # set empty response
    my $moby_exceptions = [];
    
    # Aqui escribimos las variables que necesitamos para la funcion.
    my $input_format;
    
    # Variables that will be passed to generateScoreMatrix_call
    
    my %parameters;
    my $input_data = [];
    my $queryID    = MOBY::MOBYXSLT::getInputID ($queryInput_DOM);
    my @articles   = MOBY::MOBYXSLT::getArticles($queryInput_DOM);
    my $namespace = "";

    # Get the parameters
    
    # no parameter !
    
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

	if (MOBY::MOBYXSLT::isSimpleArticle ($DOM)) {
	    # simple input is not allowed
	    
	    my $note = "Received a simple input article instead of a collection";
	    print STDERR "$note\n";
	    my $code = "201";
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      refElement => "similarity_results",
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
	
	# It's not very nice but taverna doesn't set up easily article name for input data so we let the users not setting up the article name of the input
	
	if ((MOBY::MOBYXSLT::isCollectionArticle ($DOM)) || (defined ($articleName) && ($articleName eq "similarity_results"))) {
	    
	    if ($_debug) {
		print STDERR "node ref, " . ref ($DOM) . "\n";
		print STDERR "DOM: " . $DOM->toString () . "\n";
	    }
	    
	    my @inputs_article_DOMs = MOBY::MOBYXSLT::getCollectedSimples ($DOM);
	    foreach my $input_DOM (@inputs_article_DOMs) {
		
		my $inputDataType = MOBY::MOBYXSLT::getObjectType ($input_DOM);
		my $rightType = 0;
		
		if ($inputDataType eq "Meta_Alignment_Text") {
		  $rightType++;
		}
		
		if (!$rightType) {
		    my $note = "Expecting a Meta_Alignment_Text object, and receiving a $inputDataType object";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => 'similarity_results',
									      code       => $code,
									      type       => 'error',
									      queryID    => $queryID,
									      message    => "$note",
									      );
		    push (@$moby_exceptions, $moby_exception);
		    
		    # Simple Response doesn't fit !! (the simple article is not empty as it should be!), so we need to create the string from scratch !
		    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_SIMPLE_RESPONSE ($queryID, $output_article_name);
		    return ($MOBY_RESPONSE, $moby_exceptions);
		}
		
		my $input = MOBY::MOBYXSLT::getObjectContent ($input_DOM);
		if ((not defined $input) || (length ($input) < 1)) {
		    # should be an exception!!
		    print STDERR "Error, can't parse any input!!\n";
		    exit 0;
		}
		
		if ($_debug) {
		    print STDERR "parsing content done!\n";
		    print STDERR "input, $input\n";
		}
		
		push (@$input_data, $input);
		
	    }
	}	
	
    } # Next article
    
    if (@$input_data < 1) {
	my $note = "can't parse any pairwise similarity data...\n";
	print STDERR "$note\n";
	my $code = "201";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "similarity_results",
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_SIMPLE_RESPONSE ($queryID, $output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    else {
	if ($_debug) {
	    print STDERR "\nparsed @$input_data Meta_Alignment_Text objects\n\n";
	}
    }
    
    my ($matrix_text, $moby_exceptions_tmp) = generateScoreMatrix_call (similarity_results  => $input_data, queryID => $queryID, parameters => \%parameters, debug => $_debug);
    push (@$moby_exceptions, @$moby_exceptions_tmp);
    
    # Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX 
    # nos queda encapsularla en un Objeto bioMoby. Esta operacio 
    # la podriamos realizar en una funcion a parte si fuese compleja.
    
    if ($_debug) {
	print STDERR "generated matrix:\n$matrix_text\n";
    }
    
    if (defined $matrix_text) {
	
	my $output_object;
	
	if ($output_format eq "xml") {
	    
	    if ($_debug) {
		print STDERR "returning the matrix in xml format\n";
	    }
	    
	    # Build the Moby Matrix object
	    $output_object = INB::GRIB::Utils::CommonUtilsSubs->convert_tabularScoreMatrix_into_MobyMatrix ($matrix_text, $output_moby_type, "Float");
	}
	else {
	    
	    if ($_debug) {
		print STDERR "returning the matrix in text format\n";
	    }
	    
	    $output_object = <<PRT;
<$output_moby_type namespace='$namespace' id=''>
<String namespace='' id='' articleName='content'>
<![CDATA[
$matrix_text
]]>
</String>
</$output_moby_type>
PRT

	}
	
	if ($_debug) {
	    print STDERR "Matrix moby object,\n";
	    print STDERR "$output_object\n";
	}
	
        $MOBY_RESPONSE = simpleResponse($output_object, $output_article_name, $queryID);
    }
    else {
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_SIMPLE_RESPONSE ($queryID, $output_article_name);
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


sub _do_query_MemeMotifMatrices {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby
    my $queryInput_DOM = shift @_;
    my $_output_format = shift @_;

    my $MOBY_RESPONSE   = "";     # set empty response
    my $moby_exceptions = [];
    my $output_article_name = "motif_weight_matrices";
    
    # Aqui escribimos las variables que necesitamos para la funcion.
    my $matrix_mode;

    # Variables that will be passed to meme2matrix_call
    my $meme_predictions;
    my %parameters;

    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);

    # Get the parameters

    ($matrix_mode) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "matrix mode");
    if (not defined $matrix_mode) {
	$matrix_mode = "log-likelihood";
    }
    elsif (! ((lc ($matrix_mode) eq "raw format") || (lc ($matrix_mode) eq "log-likelihood"))) {
	my $note = "matrix mode parameter, '$matrix_mode', not accepted, should be ['raw format','log-likelihood']";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "matrix mode",
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	# Return an empty moby data object, as well as an exception telling why nothing got returned
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_RESPONSE ($queryID, $output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    

    # Match the matrix object name and the type of its element depending on matrix mode, ie
    # raw format     => MatrixInteger
    # log-likelihood => MatrixFloat

    my $matrix_object_name;
    my $matrix_element_type;
    if ($matrix_mode eq "log-likelihood") {
	$matrix_object_name  = "MatrixFloat";
	$matrix_element_type = "Float";
    }
    elsif ($matrix_mode eq "raw format") {
	$matrix_object_name  = "MatrixInteger";
	$matrix_element_type = "Integer";
    }
    
    # Add the parsed parameters in a hash table
    
    $parameters{matrix_mode} = $matrix_mode;

    if ($_debug) {
	print STDERR "matrix mode, $matrix_mode\n";
    }

    # Tratamos a cada uno de los articulos
    foreach my $article (@articles) {

	# El articulo es una tupla que contiene el nombre de este
	# y su texto xml.

	my ($articleName, $DOM) = @{$article}; # get the named article

	# Si le hemos puesto nombre a los articulos del servicio,
	# podemos recoger a traves de estos nombres el valor.
	# Sino sabemos que es el input articulo porque es un simple/collection articulo

	# It's not very nice but taverna doesn't set up easily article name for input data so we let the users not setting up the article name of the input (which should be 'sequences')
	# In case of MatScan, it doesn't really matter as there is only one input anyway

	if (($articleName eq "meme_predictions") || (isSimpleArticle ($DOM))) {

	    if ($_debug) {
		print STDERR "meme_predictions...\n";
		print STDERR "Simple article DOM: " . $DOM->toString() . "\n";
	    }
	    
	    if (isCollectionArticle ($DOM)) {
		my $note = "Received a collection input article instead of a simple one";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => "meme_predictions",
									  code       => $code,
									  type       => 'error',
									  queryID    => $queryID,
									  message    => "$note",
									  );
		push (@$moby_exceptions, $moby_exception);
		
		# Return an empty moby data object, as well as an exception telling what nothing got returned
		
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_RESPONSE ($queryID, $output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
	    }
	    
	    # Validation input datatype
	    
	    my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "MEME_Text");
	    if (!$rightType) {
		my $note = "Expecting a MEME_Text object, and receiving a $inputDataType object";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => 'meme_predictions',
									  code       => $code,
									  type       => 'error',
									  queryID    => $queryID,
									  message    => "$note",
									  );
		push (@$moby_exceptions, $moby_exception);
		
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
	    }
	    
	    $meme_predictions = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "MEME_Text");

	    if ($_debug) {
		print STDERR "parsed meme_predictions,\n$meme_predictions.\n";
	    }

	}
	if (isCollectionArticle ($DOM)) {
	    my $note = "Received a collection input article instead of a simple one";
	    print STDERR "$note\n";
	    my $code = "201";
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      refElement => "meme_predictions",
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    push (@$moby_exceptions, $moby_exception);
	    
	    # Return an empty moby data object, as well as an exception telling what nothing got returned
	    
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
	    return ($MOBY_RESPONSE, $moby_exceptions);
	}
    }

    # Una vez recogido todos los parametros necesarios, llamamos a
    # la funcion que nos devuelve el report.

    # Return an array of matrices
    my ($matrices_aref, $moby_exceptions_tmp) = meme2matrix_call (meme_predictions => $meme_predictions, format => $_output_format, queryID => $queryID, parameters => \%parameters);
    push (@$moby_exceptions, @$moby_exceptions_tmp);
    
    # Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX
    # nos queda encapsularla en un Objeto bioMoby. Esta operacio
    # la podriamos realizar en una funcion a parte si fuese compleja.

    my $output_object_type  = "$_output_format";
    my $namespace = "";
    
    if ((!defined $matrices_aref) || (@$matrices_aref < 1)) {

	if ($_debug) {
	    print STDERR "no meme matrices, returning an empty collection object\n";
	}

	# Return an emtpy message !
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    my $meme_matrix_objects = [];    
    foreach my $matrix (@$matrices_aref) {
	my $moby_matrix = INB::GRIB::Utils::CommonUtilsSubs->convert_tabularPositionWeightMatrix_into_MobyMatrix ($matrix, $matrix_object_name, $matrix_element_type);
	push (@$meme_matrix_objects, $moby_matrix);
    }

    # Bien!!! ya tenemos el objeto de salida del servicio , solo nos queda
    # volver a encapsularlo en un objeto biomoby de respuesta. Pero
    # en este caso disponemos de una funcion que lo realiza. Si tuvieramos
    # una respuesta compleja (de verdad, esta era simple ;) llamariamos
    # a collection response.
    # IMPORTANTE: el identificador de la respuesta ($queryID) debe ser
    # el mismo que el de la query.
    
    $MOBY_RESPONSE = collectionResponse($meme_matrix_objects, $output_article_name, $queryID);
    return ($MOBY_RESPONSE, $moby_exceptions);
    
}

#################################
# Interface service methods
#################################

=head2 fromMetaAlignmentstoScoreMatrix

 Title   : fromMetaAlignmentstoScoreMatrix
 Usage   : Esta función está pensada para llamarla desde un cliente SOAP. 
         : No obstante, se recomienda probarla en la misma máquina, antes 
         : de instalar el servicio. Para ello, podemos llamarla de la 
         : siguiente forma:
         : 
         : my $result = GeneID("call", $in);
         : 
         : donde $in es texto que con el mensaje biomoby que contendría
         : la parte del <tag> "BODY" del mensaje soap. Es decir, un string
         : de la forma: 
         :  
         :  <?xml version='1.0' encoding='UTF-8'?>
         :   <moby:MOBY xmlns:moby='http://www.biomoby.org/moby-s'>
         :    <moby:mobyContent>
         :      <moby:mobyData queryID='1'> 
         :      ...
         :      </moby:mobyData>
         :    </moby:mobyContent>
         :   </moby:mobyContent>
         :  </moby:MOBY>
 Returns : Devuelve un string que contiene el resultado de todas las 
         : queries en GFF formate. Es decir, un mensaje xml de la forma:
         :
         : <?xml version='1.0' encoding='UTF-8'?>
         : <moby:MOBY xmlns:moby='http://www.biomoby.org/moby' 
         : xmlns='http://www.biomoby.org/moby'>
         :   <moby:mobyContent moby:authority='inb.lsi.upc.es'>
         :     <moby:mobyData moby:queryID='1'>
         :       ....
         :     </moby:mobyData>
         :     <moby:mobyData moby:queryID='2'>
         :       ....
         :     </moby:mobyData>
         :   </moby:mobyContent>
         :</moby:MOBY>

=cut

sub fromMetaAlignmentsToScoreMatrix {
    
    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_; # get the incoming MOBY query XML

    if ($_debug) {
	print STDERR "processing Moby fromMetaAlignmentsToScoreMatrix query...\n";
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
    my $moby_logger     = get_logger ("MobyServices");
    my $serviceName     = "fromMetaAlignmentsToScoreMatrix";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){
	
	# En este punto es importante recordar que el objeto $query 
	# es un XML::DOM::Node, y que si queremos trabajar con 
	# el mensaje de texto debemos llamar a: $query->toString() 
	
	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}
	
	my $output_format       = "xml";
	my $output_moby_type    = "Distance_Matrix";
	my $output_article_name = "matrix";
	
	my ($query_response, $moby_exceptions_tmp) = _do_query_generateScoreMatrixFromMetaAlignment ($queryInput, $output_format, $output_moby_type, $output_article_name);
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

=head2 fromMetaAlignmentstoTextScoreMatrix

 Title   : fromMetaAlignmentstoTextScoreMatrix
 Usage   : Esta función está pensada para llamarla desde un cliente SOAP. 
         : No obstante, se recomienda probarla en la misma máquina, antes 
         : de instalar el servicio. Para ello, podemos llamarla de la 
         : siguiente forma:
         : 
         : my $result = GeneID("call", $in);
         : 
         : donde $in es texto que con el mensaje biomoby que contendría
         : la parte del <tag> "BODY" del mensaje soap. Es decir, un string
         : de la forma: 
         :  
         :  <?xml version='1.0' encoding='UTF-8'?>
         :   <moby:MOBY xmlns:moby='http://www.biomoby.org/moby-s'>
         :    <moby:mobyContent>
         :      <moby:mobyData queryID='1'> 
         :      ...
         :      </moby:mobyData>
         :    </moby:mobyContent>
         :   </moby:mobyContent>
         :  </moby:MOBY>
 Returns : Devuelve un string que contiene el resultado de todas las 
         : queries en GFF formate. Es decir, un mensaje xml de la forma:
         :
         : <?xml version='1.0' encoding='UTF-8'?>
         : <moby:MOBY xmlns:moby='http://www.biomoby.org/moby' 
         : xmlns='http://www.biomoby.org/moby'>
         :   <moby:mobyContent moby:authority='inb.lsi.upc.es'>
         :     <moby:mobyData moby:queryID='1'>
         :       ....
         :     </moby:mobyData>
         :     <moby:mobyData moby:queryID='2'>
         :       ....
         :     </moby:mobyData>
         :   </moby:mobyContent>
         :</moby:MOBY>

=cut

sub fromMetaAlignmentsToTextScoreMatrix {
    
    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_; # get the incoming MOBY query XML

    if ($_debug) {
	print STDERR "processing Moby fromMetaAlignmentsToTextScoreMatrix query...\n";
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
    
    # Using DOM
    # my @queries = getInputs($message);  # returns XML::DOM nodes
    # Using XSLT
    my ($serviceNotes, $queries) = MOBY::MOBYXSLT::getInputs($message);  # returns XML::DOM nodes
    
    # 
    # Inicializamos la Respuesta a string vacio. Recordar que la respuesta
    # es una coleccion de respuestas a cada una de las consultas.
    my $MOBY_RESPONSE   = "";             # set empty response
    my $moby_exceptions = [];
    my $moby_logger     = get_logger ("MobyServices");
    my $serviceName     = "fromMetaAlignmentsToTextScoreMatrix";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@$queries){
	
	# En este punto es importante recordar que el objeto $query 
	# es un XML::DOM::Node, y que si queremos trabajar con 
	# el mensaje de texto debemos llamar a: $query->toString() 
	
	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}

	my $output_format       = "text";
	my $output_moby_type    = "MicroArrayData_Text";
	my $output_article_name = "microarraydata";
	
	my ($query_response, $moby_exceptions_tmp) = _do_query_generateScoreMatrixFromMetaAlignment ($queryInput, $output_format, $output_moby_type, $output_article_name);
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


=head2 parseMotifMatricesFromMEME

 Title   : parseMotifMatricesFromMEME
 Usage   : Esta función está pensada para llamarla desde un cliente SOAP.
	 : No obstante, se recomienda probarla en la misma máquina, antes
	 : de instalar el servicio. Para ello, podemos llamarla de la
	 : siguiente forma:
	 :
	 : my $result = MatScan("call", $in);
	 :
	 : donde $in es texto que con el mensaje biomoby que contendría
	 : la parte del <tag> "BODY" del mensaje soap. Es decir, un string
	 : de la forma:
	 :
	 :  <?xml version='1.0' encoding='UTF-8'?>
	 :   <moby:MOBY xmlns:moby='http://www.biomoby.org/moby-s'>
	 :    <moby:mobyContent>
	 :      <moby:mobyData queryID='1'>
	 :      ...
	 :      </moby:mobyData>
	 :    </moby:mobyContent>
	 :   </moby:mobyContent>
	 :  </moby:MOBY>
 Returns : Return a collection of text-formatted objects with MEME motif matrices

=cut

sub parseMotifMatricesFromMEME {

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
    
    my $_format = "MEME_Text";
    my $moby_logger = get_logger ("MobyServices");
    my $serviceName = "parseMotifMatricesFromMEME";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){

	# En este punto es importante recordar que el objeto $query
	# es un XML::DOM::Node, y que si queremos trabajar con
	# el mensaje de texto debemos llamar a: $query->toString()

	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}

	my ($query_response, $moby_exceptions_tmp) = _do_query_MemeMotifMatrices ($queryInput, $_format);
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

