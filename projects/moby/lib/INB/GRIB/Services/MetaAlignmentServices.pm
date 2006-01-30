# $Id: MetaAlignmentServices.pm,v 1.11 2006-01-30 18:10:55 gmaster Exp $
#
# This file is an instance of a template written
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#
# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::MetaAlignmentServices  - Implementation Package of MetaAlignment services.

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

package INB::GRIB::Services::MetaAlignmentServices;

use strict;
use warnings;
use Carp;

use INB::GRIB::Services::Factory;
use MOBY::CommonSubs qw(:all);

# Moby Exceptions
use INB::Exceptions::MobyException;

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
  &runMetaAlignment
  &runMetaAlignmentGFF
  &runMultiMetaAlignment
  &runMultiMetaAlignmentGFF
);

our $VERSION = '1.0';

my $_debug = 0;

# Preloaded methods go here.

###############################################################################

=head2 _do_query_MetaAlignment

 Title   : _do_query_MetaAlignment
	 :
	 : private function (NOT EXPORTED)
	 :
 Usage   : my $query_response = _do_query_MetaAlignment($query);
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

sub _do_query_MetaAlignment {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby
    my $queryInput_DOM = shift @_;
    # The moby output format
    my $_moby_output_format = shift @_;

    # Output definition
    my $MOBY_RESPONSE   = "";     # set empty response
    my $moby_exceptions = [];

    # Aqui escribimos las variables que necesitamos para la funcion.
    my $alpha_penalty;
    my $lambda_penalty;
    my $mu_penalty;
    my $sequenceIdentifier_1;
    my $sequenceIdentifier_2;

    # Variables that will be passed to MetaAlignment_call
    my $map1;
    my $map2;
    my %parameters;

    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);

    # Get the parameters

    ($alpha_penalty)  = getNodeContentWithArticle($queryInput_DOM, "Parameter", "alpha penalty");
    ($lambda_penalty) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "lambda penalty");
    ($mu_penalty)     = getNodeContentWithArticle($queryInput_DOM, "Parameter", "mu penalty");

    if (not defined $alpha_penalty) {
	$alpha_penalty = 0.5;
    }
    if (not defined $lambda_penalty) {
	$lambda_penalty = 0.1;
    }
    if (not defined $mu_penalty) {
	$mu_penalty = 0.1;
    }

    # Add the parsed parameters in a hash table

    if ($_debug) {
	print STDERR "alpha penalty, $alpha_penalty\n";
	print STDERR "lambda penalty, $lambda_penalty\n";
	print STDERR "mu penalty, $mu_penalty\n";
    }

    $parameters{alpha_penalty}  = $alpha_penalty;
    $parameters{lambda_penalty} = $lambda_penalty;
    $parameters{mu_penalty}     = $mu_penalty;

    $parameters{output_format} = $_moby_output_format;

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

	if ($articleName eq "map1") {

	    if ($_debug) {
		print STDERR "node ref, " . ref ($DOM) . "\n";
		print STDERR "DOM: " . $DOM->toString () . "\n";
	    }

	    ($sequenceIdentifier_1) = getSimpleArticleIDs ( [ $DOM ] );

	    if ((not defined $sequenceIdentifier_1) || (length ($sequenceIdentifier_1) == 0)) {
		print STDERR "Error, can not parsed the sequence identifier the GFF (map1) is attach to!\n";
		exit 0;
	    }

	    if ($_debug) {
		print STDERR "parsed the following sequence identifier for map1, $sequenceIdentifier_1\n";
	    }

	    $map1 = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "GFF");

	    if ($_debug) {
		print STDERR "map1, $map1\n";
	    }

	} # End parsing map1 article tag

	if ($articleName eq "map2") {

	    if ($_debug) {
		print STDERR "node ref, " . ref ($DOM) . "\n";
		print STDERR "DOM: " . $DOM->toString () . "\n";
	    }

	    ($sequenceIdentifier_2) = getSimpleArticleIDs ( [ $DOM ] );

	    if ((not defined $sequenceIdentifier_2) || (length ($sequenceIdentifier_2) == 0)) {
		print STDERR "Error, can not parsed the sequence identifier the GFF (map2) is attach to!\n";
		exit 0;
	    }

	    if ($_debug) {
		print STDERR "parsed the following sequence identifier for map2, $sequenceIdentifier_2\n";
	    }

	    $map2 = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "GFF");

	    if ($_debug) {
		print STDERR "map2, $map2\n";
	    }

	} # End parsing map2 article tag

    } # Next article

    # Check that we have parsed properly the sequences and the predictions

    if ((not defined $map1) || (not defined $map2)) {
	print STDERR "Error, can't parsed any maps...\n";
    }

    # Una vez recogido todos los parametros necesarios, llamamos a
    # la funcion que nos devuelve el report.

    my ($meta_report, $moby_exceptions_tmp) = MetaAlignment_call (map1  => $map1, map2  => $map2, queryID => $queryID, parameters => \%parameters);
    push (@$moby_exceptions, @$moby_exceptions_tmp);
    
    # Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX
    # nos queda encapsularla en un Objeto bioMoby. Esta operacio
    # la podriamos realizar en una funcion a parte si fuese compleja.

    my $output_article_name = "meta_predictions";
    my $namespace = "";

    # Build the Moby object

    my $input = <<PRT;
<moby:$_moby_output_format namespace='' id='$sequenceIdentifier_1'>
<![CDATA[
$meta_report
]]>
</moby:$_moby_output_format>
PRT

    # Bien!!! ya tenemos el objeto de salida del servicio , solo nos queda
    # volver a encapsularlo en un objeto biomoby de respuesta. Pero
    # en este caso disponemos de una funcion que lo realiza. Si tuvieramos
    # una respuesta compleja (de verdad, esta era simple ;) llamariamos
    # a collection response.
    # IMPORTANTE: el identificador de la respuesta ($queryID) debe ser
    # el mismo que el de la query.

    $MOBY_RESPONSE .= simpleResponse($input, $output_article_name, $queryID);

    return ($MOBY_RESPONSE, $moby_exceptions);
}


=head2 _do_query_MultiMetaAlignment

 Title   : _do_query_MultiMetaAlignment
	 :
	 : private function (NOT EXPORTED)
	 :
 Usage   : my $query_response = _do_query_MultiMetaAlignment($query);
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

sub _do_query_MultiMetaAlignment {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby
    my $queryInput_DOM = shift @_;
    # The moby output format
    my $_moby_output_format = shift @_;

    # Output definition
    my $moby_exceptions = [];
    my $MOBY_RESPONSE   = "";     # set empty response

    # Aqui escribimos las variables que necesitamos para la funcion.
    my $alpha_penalty;
    my $lambda_penalty;
    my $mu_penalty;

    # Variables that will be passed to MetaAlignment_call

    my %parameters;
    my $maps_gff = [];
    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);

    # Get the parameters

    ($alpha_penalty)  = getNodeContentWithArticle($queryInput_DOM, "Parameter", "alpha penalty");
    ($lambda_penalty) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "lambda penalty");
    ($mu_penalty)     = getNodeContentWithArticle($queryInput_DOM, "Parameter", "mu penalty");

    if (not defined $alpha_penalty) {
	$alpha_penalty = 0.5;
    }
    if (not defined $lambda_penalty) {
	$lambda_penalty = 0.1;
    }
    if (not defined $mu_penalty) {
	$mu_penalty = 0.1;
    }

    # Add the parsed parameters in a hash table

    if ($_debug) {
	print STDERR "alpha penalty, $alpha_penalty\n";
	print STDERR "lambda penalty, $lambda_penalty\n";
	print STDERR "mu penalty, $mu_penalty\n";
    }

    $parameters{alpha_penalty}  = $alpha_penalty;
    $parameters{lambda_penalty} = $lambda_penalty;
    $parameters{mu_penalty}     = $mu_penalty;

    $parameters{output_format} = $_moby_output_format;

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

	if ((isCollectionArticle ($DOM)) || (defined ($articleName) && ($articleName eq "maps"))) {

	    if (isSimpleArticle($DOM)) {
		print STDERR "problem, input article is simple - should be a collection!!!\n";
		exit 0;
	    }

	    if ($_debug) {
		print STDERR "node ref, " . ref ($DOM) . "\n";
		print STDERR "DOM: " . $DOM->toString () . "\n";
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

    # Make GFF pairs

    if ($_debug) {
	print STDERR "parsed " . @$maps_gff . " maps.\n";
    }

    my $output_objects = [];

    if (@$maps_gff < 2) {
	print STDERR "Error, less than two maps parsed from the input, can not run meta-alignment!!\n";
	exit 0;
    }
    else {

	# Una vez recogido todos los parametros necesarios, llamamos a
	# la funcion que nos devuelve el report.

	my $allowSelf    = 0;
	my $allowInverse = 0;

	for (my $i = 0; $i < @$maps_gff; $i++) {
	    for (my $j = 0; $j < @$maps_gff; $j++) {
		if (($i != $j || $allowSelf) && ($i <= $j || $allowInverse)) {

		    my $map1 = $maps_gff->[$i];
		    my $map2 = $maps_gff->[$j];

		    if (not defined $map1) {
			print STDERR "pb, map1 not defined !!\n";
			print STDERR "i: $i\n";
			print STDERR "j: $j\n";
		    }

		    if (not defined $map2) {
			print STDERR "pb, map2 not defined !!\n";
			print STDERR "i: $i\n";
			print STDERR "j: $j\n";
		    }

		    my ($meta_report, $moby_exceptions_tmp) = MetaAlignment_call (map1  => $map1, map2  => $map2, queryID => $queryID, parameters => \%parameters);
		    push (@$moby_exceptions, @$moby_exceptions_tmp);

		    # Leave it for now on, because for such collection, i don't know how to report properly exceptions

		    # Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX
		    # nos queda encapsularla en un Objeto bioMoby. Esta operacio
		    # la podriamos realizar en una funcion a parte si fuese compleja.

		    my $namespace = "";

		    # Build the Moby object

		    my $output_object = <<PRT;
<moby:$_moby_output_format namespace='' id=''>
<![CDATA[
$meta_report
]]>
</moby:$_moby_output_format>
PRT

		    push (@$output_objects, $output_object);
		}
	    }
	}
    }

    # Bien!!! ya tenemos el objeto de salida del servicio , solo nos queda
    # volver a encapsularlo en un objeto biomoby de respuesta. Pero
    # en este caso disponemos de una funcion que lo realiza. Si tuvieramos
    # una respuesta compleja (de verdad, esta era simple ;) llamariamos
    # a collection response.
    # IMPORTANTE: el identificador de la respuesta ($queryID) debe ser
    # el mismo que el de la query.

    my $output_article_name = "meta_predictions";
    $MOBY_RESPONSE .= collectionResponse($output_objects, $output_article_name, $queryID);

    return ($MOBY_RESPONSE, $moby_exceptions);
}


=head2 runMetaAlignment

 Title   : runMetaAlignment
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

sub runMetaAlignment {

    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_;        # get the incoming MOBY query XML

    my $_output_format = "text-formatted";

    if ($_debug) {
	print STDERR "processing Moby runMetaAlignment query...\n";
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

	my ($query_response, $moby_exceptions_tmp) = _do_query_MetaAlignment ($queryInput, $_output_format);
	push (@$moby_exceptions, @$moby_exceptions_tmp);

	# $query_response es un string que contiene el codigo xml de
	# la respuesta.  Puesto que es un codigo bien formado, podemos
	# encadenar sin problemas una respuesta con otra.
	$MOBY_RESPONSE .= $query_response;
    }
    # Una vez tenemos la coleccion de respuestas, debemos encapsularlas
    # todas ellas con una cabecera y un final. Esto lo podemos hacer
    # con las llamadas de la libreria Common de BioMoby.
    if (@$moby_exceptions > 0) {
	# build the moby exception response
	my $moby_exception_response = "";
	foreach my $moby_exception (@$moby_exceptions) {
	    $moby_exception_response .= $moby_exception->retrieveExceptionResponse() . "\n";
	}

	return responseHeader(
			       -authority => "genome.imim.es",
			       -note      => "$moby_exception_response"
			       )
	    . $MOBY_RESPONSE . responseFooter;
    }
    else {
	my $note = "Service execution succeeded";
	return responseHeader (
			       -authority => "genome.imim.es",
			       -note      => "<Notes>$note</Notes>"
			       )
	    . $MOBY_RESPONSE . responseFooter;
    }
}


=head2 runMetaAlignmentGFF

 Title   : runMetaAlignmentGFF
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

sub runMetaAlignmentGFF {

    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_;        # get the incoming MOBY query XML

    my $_output_format  = "GFF";

    if ($_debug) {
	print STDERR "processing Moby runMetaAlignmentGFF query...\n";
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

	my ($query_response, $moby_exceptions_tmp) = _do_query_MetaAlignment ($queryInput, $_output_format);
	push (@$moby_exceptions, @$moby_exceptions_tmp);

	# $query_response es un string que contiene el codigo xml de
	# la respuesta.  Puesto que es un codigo bien formado, podemos
	# encadenar sin problemas una respuesta con otra.
	$MOBY_RESPONSE .= $query_response;
    }
    # Una vez tenemos la coleccion de respuestas, debemos encapsularlas
    # todas ellas con una cabecera y un final. Esto lo podemos hacer
    # con las llamadas de la libreria Common de BioMoby.
    if (@$moby_exceptions > 0) {
	# build the moby exception response
	my $moby_exception_response = "";
	foreach my $moby_exception (@$moby_exceptions) {
	    $moby_exception_response .= $moby_exception->retrieveExceptionResponse() . "\n";
	}

	return responseHeader(
			       -authority => "genome.imim.es",
			       -note      => "$moby_exception_response"
			       )
	    . $MOBY_RESPONSE . responseFooter;
    }
    else {
	my $note = "Service execution succeeded";
	return responseHeader (
			       -authority => "genome.imim.es",
			       -note      => "<Notes>$note</Notes>"
			       )
	    . $MOBY_RESPONSE . responseFooter;
    }
}

=head2 runMultiMetaAlignment

 Title   : runMultiMetaAlignment
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

sub runMultiMetaAlignment {

    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_;        # get the incoming MOBY query XML

    my $_output_format = "text-formatted";

    if ($_debug) {
	print STDERR "processing Moby runMultiMetaAlignment query...\n";
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

	my ($query_response, $moby_exceptions_tmp) = _do_query_MultiMetaAlignment ($queryInput, $_output_format);
	push (@$moby_exceptions, @$moby_exceptions_tmp);

	# $query_response es un string que contiene el codigo xml de
	# la respuesta.  Puesto que es un codigo bien formado, podemos
	# encadenar sin problemas una respuesta con otra.
	$MOBY_RESPONSE .= $query_response;
    }
    # Una vez tenemos la coleccion de respuestas, debemos encapsularlas
    # todas ellas con una cabecera y un final. Esto lo podemos hacer
    # con las llamadas de la libreria Common de BioMoby.
    if (@$moby_exceptions > 0) {
	# build the moby exception response
	my $moby_exception_response = "";
	foreach my $moby_exception (@$moby_exceptions) {
	    $moby_exception_response .= $moby_exception->retrieveExceptionResponse() . "\n";
	}

	return responseHeader(
			      -authority => "genome.imim.es",
			      -note      => $moby_exception_response
			      )
	    . $MOBY_RESPONSE . responseFooter;
    }
    else {
	my $note = "Service execution succeeded";
	return responseHeader (
			       -authority => "genome.imim.es",
			       -note      => "<Notes>$note</Notes>"
			      )
	    . $MOBY_RESPONSE . responseFooter;
    }
}


=head2 runMultiMetaAlignmentGFF

 Title   : runMultiMetaAlignmentGFF
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

sub runMultiMetaAlignmentGFF {

    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_;        # get the incoming MOBY query XML

    my $_output_format = "GFF";

    if ($_debug) {
	print STDERR "processing Moby runMultiMetaAlignmentGFF query...\n";
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
    my $MOBY_RESPONSE = "";             # set empty response
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

	my ($query_response, $moby_exceptions_tmp) = _do_query_MultiMetaAlignment ($queryInput, $_output_format);
	push (@$moby_exceptions, @$moby_exceptions_tmp);

	# $query_response es un string que contiene el codigo xml de
	# la respuesta.  Puesto que es un codigo bien formado, podemos
	# encadenar sin problemas una respuesta con otra.
	$MOBY_RESPONSE .= $query_response;
    }
    # Una vez tenemos la coleccion de respuestas, debemos encapsularlas
    # todas ellas con una cabecera y un final. Esto lo podemos hacer
    # con las llamadas de la libreria Common de BioMoby.
    if (@$moby_exceptions > 0) {
	# build the moby exception response
	my $moby_exception_response = "";
	foreach my $moby_exception (@$moby_exceptions) {
	    $moby_exception_response .= $moby_exception->retrieveExceptionResponse() . "\n";
	}

	return responseHeader(
			      -authority => "genome.imim.es",
			      -note      => $moby_exception_response
			      )
	    . $MOBY_RESPONSE . responseFooter;
    }
    else {
	my $note = "Service execution succeeded";
	return responseHeader (
			       -authority => "genome.imim.es",
			       -note      => "<Notes>$note</Notes>"
			      )
	    . $MOBY_RESPONSE . responseFooter;
    }
}

1;

__END__
