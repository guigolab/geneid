# $Id: ClusteringServices.pm,v 1.6 2007-03-13 21:59:19 gmaster Exp $
#
# This file is an instance of a template written
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#
# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::ClusteringServices  - Implementation Package of gene clustering services.

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

package INB::GRIB::Services::ClusteringServices;

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
  &runKMeansClustering
  &runSOTAClustering
);

our $VERSION = '1.0';

my $_debug = 0;

# Preloaded methods go here.

###############################################################################

=head2 _do_query_KMeans

 Title   : _do_query_KMeans
	 :
	 : private function (NOT EXPORTED)
	 :
 Usage   : my $query_response = _do_query_KMeans($query);
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

sub _do_query_KMeans {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby
    my $queryInput_DOM      = shift @_;
    # The moby output format
    my $_moby_output_format = shift @_;
    
    # Output definition
    my $MOBY_RESPONSE   = "";     # set empty response
    my $moby_exceptions = [];
    my $output_article_name = "gene_clusters";
    
    # Aqui escribimos las variables que necesitamos para la funcion.
    my $gene_centering;
    my $iteration_number;
    my $cluster_number;
    my $sequenceIdentifier;
    
    # Variables that will be passed to KMeans_call
    my $gene_matrix;
    my %parameters;
    
    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);
    
    # Get the parameters
    
    ($gene_centering)   = getNodeContentWithArticle($queryInput_DOM, "Parameter", "gene centering");
    ($iteration_number) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "iterations clusters");
    ($cluster_number)   = getNodeContentWithArticle($queryInput_DOM, "Parameter", "clusters number");

    if (not defined $gene_centering) {
	$gene_centering = "None";
    }
    elsif (! (($gene_centering =~ /^none$/i) || ($gene_centering =~ /^k-means$/i) || ($gene_centering =~ /^k-medians$/i))) {
	my $note = "'gene centering' parameter, '$gene_centering', not accepted should be part of ['None', 'K-means', 'K-Medians']";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "gene centering",
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	# Return an empty moby data object, as well as an exception telling why nothing got returned
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    if (not defined $iteration_number) {
	$iteration_number = 200;
    }
    elsif (! ($iteration_number =~ /\d+/)) {
	my $note = "iterations number parameter, '$iteration_number', not accepted as it is not an integer";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "iterations number",
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	# Return an empty moby data object, as well as an exception telling why nothing got returned
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    if (not defined $cluster_number) {
	$cluster_number = 10;
    }
    elsif (! ($cluster_number =~ /\d+/)) {
	my $note = "clusters number parameter, '$cluster_number', not accepted as it is not an integer";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "clusters number",
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	# Return an empty moby data object, as well as an exception telling why nothing got returned
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    # Add the parsed parameters in a hash table
    
    if ($_debug) {
	print STDERR "gene centering, $gene_centering\n";
	print STDERR "number of iterations, $iteration_number\n";
	print STDERR "number of clusters, $cluster_number\n";
    }
    
    $parameters{gene_centering}   = $gene_centering;
    $parameters{iteration_number} = $iteration_number;
    $parameters{cluster_number}   = $cluster_number;
    
    $parameters{output_format} = $_moby_output_format;
    
    # Tratamos a cada uno de los articulos
    foreach my $article (@articles) {

	# El articulo es una tupla que contiene el nombre de este
	# y su texto xml.
	
	my ($articleName, $DOM) = @{$article}; # get the named article
	
	if ($_debug) {
	    print STDERR "processing article, $articleName...\n";
	}
	
	if (isCollectionArticle($DOM)) {
	    
	    # not allowed
	    
	    my $note = "Received a collection input article instead of a simple";
	    print STDERR "$note\n";
	    my $code = "201";
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      refElement => "gene_score_matrix",
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
	
	if (($articleName eq "gene_score_matrix") || isSimpleArticle ($DOM)) {
	    
	    if ($_debug) {
		print STDERR "node ref, " . ref ($DOM) . "\n";
		print STDERR "DOM: " . $DOM->toString () . "\n";
	    }
	    
	    # Validate the type first
		
	    my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "MicroArrayData_Text");
	    if (!$rightType) {
		my $note = "Expecting a MicroArrayData_Text object, and receiving a $inputDataType object";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => "gene_score_matrix",
									  code       => $code,
									  type       => 'error',
									  queryID    => $queryID,
									  message    => "$note",
									  );
		push (@$moby_exceptions, $moby_exception);
		
		# Simple Response doesn't fit !! (the simple article is not empty as it should be!), so we need to create the string from scratch !
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
	    }
	    
	    ($sequenceIdentifier) = getSimpleArticleIDs ( [ $DOM ] );
	    
	    if ($_debug) {
		print STDERR "parsed the following sequence identifier, $sequenceIdentifier\n";
	    }
	    
	    $gene_matrix = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "MicroArrayData_Text");
	    
	    if ($_debug) {
		print STDERR "gene matrix, $gene_matrix\n";
	    }
	    
	} # End parsing the gene matrix article tag
	
    } # Next article
    
    # Check that we have parsed properly the sequences and the predictions
    
    if ((not defined $gene_matrix) || (length $gene_matrix < 1)) {
	my $note = "could not parse any data in 'gene_score_matrix' article";
	print STDERR "$note\n";
	my $code = "201";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "gene_score_matrix",
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
    
    # Una vez recogido todos los parametros necesarios, llamamos a
    # la funcion que nos devuelve el report.
    
    my ($gene_clusters_aref, $moby_exceptions_tmp) = KMeans_call (gene_matrix  => $gene_matrix, queryID => $queryID, debug => $_debug, parameters => \%parameters);
    push (@$moby_exceptions, @$moby_exceptions_tmp);
    
    # Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX
    # nos queda encapsularla en un Objeto bioMoby. Esta operacio
    # la podriamos realizar en una funcion a parte si fuese compleja.
    
    if (defined $gene_clusters_aref) {
	
	if ($_debug) {
	    print STDERR "gene clusters output defined\n";
	}
	
	my $namespace = "";
	my $inputs = [];
	
	foreach my $gene_list (@$gene_clusters_aref) {
	    
	    # Build the Moby object
	    
	    my $input = <<PRT;
<moby:$_moby_output_format namespace='' id='$sequenceIdentifier'>
<String namespace='' id='' articleName='content'>
<![CDATA[
$gene_list
]]>
</String>
</moby:$_moby_output_format>
PRT

            push (@$inputs, $input);

        }
	
        $MOBY_RESPONSE = collectionResponse($inputs, $output_article_name, $queryID);
    }
    else {
	
	if ($_debug) {
	    print STDERR "no gene clusters defined\n";
	}
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
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

=head2 _do_query_SOTA

 Title   : _do_query_SOTA
	 :
	 : private function (NOT EXPORTED)
	 :
 Usage   : my $query_response = _do_query_KMeans($query);
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

sub _do_query_SOTA {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby
    my $queryInput_DOM      = shift @_;
    # The moby output format
    my $_moby_output_format = shift @_;
    
    # Output definition
    my $MOBY_RESPONSE   = "";     # set empty response
    my $moby_exceptions = [];
    my $output_article_name = "gene_clusters";
    
    # Aqui escribimos las variables que necesitamos para la funcion.
    my $distance;
    my $sequenceIdentifier;
    
    # Variables that will be passed to SOTA_call
    my $gene_matrix;
    my %parameters;
    
    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);
    
    # Get the parameters
    
    ($distance)   = getNodeContentWithArticle($queryInput_DOM, "Parameter", "distance");

    if (not defined $distance) {
	$distance = "euclidean";
    }
    elsif (! (($distance =~ /^euclidean$/i) || ($distance =~ /^square$/i) || ($distance =~ /^correlation$/i) || ($distance =~ /^offset$/i) || ($distance =~ /^spearman$/i) || ($distance =~ /^jackknife$/i))) {
	my $note = "'gene centering' parameter, '$distance', not accepted should be part of ['euclidean', 'square', 'correlation', 'offset', 'spearman', 'jackknife']";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "distance",
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	# Return an empty moby data object, as well as an exception telling why nothing got returned
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    # Add the parsed parameters in a hash table
    
    if ($_debug) {
	print STDERR "distance, $distance\n";
    }
    
    $parameters{distance}   = $distance;
    
    $parameters{output_format} = $_moby_output_format;
    
    # Tratamos a cada uno de los articulos
    foreach my $article (@articles) {

	# El articulo es una tupla que contiene el nombre de este
	# y su texto xml.
	
	my ($articleName, $DOM) = @{$article}; # get the named article
	
	if ($_debug) {
	    print STDERR "processing article, $articleName...\n";
	}
	
	if (isCollectionArticle($DOM)) {
	    
	    # not allowed
	    
	    my $note = "Received a collection input article instead of a simple";
	    print STDERR "$note\n";
	    my $code = "201";
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      refElement => "gene_score_matrix",
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
	
	if (($articleName eq "gene_score_matrix") || isSimpleArticle ($DOM)) {
	    
	    if ($_debug) {
		print STDERR "node ref, " . ref ($DOM) . "\n";
		print STDERR "DOM: " . $DOM->toString () . "\n";
	    }
	    
	    # Validate the type first
		
	    my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "MicroArrayData_Text");
	    if (!$rightType) {
		my $note = "Expecting a MicroArrayData_Text object, and receiving a $inputDataType object";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => "gene_score_matrix",
									  code       => $code,
									  type       => 'error',
									  queryID    => $queryID,
									  message    => "$note",
									  );
		push (@$moby_exceptions, $moby_exception);
		
		# Simple Response doesn't fit !! (the simple article is not empty as it should be!), so we need to create the string from scratch !
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
	    }
	    
	    ($sequenceIdentifier) = getSimpleArticleIDs ( [ $DOM ] );
	    
	    if ($_debug) {
		print STDERR "parsed the following sequence identifier, $sequenceIdentifier\n";
	    }
	    
	    $gene_matrix = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "MicroArrayData_Text");
	    
	    if ($_debug) {
		print STDERR "gene matrix, $gene_matrix\n";
	    }
	    
	} # End parsing the gene matrix article tag
	
    } # Next article
    
    # Check that we have parsed properly the sequences and the predictions
    
    if ((not defined $gene_matrix) || (length $gene_matrix < 1)) {
	my $note = "could not parse any data in 'gene_score_matrix' article";
	print STDERR "$note\n";
	my $code = "201";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "gene_score_matrix",
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
    
    # Una vez recogido todos los parametros necesarios, llamamos a
    # la funcion que nos devuelve el report.
    
    my ($gene_clusters_aref, $gene_tree, $moby_exceptions_tmp) = SOTA_call (gene_matrix  => $gene_matrix, queryID => $queryID, debug => $_debug, parameters => \%parameters);
    push (@$moby_exceptions, @$moby_exceptions_tmp);
    
    # Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX
    # nos queda encapsularla en un Objeto bioMoby. Esta operacio
    # la podriamos realizar en una funcion a parte si fuese compleja.
    
    if (defined $gene_clusters_aref) {
	
	if ($_debug) {
	    print STDERR "gene clusters output defined\n";
	}
	
	my $namespace = "";
	my $inputs = [];
	
	foreach my $gene_list (@$gene_clusters_aref) {
	    
	    # Build the Moby object
	    
	    my $input = <<PRT;
<moby:$_moby_output_format namespace='' id='$sequenceIdentifier'>
<String namespace='' id='' articleName='content'>
<![CDATA[
$gene_list
]]>
</String>
</moby:$_moby_output_format>
PRT

            push (@$inputs, $input);

        }
	
        $MOBY_RESPONSE = collectionResponse($inputs, $output_article_name, $queryID);
    }
    else {
	
	if ($_debug) {
	    print STDERR "no gene clusters defined\n";
	}
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
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


##############################################################################

sub runKMeansClustering {
    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_; # get the incoming MOBY query XML
    
    my $_output_format = "List_Text";
    my $moby_logger    = get_logger ("MobyServices");
    my $serviceName    = "runKMeansClustering";
    
    if ($_debug) {
	print STDERR "processing Moby runKMeansClustering query...\n";
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

	my ($query_response, $moby_exceptions_tmp) = _do_query_KMeans ($queryInput, $_output_format);
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


sub runSOTAClustering {
    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_; # get the incoming MOBY query XML
    
    my $_output_format = "List_Text";
    my $moby_logger    = get_logger ("MobyServices");
    my $serviceName    = "runSOTAClustering";
    
    if ($_debug) {
	print STDERR "processing Moby runSOTAClustering query...\n";
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

	my ($query_response, $moby_exceptions_tmp) = _do_query_SOTA ($queryInput, $_output_format);
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
