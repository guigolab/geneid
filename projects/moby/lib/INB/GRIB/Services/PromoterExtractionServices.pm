# $Id: PromoterExtractionServices.pm,v 1.15 2006-07-22 14:41:33 gmaster Exp $
#
#
# This file is an instance of a template written 
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#
# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::PromoterExtractionServices  - Package for parser the Moby message to call the PromoterExtraction program .

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

Copyright (c) 2005, Arnaud Kerhornou and INB - Nodo Vertical 1 GRIB/IMIM.
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

package INB::GRIB::Services::PromoterExtractionServices;

use strict;
use warnings;
use Carp;

# INB
use INB::GRIB::Services::Factory; 
use INB::GRIB::Utils::CommonUtilsSubs;

use MOBY::CommonSubs qw(:all);

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
  &getUpstreamSeqFromEnsembl
  &getOrthologousUpstreamSeqFromEnsembl
);

our $VERSION = '1.0';

my $_debug = 0;

# Preloaded methods go here.

###############################################################################

=head2 _do_query_PromoterExtraction

 Title   : _do_query_PromoterExtraction
         : 
         : private function (NOT EXPORTED)
         : 
 Usage   : my $query_response = _do_query_PromoterExtraction($query);
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

sub _do_query_PromoterExtraction {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby 
    my $queryInput_DOM   = shift @_;
    my $orthologous_mode = shift @_;
    
    # The latest that is working !
    my $dbrelease = 38;
    
    my $MOBY_RESPONSE   = ""; # set empty response
    my $moby_exceptions = [];
    my $output_object_type  = "CommentedDNASequence";
    my $output_article_name = "upstream_sequences";
    my $namespace = "ENSEMBL";
    
    # Variables that will be passed to PromoterExtraction_call
    my @genes;
    my %parameters;
    
    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);
    
    # Get the parameters
    
    my ($organism)          = getNodeContentWithArticle($queryInput_DOM, "Parameter", "organism");
    my ($upstream_length)   = getNodeContentWithArticle($queryInput_DOM, "Parameter", "upstream length");
    my ($downstream_length) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "downstream length");
    my ($intergenic_only)   = getNodeContentWithArticle($queryInput_DOM, "Parameter", "intergenic only");
    
    (not defined $organism)          and $organism          = "Homo sapiens";
    (not defined $upstream_length)   and $upstream_length   = 2000;
    (not defined $downstream_length) and $downstream_length = 0;
    (not defined $intergenic_only)   and $intergenic_only   = "False";
    
    # Add the parsed parameters in a hash table
    
    $parameters{organism}          = $organism;
    $parameters{dbrelease}         = $dbrelease;
    $parameters{upstream_length}   = $upstream_length;
    $parameters{downstream_length} = $downstream_length;
    $parameters{intergenic_only}   = $intergenic_only;
    $parameters{orthologous_mode}  = $orthologous_mode;

    # Tratamos a cada uno de los articulos
    foreach my $article (@articles) {
	
	# El articulo es una tupla que contiene el nombre de este 
	# y su texto xml. 
	
	my ($articleName, $DOM) = @{$article}; # get the named article
	
	if (isCollectionArticle ($DOM)) {
	    
	    # not allowed
	    
	    my $note = "Received a collection input article instead of a simple";
	    print STDERR "$note\n";
	    my $code = "201";
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      refElement => "genes",
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
	
	if (($articleName =~ /gene/) || isSimpleArticle ($DOM)) { 
	    if ($_debug) {
		print STDERR "\"genes\" tag is a simple article...\n";
		print STDERR "node ref, " . ref ($DOM) . "\n";
		print STDERR "DOM: " . $DOM->toString () . "\n";
	    }
	    
	    # Validate the type first
	    
	    my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "List_Text");
	    if (!$rightType) {
		my $note = "Expecting a List_Text object, and receiving a $inputDataType object";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => "genes",
									  code       => $code,
									  type       => 'error',
									  queryID    => $queryID,
									  message    => "$note",
									  );
		push (@$moby_exceptions, $moby_exception);
		
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
	    }
	    
	    my $genes_lst = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "List_Text");
	    my @genes_tmp = split ("\n", $genes_lst);
	    
	    if ($_debug) {
		print STDERR "got a list before filtering of " . @genes_tmp . " genes\n";
	    }
	    
	    # Filter the list from empty elements
	    
	    @genes = map { /(\S+)/ } @genes_tmp;
	    
	    if ($_debug) {
		print STDERR "got a list of " . @genes . " genes\n";
	    }
	    
	    if ($_debug) {
		print STDERR "genes_lst,\n@genes.\n";
	    }
	} # End parsing genes tag
    } # End parsing articles
    
    if (@genes < 1) {
	my $note = "Error parsing input element, 'genes', the article doesn't contain any gene identifier!";
	print STDERR "$note\n";
	my $code = "201";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "genes",
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    # Una vez recogido todos los parametros necesarios, llamamos a 
    # la funcion que nos devuelve el report. 	
    
    my ($fasta_sequences, $moby_exceptions_tmp) = PromoterExtraction_call (genes => \@genes, debug => $_debug, queryID => $queryID, parameters => \%parameters);
    push (@$moby_exceptions, @$moby_exceptions_tmp);
    
    # Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX 
    # nos queda encapsularla en un Objeto bioMoby. Esta operacio 
    # la podriamos realizar en una funcion a parte si fuese compleja.  
    
    if (not defined $fasta_sequences) {
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
    }
    else {
	my $sequence_objects  = INB::GRIB::Utils::CommonUtilsSubs->createSequenceObjectsFromFASTA ($fasta_sequences, $output_object_type, $namespace);
	
	# Bien!!! ya tenemos el objeto de salida del servicio , solo nos queda
	# volver a encapsularlo en un objeto biomoby de respuesta. Pero 
	# en este caso disponemos de una funcion que lo realiza. Si tuvieramos 
	# una respuesta compleja (de verdad, esta era simple ;) llamariamos 
	# a collection response. 
	# IMPORTANTE: el identificador de la respuesta ($queryID) debe ser 
	# el mismo que el de la query. 
	
	$MOBY_RESPONSE .= collectionResponse($sequence_objects, $output_article_name, $queryID);
    }
    
    return ($MOBY_RESPONSE, $moby_exceptions);
}

##################################################################################################################
#
# Public methods
#
##

=head2 getUpstreamSeqFromEnsembl

 Title   : getUpstreamSeqFromEnsembl
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
         : queries. Es decir, un mensaje xml de la forma:
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

sub getUpstreamSeqFromEnsembl {

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
        my $MOBY_RESPONSE = "";             # set empty response
	my $moby_exceptions = [];
	
	my $moby_logger = get_logger ("MobyServices");
	my $serviceName = "getUpstreamSeqFromEnsembl";
	
	# Para cada query ejecutaremos el _execute_query.
        foreach my $queryInput(@queries) {
	    
	    if ($_debug) {
		print STDERR "upstream sequence retrieval service queryInput,\n";
		print STDERR $queryInput->toString() . "\n";
	    }
	    
	    # En este punto es importante recordar que el objeto $query 
	    # es un XML::DOM::Node, y que si queremos trabajar con 
	    # el mensaje de texto debemos llamar a: $query->toString() 
	    my ($query_response, $moby_exceptions_tmp) = _do_query_PromoterExtraction ($queryInput, "False");
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

sub getOrthologousUpstreamSeqFromEnsembl {

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
        my $MOBY_RESPONSE = "";             # set empty response
	my $moby_exceptions = [];
	
	my $moby_logger = get_logger ("MobyServices");
	my $serviceName = "getOrthologousUpstreamSeqFromEnsembl";
	
	# Para cada query ejecutaremos el _execute_query.
        foreach my $queryInput(@queries) {
	    
	    if ($_debug) {
		print STDERR "upstream sequence retrieval service queryInput,\n";
		print STDERR $queryInput->toString() . "\n";
	    }
	    
	    # En este punto es importante recordar que el objeto $query 
	    # es un XML::DOM::Node, y que si queremos trabajar con 
	    # el mensaje de texto debemos llamar a: $query->toString() 
	    my $query_response = _do_query_PromoterExtraction ($queryInput, "True");
	    
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

