# $Id: MemeServices.pm,v 1.12 2006-01-31 10:37:05 gmaster Exp $
#
# This file is an instance of a template written
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#
# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::MemeServices  - Package for parser the Moby message to call the MEME motif search progam.

=head1 SYNOPSIS



=head1 DESCRIPTION



=head1 AUTHOR

Arnaud Kerhornou, akerhornou@imim.es

=head1 COPYRIGHT

Copyright (c) 2005, Arnaud Kerhornou and INB - Nodo INB 1 - GRIB/IMIM.
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

package INB::GRIB::Services::MemeServices;

use strict;
use warnings;
use Carp;

use INB::GRIB::Services::Factory;
use INB::GRIB::Utils::CommonUtilsSubs;
use MOBY::CommonSubs qw(:all);

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
  &runMemeHTML
  &runMemeText
  &parseMotifMatricesfromMEME
);

our $VERSION = '1.0';

my $_debug = 0;

# Preloaded methods go here.

###############################################################################

=head2 _do_query_Meme

 Title   : _do_query_Meme
	 :
	 : private function (NOT EXPORTED)
	 :
 Usage   : my $query_response = _do_query_Meme($query, $output_format);
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

sub _do_query_Meme {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby
    my $queryInput_DOM = shift @_;
    # $_format is the type of output that returns Meme (e.g. HTML or Text)
    my $_output_format        = shift @_;

    my $MOBY_RESPONSE   = "";     # set empty response
    my $moby_exceptions = [];
    
    # Aqui escribimos las variables que necesitamos para la funcion.
    my $motif_distribution;
    my $maximum_number_motifs;
    my $minimum_number_sites;
    my $maximum_number_sites;
    my $minimum_motif_width;
    my $maximum_motif_width;
    my $e_value_cutoff;
    my $background_order;

    # Variables that will be passed to MatScan_call
    my %sequences;
    my %parameters;

    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);

    # Get the parameters

    ($motif_distribution) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "motif distribution");
    if (not defined $motif_distribution) {
	$motif_distribution = "zero or one";
    }

    ($maximum_number_motifs) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "maximum number of motifs");
    if (not defined $maximum_number_motifs) {
	$maximum_number_motifs = 8;
    }

    ($minimum_number_sites) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "minimum sites for each motif");
    if (not defined $minimum_number_sites) {
	# no default
    }
    # Check for constraint
    if (defined $minimum_number_sites && ($minimum_number_sites < 2)) {
	print STDERR "Error, the minimum number of sites should be at least 2\n";
	exit 1;
    }

    ($maximum_number_sites) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "maximum sites for each motif");
    if (not defined $maximum_number_sites) {
	# no default
    }
    # Check for constraint
    if (defined $maximum_number_sites && ($maximum_number_sites > 300)) {
	print STDERR "Error, the maximum number of sites is too big, should be at most 300\n";
	exit 1;
    }

    ($minimum_motif_width) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "minimum optimum width");
    if (not defined $minimum_motif_width) {
	$minimum_motif_width = 6;
    }
    # Check for constraint
    if (defined $minimum_motif_width && ($minimum_motif_width < 2)) {
	print STDERR "Error, the minimum width should be at least 2\n";
	exit 1;
    }

    ($maximum_motif_width) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "maximum optimum width");
    if (not defined $maximum_motif_width) {
	$maximum_motif_width = 15;
    }
    # Check for constraint
    if (defined $maximum_motif_width && ($maximum_motif_width > 300)) {
	print STDERR "Error, the maximum width is too big, should be at most 300\n";
	exit 1;
    }

    ($e_value_cutoff) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "motif E-value cutoff");
    if (not defined $e_value_cutoff) {
	$e_value_cutoff = "1";
    }

    ($background_order) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "background markov model training (value is the model order)");
    if (not defined $background_order) {
	$background_order = "1";
    }

    # Add the parsed parameters in a hash table

    $parameters{motif_distribution}    = $motif_distribution;
    $parameters{maximum_number_motifs} = $maximum_number_motifs;
    $parameters{minimum_number_sites}  = $minimum_number_sites;
    $parameters{maximum_number_sites}  = $maximum_number_sites;
    $parameters{minimum_motif_width}   = $minimum_motif_width;
    $parameters{maximum_motif_width}   = $maximum_motif_width;
    $parameters{e_value_cutoff}        = $e_value_cutoff;
    $parameters{background_order}      = $background_order;

    if ($_debug) {
	print STDERR "motif distribution, $motif_distribution\n";
	print STDERR "maximum number motifs, $maximum_number_motifs\n";
	print STDERR "minimum number sites, $minimum_number_sites\n";
	print STDERR "maximum number sites, $maximum_number_sites\n";
	print STDERR "minimum motif width, $minimum_motif_width\n";
	print STDERR "maximum motif width, $maximum_motif_width\n";
	print STDERR "E-value cutoff, $e_value_cutoff\n";
	print STDERR "background model order, $background_order\n";
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

	if (($articleName eq "upstream_sequences") || (isSimpleArticle ($DOM) || (isCollectionArticle ($DOM)))) {

	    if (isSimpleArticle ($DOM)) {

		print STDERR "sequences tag is a simple article, should be a collection...\n";
		exit 1;
	    }
	    elsif (isCollectionArticle ($DOM)) {

		if ($_debug) {
		    print STDERR "sequences is a collection article...\n";
		    print STDERR "Collection DOM: " . $DOM->toString() . "\n";
		}

		my @sequence_articles_DOM = getCollectedSimples ($DOM);

		foreach my $sequence_article_DOM (@sequence_articles_DOM) {
		    %sequences = INB::GRIB::Utils::CommonUtilsSubs->parseMobySequenceObjectFromDOM ($sequence_article_DOM, \%sequences);
		}
	    }
	    else {
		print STDERR "It is not a simple or collection article...\n";
		print STDERR "DOM: " . $DOM->toString() . "\n";
	    }

	} # End parsing sequences article tag

    } # Next article

    # Check that we have parsed properly the sequences

    if ((keys (%sequences)) == 0) {
	print STDERR "Error, can't parsed any sequences...\n";
    }

    # Una vez recogido todos los parametros necesarios, llamamos a
    # la funcion que nos devuelve el report.

    my $output_article_name = "meme_predictions";

    if ($_debug) {
	print STDERR "making a simple response\n";
    }

    my ($report, $moby_exceptions_tmp) = MEME_call (sequences  => \%sequences, format => $_output_format, debug => $_debug, queryID => $queryID, parameters => \%parameters);
    push (@$moby_exceptions, @$moby_exceptions_tmp);
    
    # Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX
    # nos queda encapsularla en un Objeto bioMoby. Esta operacio
    # la podriamos realizar en una funcion a parte si fuese compleja.

    my $input = <<PRT;
<moby:$_output_format namespace='' id=''>
<![CDATA[
$report
]]>
</moby:$_output_format>
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

sub _do_query_MemeMotifMatrices {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby
    my $queryInput_DOM = shift @_;
    my $_output_format        = shift @_;

    my $MOBY_RESPONSE   = "";     # set empty response
    my $moby_exceptions = [];
    
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

    # Add the parsed parameters in a hash table

    $parameters{matrix_mode}    = $matrix_mode;

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

	    $meme_predictions = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "text-formatted");

	    if ($_debug) {
		print STDERR "parsed meme_predictions,\n$meme_predictions.\n";
	    }

	}
	if (isCollectionArticle ($DOM)) {
	    print STDERR "Error, collection input are not expected!\n";
	    exit 0;
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
    my $output_article_name = "meme_matrices";
    my $namespace = "";
    
    if (not defined $matrices_aref) {
	# Return an emtpy message !
	$MOBY_RESPONSE .= simpleResponse (undef, $output_article_name, $queryID);
	return $MOBY_RESPONSE;
    }
    
    # Concatenate all matrices of the array into one string
    my $matrices = join ('', @$matrices_aref);
    my $meme_matrix_object = <<PRT;
<moby:$_output_format namespace='$namespace' id=''>
<![CDATA[
$matrices
]]>
</moby:$_output_format>
PRT

    # Bien!!! ya tenemos el objeto de salida del servicio , solo nos queda
    # volver a encapsularlo en un objeto biomoby de respuesta. Pero
    # en este caso disponemos de una funcion que lo realiza. Si tuvieramos
    # una respuesta compleja (de verdad, esta era simple ;) llamariamos
    # a collection response.
    # IMPORTANTE: el identificador de la respuesta ($queryID) debe ser
    # el mismo que el de la query.
    
    $MOBY_RESPONSE .= simpleResponse($meme_matrix_object, $output_article_name, $queryID);
    return ($MOBY_RESPONSE, $moby_exceptions);
    
}

=head2 runMemeHTML

 Title   : runMemeHTML
 Usage   : Esta función está pensada para llamarla desde un cliente SOAP.
	 : No obstante, se recomienda probarla en la misma máquina, antes
	 : de instalar el servicio. Para ello, podemos llamarla de la
	 : siguiente forma:
	 :
	 : my $result = Meme("call", $in);
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
 Returns : Return a text-html object with MEME results

=cut

sub runMemeHTML {

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
    
    #
    # The moby output format for this service is text-html
    #
    my $_moby_output_format   = "text-html";

    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){

	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}

	# En este punto es importante recordar que el objeto $query
	# es un XML::DOM::Node, y que si queremos trabajar con
	# el mensaje de texto debemos llamar a: $query->toString()
	my ($query_response, $moby_exceptions_tmp) = _do_query_Meme ($queryInput, $_moby_output_format);
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

=head2 parseMotifMatricesfromMEME

 Title   : parseMotifMatricesfromMEME
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

sub parseMotifMatricesfromMEME {

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
    
    my $_format = "text-formatted";

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


=head2

 Title   : runMemeText
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
 Returns : Return a text-formatted object with MEME results

=cut

sub runMemeText {

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
    
    my $_format = "text-formatted";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){
	
	# En este punto es importante recordar que el objeto $query
	# es un XML::DOM::Node, y que si queremos trabajar con
	# el mensaje de texto debemos llamar a: $query->toString()

	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}
	
	my ($query_response, $moby_exceptions_tmp) = _do_query_Meme ($queryInput, $_format);
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

1;

__END__
