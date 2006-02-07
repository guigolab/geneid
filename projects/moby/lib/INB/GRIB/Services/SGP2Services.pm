# $Id: SGP2Services.pm,v 1.13 2006-02-07 12:13:00 gmaster Exp $
#
# This file is an instance of a template written 
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#
# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::SGP2Services  - Package for parser the Moby message to call the SGP2 gene finding program .

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

Copyright (c) 2005, Arnaud Kerhornou and INB - Nodo Genomico GRIB/IMIM.
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

package INB::GRIB::Services::SGP2Services;

use strict;
use warnings;
use Carp;

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
  &runSGP2GFF 
);

our $VERSION = '1.0';

my $_debug = 0;

# Preloaded methods go here.

###############################################################################

=head2 _do_query_SGP2

 Title   : _do_query_SGP2
         : 
         : private function (NOT EXPORTED)
         : 
 Usage   : my $query_response = _do_query_SGP2($query);
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

sub _do_query_SGP2 {
	# $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby 
	my $queryInput_DOM = shift @_;
	# $_format is the type of output that returns SGP2 (e.g. GFF)
	my $_format        = shift @_;
	
        my $MOBY_RESPONSE   = "";     # set empty response
	my $moby_exceptions = [];
	my $output_article_name  = "geneid_predictions";
	
	# Variables that will be passed to SGP2_call
	my %sequences;
	my $tblastx_output;
	my %parameters;

	my $queryID  = getInputID ($queryInput_DOM);
	my @articles = getArticles($queryInput_DOM);

	# Get the parameters

	# No parameters available yet...

	# Tratamos a cada uno de los articulos
	foreach my $article (@articles) {       
	    
	    # El articulo es una tupla que contiene el nombre de este 
	    # y su texto xml. 
	    
	    my ($articleName, $DOM) = @{$article}; # get the named article
	    
	    # Si le hemos puesto nombre a los articulos del servicio,  
	    # podemos recoger a traves de estos nombres el valor. 
	    
	    # Don't allow the use of collections...
	    
	    if (isCollectionArticle ($DOM)) {
		my $note = "Received a collection input article instead of a simple";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
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
	    
	    if ($articleName eq "sequences") { 
		
		if ($_debug) {
		    print STDERR "parsing the article \"sequences\"...\n";
		}
		# Validate the type first
		
		my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "NucleotideSequence");
		if (!$rightType) {
		    my $note = "Expecting a NucleotideSequence object, and receiving a $inputDataType object";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => "sequences",
									      code       => $code,
									      type       => 'error',
									      queryID    => $queryID,
									      message    => "$note",
									      );
		    push (@$moby_exceptions, $moby_exception);
		    
		    # Simple Response doesn't fit !! (the simple article is not empty as it should be!), so we need to create the string from scratch !
		    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_RESPONSE ($queryID, $output_article_name);
		    return ($MOBY_RESPONSE, $moby_exceptions);
		}
		else {
		    # print STDERR "type is fine\n";
		}
		
		# Get the sequence
		%sequences = INB::GRIB::Utils::CommonUtilsSubs->parseMobySequenceObjectFromDOM ($DOM, \%sequences);
		
		if (keys (%sequences) == 0) {
		    my $note = "Couldn't parse any sequence";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => "sequences",
									      code       => $code,
									      type       => 'error',
									      queryID    => $queryID,
									      message    => "$note",
									      );
		    push (@$moby_exceptions, $moby_exception);
		    
		    # Simple Response doesn't fit !! (the simple article is not empty as it should be!), so we need to create the string from scratch !
		    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_RESPONSE ($queryID, $output_article_name);
		    return ($MOBY_RESPONSE, $moby_exceptions);
		}
		
	    }
	    elsif ($articleName eq "tblastx") {
		
		if ($_debug) {
		    print STDERR "parsing the article \"tblastx\"...\n";
		}
		
		# Validate the type first

		my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "Blast-Text");
		if (!$rightType) {
		    my $note = "Expecting a Blast-Text object, and receiving a $inputDataType object";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => "tblastx",
									      code       => $code,
									      type       => 'error',
									      queryID    => $queryID,
									      message    => "$note",
									      );
		    push (@$moby_exceptions, $moby_exception);
		    
		    # Simple Response doesn't fit !! (the simple article is not empty as it should be!), so we need to create the string from scratch !
		    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_RESPONSE ($queryID, $output_article_name);
		    return ($MOBY_RESPONSE, $moby_exceptions);
		}
		else {
		    # print STDERR "type is fine\n";
		}
		
	    	$tblastx_output = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "Blast-Text");
		
		if (length ($tblastx_output) < 1) {
		    my $note = "Couldn't parse any tblastx data";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => "tblastx",
									      code       => $code,
									      type       => 'error',
									      queryID    => $queryID,
									      message    => "$note",
									      );
		    push (@$moby_exceptions, $moby_exception);
		    
		    # Simple Response doesn't fit !! (the simple article is not empty as it should be!), so we need to create the string from scratch !
		    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_RESPONSE ($queryID, $output_article_name);
		    return ($MOBY_RESPONSE, $moby_exceptions);
		}
	    }
	} # Next article
	
	# Una vez recogido todos los parametros necesarios, llamamos a 
	# la funcion que nos devuelve el report. 	

	# Set up the temporary file here and give it as a parameter or give a hash of sequence entries, so we can submit multiple sequences...
	
	my ($report, $moby_exceptions_tmp) = SGP2_call (sequences  => \%sequences, tblastx_output => $tblastx_output, format => $_format, queryID => $queryID, parameters => \%parameters);
	push (@$moby_exceptions, @$moby_exceptions_tmp);
	
	# Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX 
	# nos queda encapsularla en un Objeto bioMoby. Esta operacio 
	# la podriamos realizar en una funcion a parte si fuese compleja.  
	
	
	my $input = undef;
	if (defined $report) {
	    
	    # Quick hack to add the sequence identifier
	    # Anyway even if the parsing code handles input collection, the output report code doesn't and the runGeneIDGFF service registration specs tell that the input and the output are simple articles !!
	    my ($sequenceIdentifier) = keys (%sequences);
	    $input = <<PRT;
<moby:$_format namespace='' id='$sequenceIdentifier'>
<![CDATA[
$report
]]>
</moby:$_format>
PRT

            $MOBY_RESPONSE = simpleResponse($input, $output_article_name, $queryID);
        }
	else {
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_RESPONSE ($queryID, $output_article_name);
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




=head2 runSGP2GFF

 Title   : runSGP2GFF 
 Usage   : Esta función está pensada para llamarla desde un cliente SOAP. 
         : No obstante, se recomienda probarla en la misma máquina, antes 
         : de instalar el servicio. Para ello, podemos llamarla de la 
         : siguiente forma:
         : 
         : my $result = SGP2("call", $in);
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

sub runSGP2GFF {

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
    # The output format for this service is GFF
    #
    my $_format = "GFF";
    my $moby_logger = get_logger ("MobyServices");
    my $serviceName = "runSGP2GFF";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){
	
	# En este punto es importante recordar que el objeto $query 
	# es un XML::DOM::Node, y que si queremos trabajar con 
	# el mensaje de texto debemos llamar a: $query->toString() 
	
	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}
	
	my ($query_response, $moby_exceptions_tmp) = _do_query_SGP2 ($queryInput, $_format);
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
	$moby_logger->info ("$serviceName terminated successfully");
	$moby_logger->info ("Exception code, 700");

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

