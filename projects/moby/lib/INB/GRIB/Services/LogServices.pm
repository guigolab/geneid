# $Id: LogServices.pm,v 1.5 2007-07-18 21:36:16 arnau Exp $
#
# This file is an instance of a template written
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#
# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::LogServices  - Package for parser the Moby message to call the Log Services reporting .

=head1 SYNOPSIS

Con este package podremos parsear las llamadas xml de BioMoby para
poder llamar al servicio getStatisticalLog.

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
  my $result = Masking("call", $in);


  # Esta llamada/s nos devuelve una variable que contiene el texto con la
  # salida del programa Geneid encapsulado en objetos Moby.

=head1 DESCRIPTION

Este package sirve para parsear las llamadas BioMoby de entrada y salida.
De esta forma hace de puente entre el cliente que llama el servicio y la
aplicacion Masking.

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

Copyright (c) 2007, Arnaud Kerhornou and INB - Nodo INB 1 - GRIB/CRG.
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

package INB::GRIB::Services::LogServices;

use strict;
use warnings;
use Carp;

use INB::GRIB::Services::FactoryLog;
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
  &getStatisticalLog
);

our $VERSION = '1.0';

my $_debug = 1;

# Preloaded methods go here.

###############################################################################

sub _do_query_getStatisticalLog {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby
    my $queryInput_DOM = shift @_;
    
    if ($_debug) {
	print STDERR "getStatisticalLog...\n";
    }
    
    my $MOBY_RESPONSE   = "";     # set empty response
    my $moby_exceptions = [];
    my $output_article_name = "logReport";
    my $namespace       = "";
    my $node = "genome.imim.es";
    
    # Aqui escribimos las variables que necesitamos para la funcion.
    my $includeTests;
    my $startTime;
    my $endTime;
    
    # Variables that will be passed to MatScan_call
    my %parameters;

    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);

    # Get the parameters

    ($includeTests) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "includeTests");
    if (not defined $includeTests) {
	$includeTests = 0;
    }
    elsif (! ($includeTests =~ /true|false/)) {
	my $note = "includeTests parameter, '$includeTests', not accepted, should be ['true', 'false']";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "includeTests",
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
    elsif ($includeTests eq "true") {
      $includeTests = 1;
    }
    else {
      $includeTests = 0;
    }
    
    # Add the parsed parameters in a hash table
    
    $parameters{includeTests} = $includeTests;
    
    if ($_debug) {
	print STDERR "includeTests, $includeTests\n";
    }
    
    # Tratamos a cada uno de los articulos
    foreach my $article (@articles) {
	
	# El articulo es una tupla que contiene el nombre de este
	# y su texto xml.
	
	my ($articleName, $DOM) = @{$article}; # get the named article

	# Si le hemos puesto nombre a los articulos del servicio,
	# podemos recoger a traves de estos nombres el valor.
	# Sino sabemos que es el input articulo porque es un simple/collection articulo
	
	if ($articleName eq "start_time") {

	    if (isCollectionArticle ($DOM)) {
		
		    # not allowed
		    
		    my $note = "Received a collection input article instead of a simple";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => 'start_time',
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
		
            # Validate the type
            my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "DateTime");
            if (!$rightType) {
                my $note = "Expecting a DateTime object, and receiving $inputDataType objects";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
								              refElement => "start_time",
									      code       => $code,
									      type       => 'error',
									      queryID    => $queryID,
									      message    => "$note",
                                                                          );
	        push (@$moby_exceptions, $moby_exception);
		    
		# Empty response
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_RESPONSE ($queryID, $output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
            }
		
            # check that's correct ?
            $startTime = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "DateTime");
            
            if ($_debug) {
              print STDERR "parsed start_time value, $startTime\n";
            }
            
            $parameters{startTime} = $startTime;
        }
        if ($articleName eq "end_time") {

	    if (isCollectionArticle ($DOM)) {
		
		    # not allowed
		    
		    my $note = "Received a collection input article instead of a simple";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => 'end_time',
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
		
            # Validate the type
            my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "DateTime");
            if (!$rightType) {
                my $note = "Expecting a DateTime object, and receiving $inputDataType objects";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
								              refElement => "end_time",
									      code       => $code,
									      type       => 'error',
									      queryID    => $queryID,
									      message    => "$note",
                                                                          );
	        push (@$moby_exceptions, $moby_exception);
		    
		# Empty response
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_RESPONSE ($queryID, $output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
            }
		
            $endTime = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "DateTime");
            
            if ($_debug) {
              print STDERR "parsed end_time, $endTime\n";
            }
            
            $parameters{endTime} = $endTime;
        }
        
	
    } # Next article
    
    # Check that we have parsed properly the two input articles
    
    if ((! defined $startTime) && (! defined $endTime)) {
          my $note = "input without start_time and end_time articles. These two article inputs are required";
	  print STDERR "$note\n";
	  my $code = "201";
	  my $moby_exception = INB::Exceptions::MobyException->new (
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
    
    # Una vez recogido todos los parametros necesarios, llamamos a
    # la funcion que nos devuelve el report.
    
	my ($reports_href, $number_of_events, $moby_exceptions_tmp) = LogReport_call (queryID => $queryID, parameters => \%parameters, debug => $_debug);
	push (@$moby_exceptions, @$moby_exceptions_tmp);
	
	# Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX
	# nos queda encapsularla en un Objeto bioMoby. Esta operacio
	# la podriamos realizar en una funcion a parte si fuese compleja.
	
	my $logReport_xml = "<LogReport namespace='' id='' articleName=''>\n" . 
	  "<String namespace='' id='' articleName='node'>$node</String>\n" .
	  "<DateTime namespace='' id='' articleName='start'>$startTime</DateTime>\n" .
	  "<DateTime namespace='' id='' articleName='end'>$endTime</DateTime>\n" .
	  "<Integer namespace='' id='' articleName='numberEvents'>$number_of_events</Integer>\n";
	
	# foreach hash entry, make a LogEvent object
	
	foreach my $id (keys (%$reports_href)) {
	
	  my $logEvent_href = $reports_href->{$id};
	  
	  if (defined $logEvent_href) {
	  
	    my $start = $logEvent_href->{START};
	    my $end = $logEvent_href->{END};
	  
	    # format START and END
	  
	    my $formatted_start = "";
	    my $formatted_end = "";
	    # ...
	
	    if (($includeTests) || (!$includeTests && ($logEvent_href->{IS_TEST} eq "false"))) {

	      # sometimes status is not set, due to a bug in some of the services
	      # set the status to false then
	      # this bug is to fix
	      
	      # is that correct 1 ?
	      
	      my $status = $logEvent_href->{STATUS};
	      if (!defined $status) {
	        $status = 1;
	      }

	      if ($_debug) {
	        (! defined $logEvent_href->{SERVICE}) && print STDERR "service not defined for id, $id\n";
	        (! defined $logEvent_href->{IP}) && print STDERR "ip not defined for id, $id\n";
	        (! defined $logEvent_href->{STATUS}) && print STDERR "status not defined for id, $id - set to 1!\n";
	        (! defined $logEvent_href->{NUMBER_CPUs}) && print STDERR "nb CPUs not defined for id, $id\n";
	        (! defined $logEvent_href->{IS_TEST}) && print STDERR "is test not defined for id, $id\n";
	      }
	
	      my $logEvent = "<LogEvent namespace='' id='' articleName=''>\n" .
	        "<String namespace='' id='' articleName='id'>$id</String>\n" .
	        "<String namespace='' id='' articleName='serviceName'>" . $logEvent_href->{SERVICE} . "</String>\n" .
	        "<String namespace='' id='' articleName='ip'>" . $logEvent_href->{IP} . "</String>\n" .
	        "<DateTime namespace='' id='' articleName='start'>" . $formatted_start . "</DateTime>\n" .
	        "<DateTime namespace='' id='' articleName='end'>" . $formatted_end . "</DateTime>\n" .
	        "<Integer namespace='' id='' articleName='status'>" . $status . "</Integer>\n" .
	        "<Integer namespace='' id='' articleName='numberCPUs'>" . $logEvent_href->{NUMBER_CPUs} . "</Integer>\n" .
	        "<Boolean namespace='' id='' articleName='isTest'>" . $logEvent_href->{IS_TEST} . "</Boolean>\n" .
	        "</LogEvent>\n";
	  
              $logReport_xml .= $logEvent;
            }
          }
        } 
        
        $logReport_xml .= "</LogReport>\n";
	
	# Bien!!! ya tenemos el objeto de salida del servicio , solo nos queda
	# volver a encapsularlo en un objeto biomoby de respuesta. Pero
	# en este caso disponemos de una funcion que lo realiza. Si tuvieramos
	# una respuesta compleja (de verdad, esta era simple ;) llamariamos
	# a collection response.
	# IMPORTANTE: el identificador de la respuesta ($queryID) debe ser
	# el mismo que el de la query.
	
	$MOBY_RESPONSE = simpleResponse($logReport_xml, $output_article_name, $queryID);
	
	return ($MOBY_RESPONSE, $moby_exceptions);

}

##########################################################################

=head2 getStatisticalLog

 Title   : getStatisticalLog
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

sub getStatisticalLog {
	
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
    
    my $moby_logger = get_logger ("MobyServices");
    my $serviceName = "getStatisticalLog";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput(@queries){

	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}

	# En este punto es importante recordar que el objeto $query
	# es un XML::DOM::Node, y que si queremos trabajar con
	# el mensaje de texto debemos llamar a: $query->toString()
	my ($query_response, $moby_exceptions_tmp) = _do_query_getStatisticalLog ($queryInput);
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
