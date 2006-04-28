# $Id: VectorScreeningServices.pm,v 1.1 2006-04-28 13:32:57 gmaster Exp $
#
# This file is an instance of a template written
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#
# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::VectorScreeningServices  - Package for parser the Moby message to call the Vector Screening tools.

=head1 SYNOPSIS

Con este package podremos parsear las llamadas xml de BioMoby para
poder llamar al servicio runVectorScreeningGFF. Una vez que tengamos la salida de llamando
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
  my $result = VectorScreening("call", $in);


  # Esta llamada/s nos devuelve una variable que contiene el texto con la
  # salida del programa Geneid encapsulado en objetos Moby.

=head1 DESCRIPTION

Este package sirve para parsear las llamadas BioMoby de entrada y salida.
De esta forma hace de puente entre el cliente que llama el servicio y la
aplicacion VectorScreening.

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

package INB::GRIB::Services::VectorScreeningServices;

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
  &runCrossMatchToScreenVector
  &runCrossMatchToScreenVectorCollection
);

our $VERSION = '1.0';

my $_debug = 0;

# Preloaded methods go here.

###############################################################################

sub _do_query_CrossMatchToScreenVector {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby
    my $queryInput_DOM = shift @_;
    # $_format is the type of output
    my $output_format = shift @_;
    # $input_type tells if the service specifies a collection of input objects or just a simple input object
    # We need this because the two services, runDust and runDustCollection call both this method, so they share the same code
    # But the output type will have to be different
    my $input_type     = shift @_;
    
    if ($_debug) {
	print STDERR "CrossMatch To Screen Vector execution...\n";
    }
    
    if (($input_type ne "simple") && ($input_type ne "collection")) {
	# Don't know the type, don't know what to return !
	print STDERR "input_type unknown (should be setup as being simple or collection)\n";
	exit 0;
    }

    my $MOBY_RESPONSE   = "";     # set empty response
    my $moby_exceptions = [];
    my $output_article_name = "screened_sequence";
    if ($input_type eq "collection") {
	$output_article_name    = "screened_sequences";
    }
    my $namespace       = "";
    
    # Aqui escribimos las variables que necesitamos para la funcion.
    my $minmatch;
    my $minscore;
    
    # Variables that will be passed to MatScan_call
    my %sequences;
    my %parameters;
    
    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);

    # Get the parameters
    
    ($minmatch) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "minmatch");
    if (not defined $minmatch) {
	$minmatch = 12;
    }
    elsif (! ($minmatch >= 0)) {
	my $note = "minmatch parameter, '$minmatch', not accepted an integer greater or equal than 0\n";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "minmatch",
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	# Return an empty moby data object, as well as an exception telling what nothing got returned
	
	if ($input_type =~ /simple/i) {
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_RESPONSE ($queryID, $output_article_name);
	}
	else {
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
	}
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    ($minscore) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "minscore");
    if (not defined $minscore) {
	$minscore = 20;
    }
    elsif (! ($minscore >= 0)) {
	my $note = "minscore parameter, '$minscore', not accepted an integer greater or equal than 0\n";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "minscore",
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	# Return an empty moby data object, as well as an exception telling what nothing got returned
	
	if ($input_type =~ /simple/i) {
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_RESPONSE ($queryID, $output_article_name);
	}
	else {
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
	}
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    # Add the parsed parameters in a hash table
    
    $parameters{minmatch} = $minmatch;
    $parameters{minscore} = $minscore;
    
    if ($_debug) {
	print STDERR "minmatch, $minmatch\n";
	print STDERR "minscore, $minscore\n";
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

	if (($articleName =~ /sequence/i) || (isSimpleArticle ($DOM) || (isCollectionArticle ($DOM)))) {

	    if (isSimpleArticle ($DOM)) {
		
		if ($_debug) {
		    print STDERR "$articleName tag is a simple article...\n";
		}
		
		if ($input_type eq "collection") {
		    
		    # not allowed
		    
		    my $note = "Received a collection input article instead of a simple";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => 'sequence',
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
		
		# Validate the type
		my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "DNASequence");
		if (!$rightType) {
		    my $note = "Expecting DNASequence objects, and receiving $inputDataType objects";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => "sequence",
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
		
		%sequences = INB::GRIB::Utils::CommonUtilsSubs->parseMobySequenceObjectFromDOM ($DOM, \%sequences);
	    }
	    elsif (isCollectionArticle ($DOM)) {

		if ($_debug) {
		    print STDERR "$articleName is a collection article...\n";
		    print STDERR "Collection DOM: " . $DOM->toString() . "\n";
		}
		
		if ($input_type eq "simple") {
		    # is it allowed ? - check this out
		    # at the moment - disallow this !!
		    # Anyway the client is supposed to deal with this before the execution of the service
		    
		    my $note = "Received a simple input article instead of a collection";
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
		    
		    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
		    return ($MOBY_RESPONSE, $moby_exceptions);
		}
		
		# Validate the type of the simples in the collection - should all be DNASequence objects
		my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "DNASequence");
		if (!$rightType) {
		    my $note = "Expecting a DNASequence object, and receiving a $inputDataType object";
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

		    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
		    return ($MOBY_RESPONSE, $moby_exceptions);
		}
		
		my @sequence_articles_DOM = getCollectedSimples ($DOM);
		foreach my $sequence_article_DOM (@sequence_articles_DOM) {
		   %sequences = INB::GRIB::Utils::CommonUtilsSubs->parseMobySequenceObjectFromDOM ($sequence_article_DOM, \%sequences);
		}
	    }
	    else {
		print STDERR "not supposed to ge here!\n";
		print STDERR "not a simple or a collection, what is it there then!!\n";
		print STDERR "DOM: " . $DOM->toString() . "\n";
	    }
	} # End parsing sequences article tag
	
    } # Next article
    
    # Check that we have parsed properly the sequences
    
    if ((keys (%sequences)) == 0) {
	my $note = "can't parsed any sequences...\n";
	print STDERR "$note\n";
	my $code = "201";

	my $refElement;
	if ($input_type eq "simple") {
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
	    $refElement = "sequence";
	}
	else {
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_SIMPLE_RESPONSE ($queryID, $output_article_name);
	    $refElement = "sequences"
	}
	
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => $refElement,
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    # Una vez recogido todos los parametros necesarios, llamamos a
    # la funcion que nos devuelve el report.
    
    if ($input_type eq "simple") {
	
	if ($_debug) {
	    print STDERR "making a simple response\n";
	}
	
	my ($screened_fasta_seqs, $moby_exceptions_tmp) = CrossMatchToScreenVector_call (sequences  => \%sequences, queryID => $queryID, parameters => \%parameters, debug => $_debug);
	push (@$moby_exceptions, @$moby_exceptions_tmp);
	
	# Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX
	# nos queda encapsularla en un Objeto bioMoby. Esta operacio
	# la podriamos realizar en una funcion a parte si fuese compleja.
	
	if (defined $screened_fasta_seqs) {
	    my ($sequenceIdentifier) = keys (%sequences);
	    
	    # Convert FASTA into a DNASequence object
	    my $moby_seqobjs = INB::GRIB::Utils::CommonUtilsSubs->createSequenceObjectsFromFASTA ($screened_fasta_seqs, $output_format, $namespace);
	    if (@$moby_seqobjs != 1) {
		my $note = "Can't parse any sequences, Please check the syntax of your fasta sequence.\n";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => "sequence",
									  code       => $code,
									  type       => 'error',
									  queryID    => $queryID,
									  message    => "$note",
									  );
		push (@$moby_exceptions, $moby_exception);
		
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_SIMPLE_RESPONSE ($queryID, $output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
	    }
	    
	    my $moby_seqobj = $moby_seqobjs->[0];
            $MOBY_RESPONSE  = simpleResponse($moby_seqobj, $output_article_name, $queryID);   
	}
	else {
	    my $note = "VectorScreening failed, Please check the syntax of your fasta input sequence.\n";
	    print STDERR "$note\n";
	    my $code = "201";
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      refElement => "sequence",
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    push (@$moby_exceptions, $moby_exception);
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
    elsif ($input_type eq "collection") {
	
	if ($_debug) {
	    print STDERR "making a collection response\n";
	}
	
	# The input was a collection of Sequences, so we have to return a collection of GFF objects
	
	my ($screened_seqs_fasta, $moby_exceptions_tmp) = CrossMatchToScreenVector_call (sequences  => \%sequences, format => $output_format, parameters => \%parameters, debug => $_debug);
	push (@$moby_exceptions, @$moby_exceptions_tmp);
	
	if (defined $screened_seqs_fasta) {
	    # Convert FASTA sequences into DNASequences
	    my $moby_seqobjs = INB::GRIB::Utils::CommonUtilsSubs->createSequenceObjectsFromFASTA ($screened_seqs_fasta, $output_format, $namespace);
	    $MOBY_RESPONSE   = collectionResponse($moby_seqobjs, $output_article_name, $queryID);
	}
	else {
	    my $note = "VectorScreening failed, Please check the syntax of your fasta input sequences.\n";
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
}


=head2 runCrossMatchToScreenVector

 Title   : runCrossMatchToScreenVector
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

sub runCrossMatchToScreenVector {
	
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
    # (The MatScan output format for this service is by default GFF - right now it is hardcoded)
    #

    my $_moby_output_format   = "DNASequence";
    my $moby_logger = get_logger ("MobyServices");
    my $serviceName = "runCrossMatchToScreenVector";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput(@queries){

	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}

	# En este punto es importante recordar que el objeto $query
	# es un XML::DOM::Node, y que si queremos trabajar con
	# el mensaje de texto debemos llamar a: $query->toString()
	my ($query_response, $moby_exceptions_tmp) = _do_query_CrossMatchToScreenVector ($queryInput, $_moby_output_format, "simple");
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

=head2 runCrossMatchToScreenVectorCollection

 Title   : runCrossMatchToScreenVectorCollection
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

sub runCrossMatchToScreenVectorCollection {
    
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
    
    my $_moby_output_format = "DNASequence";
    my $moby_logger = get_logger ("MobyServices");
    my $serviceName = "runCrossMatchToScreenVector";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput(@queries){

	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}

	# En este punto es importante recordar que el objeto $query
	# es un XML::DOM::Node, y que si queremos trabajar con
	# el mensaje de texto debemos llamar a: $query->toString()
	my ($query_response, $moby_exceptions_tmp) = _do_query_CrossMatchToScreenVector ($queryInput, $_moby_output_format, "collection");
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
