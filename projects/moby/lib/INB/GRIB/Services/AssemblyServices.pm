# $Id: AssemblyServices.pm,v 1.4 2007-03-13 18:20:07 gmaster Exp $
#
# This file is an instance of a template written
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#
# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::AssemblyServices  - Package for parser the Moby message to call the Assembly program .

=head1 SYNOPSIS

Con este package podremos parsear las llamadas xml de BioMoby para
poder llamar al servicio runAssemblyGFF. Una vez que tengamos la salida de llamando
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

package INB::GRIB::Services::AssemblyServices;

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
  &runPhrap
  &runPhrapWithQualityData
);

our $VERSION = '1.0';

my $_debug = 0;

# Preloaded methods go here.

###############################################################################

=head2 _do_query_Phrap

 Title   : _do_query_Phrap
	 :
	 : private function (NOT EXPORTED)
	 :
 Usage   : my $query_response = _do_query_Phrap($query);
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

sub _do_query_Phrap {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby
    my $queryInput_DOM   = shift @_;
    my $sequences_input_format  = shift @_;
    my $sequences_output_format = shift @_;
    my $phrap_output_format     = shift @_;
    # ServiceName tells us whether we receive quality data or just DNA sequences 
    my $serviceName      = shift @_;
    
    if ($_debug) {
      print STDERR "Phrap...\n";
    }
    
    my $MOBY_RESPONSE   = "";
    my $moby_exceptions = [];
    my $sequences_output_article_name = "contig_and_singlet_sequences";
    my $phrap_output_article_name     = "assembly";
    my $namespace = "";
    
    # Aqui escribimos las variables que necesitamos para la funcion.
    my $node_space;
    my $node_seg;
    
    # Variables that will be passed to Phrap_call
    my $fasta_quality_data_str;
    my $fasta_sequences_str;
    my %parameters;
    
    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);
    
    # Get the parameters
    
    ($node_seg) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "node_seg");
    if (not defined $node_seg) {
	$node_seg = 8;
    }
    elsif (! ($node_seg =~ /\d+/)) {
	my $note = "node_seg parameter, '$node_seg', not numerical";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "node_seg",
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	# Return an empty moby data object, as well as an exception telling what nothing got returned
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $phrap_output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    ($node_space) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "node_space");
    if (not defined $node_space) {
	$node_space = 4;
    }
    elsif (! ($node_space =~ /\d+/)) {
	my $note = "node_space parameter, '$node_space', not numerical";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "node_space",
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	# Return an empty moby data object, as well as an exception telling what nothing got returned
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $phrap_output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    # Add the parsed parameters in a hash table
    
    $parameters{node_seg}   = $node_seg;
    $parameters{node_space} = $node_space;
    
    if ($_debug) {
	print STDERR "node_seg, $node_seg\n";
	print STDERR "node_space, $node_space\n";
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
		
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $phrap_output_article_name);
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
		
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $phrap_output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
	    }
	    
	    if ($_debug) {
		print STDERR "Parsing the fasta sequences string...\n";
	    }
	    
	    $fasta_sequences_str  = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "String");
	    # Testing compression
	    # $fasta_sequences = decode_base64 ($fasta_sequences);
	    	    
	} # End parsing sequences article tag
	
	# Parse the quality data if any - check the serviceName, if not "runPhrapWithQualityData" - tell that it is the wrong service !!
	
	if ($articleName =~ /base_quality_data/i) {
	    
	    # Chech also that the service is not "runPhrap", because not allowed with that one !!
	    
	    if ($serviceName eq "runPhrapWithQualityData") {
		
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
		    
		    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $phrap_output_article_name);
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
		my $note = "Found an article with quality data, but the service, $serviceName, doesn't require them, so they won't be taken into account. Use runPhrapWithQualityData instead.";
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
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $phrap_output_article_name);
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
    
    if ($serviceName eq "runPhrapWithQualityData") {
	if (length ($fasta_quality_data_str) < 1) {
	    my $note = "can't parse any quality data...\n";
	    print STDERR "$note\n";
	    my $code = "201";
	    
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $phrap_output_article_name);
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
	
	my $nb_sequences = -1;
	my $nb_quality_sequences = -1;
	
	# nb occurrences of '>' to extract the number of sequences
	$nb_sequences = ($fasta_sequences_str =~ tr/>//);
	
	print STDERR "Number of sequences submitted for $serviceName execution is $nb_sequences\n";
	
	if ($nb_sequences > 5000) {
	    my $note = "The number of DNA sequences, $nb_sequences is not supported, the maximum limit of input sequences is 5000\n";
	    print STDERR "$note\n";
	    my $code = "201";
	    
	    $MOBY_RESPONSE     = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $phrap_output_article_name);
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code    => $code,
								      type    => 'error',
								      queryID => $queryID,
								      message => "$note",
								      );
	    push (@$moby_exceptions, $moby_exception);
	    
	    return ($MOBY_RESPONSE, $moby_exceptions);
	}
	
	$nb_quality_sequences = ($fasta_quality_data_str =~ tr/>//);
	
	if ($nb_sequences != $nb_quality_sequences) {
	    my $note = "The number of DNA sequences, $nb_sequences, doesn't match the number of quality data sequences, $nb_quality_sequences...\n";
	    print STDERR "$note\n";
	    my $code = "201";
	    
	    $MOBY_RESPONSE     = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $phrap_output_article_name);
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
	print STDERR "parsing of the Moby input data done, preparing now Phrap execution...\n";
    }
    
    my ($assembled_seqs_fasta, $ace_data, $moby_exceptions_tmp) = Phrap_call (sequences  => $fasta_sequences_str, quality_data => $fasta_quality_data_str, parameters => \%parameters, queryID => $queryID, debug => $_debug);
    push (@$moby_exceptions, @$moby_exceptions_tmp);
    
    undef $fasta_sequences_str;
    undef $fasta_quality_data_str;
    
    if (defined $assembled_seqs_fasta) {
	my $sequences_moby_seqobj = "<$sequences_output_format namespace='$namespace' id='Default'>\n<String id='' namespace='' articleName='content'><![CDATA[$assembled_seqs_fasta]]></String>\n</$sequences_output_format>\n";
	my $ace_moby_seqobj = "<$phrap_output_format namespace='$namespace' id='Default'>\n<String id='' namespace='' articleName='content'><![CDATA[$ace_data]]></String>\n</$phrap_output_format>\n";
	
	$MOBY_RESPONSE  = INB::GRIB::Utils::CommonUtilsSubs->MOBY_DOUBLE_SIMPLE_RESPONSE ($sequences_moby_seqobj, $sequences_output_article_name, $ace_moby_seqobj, $phrap_output_article_name, $queryID);
    }
    else {
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $phrap_output_article_name);
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

=head2 runPhrap

 Title   : runPhrap
 Usage   : Esta función está pensada para llamarla desde un cliente SOAP.
	 : No obstante, se recomienda probarla en la misma máquina, antes
	 : de instalar el servicio. Para ello, podemos llamarla de la
	 : siguiente forma:
	 :
	 : my $result = Phrap("call", $in);
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
 Returns : Return a collection of GFF objects, one object for each input sequence object

=cut

sub runPhrap {

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
    
    my $_sequences_input_format  = "FASTA_NA_multi";
    my $_sequences_output_format = "FASTA_NA_multi";
    my $_phrap_output_format     = "Ace_Text";
    my $moby_logger       = get_logger ("MobyServices");
    my $serviceName       = "runPhrap";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){
	
	# En este punto es importante recordar que el objeto $query
	# es un XML::DOM::Node, y que si queremos trabajar con
	# el mensaje de texto debemos llamar a: $query->toString()
	
	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}

	my ($query_response, $moby_exceptions_tmp) = _do_query_Phrap ($queryInput, $_sequences_input_format, $_sequences_output_format, $_phrap_output_format, $serviceName);
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


sub runPhrapWithQualityData {

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
    
    my $_sequences_input_format  = "FASTA_NA_multi";
    my $_sequences_output_format = "FASTA_NA_multi";
    my $_phrap_output_format     = "Ace_Text";
    my $moby_logger       = get_logger ("MobyServices");
    my $serviceName       = "runPhrapWithQualityData";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){
	
	# En este punto es importante recordar que el objeto $query
	# es un XML::DOM::Node, y que si queremos trabajar con
	# el mensaje de texto debemos llamar a: $query->toString()
	
	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}
	
	my ($query_response, $moby_exceptions_tmp) = _do_query_Phrap ($queryInput, $_sequences_input_format, $_sequences_output_format, $_phrap_output_format, $serviceName);
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
