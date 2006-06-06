# $Id: BaseCallingServices.pm,v 1.3 2006-06-06 09:59:13 gmaster Exp $
#
# This file is an instance of a template written
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#
# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::BaseCallingServices  - Package for parser the Moby message to call the Phred Base Calling program .

=head1 SYNOPSIS

Con este package podremos parsear las llamadas xml de BioMoby para
poder llamar al servicio runBaseCallingGFF. Una vez que tengamos la salida de llamando
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
  my $result = BaseCalling("call", $in);


  # Esta llamada/s nos devuelve una variable que contiene el texto con la
  # salida del programa Geneid encapsulado en objetos Moby.

=head1 DESCRIPTION

Este package sirve para parsear las llamadas BioMoby de entrada y salida.
De esta forma hace de puente entre el cliente que llama el servicio y la
aplicacion BaseCalling.

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

package INB::GRIB::Services::BaseCallingServices;

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
  &runPhred
  &runPhredCollection
);

our $VERSION = '1.0';

my $_debug = 0;

# Preloaded methods go here.

sub _do_query_Phred {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby
    my $queryInput_DOM = shift @_;
    # formats
    my $sequences_moby_output_format = shift @_;
    my $quality_moby_output_format   = shift @_;
    my $input_type     = shift @_;

    if ($_debug) {
	print STDERR "Executing Phred service...\n";
    }
    
    if (($input_type ne "simple") && ($input_type ne "collection")) {
	# Don't know the type, don't know what to return !
	print STDERR "input_type unknown (should be setup as being simple or collection)\n";
	exit 1;
    }

    my $namespace       = "";
    my $MOBY_RESPONSE   = ""; # set empty response
    my $moby_exceptions = [];
    my $traces_input_article_name     = "trace";
    ($input_type eq "collection") and ($traces_input_article_name = "traces");
    my $sequences_output_article_name = "sequence";
    ($input_type eq "collection") and $sequences_output_article_name = "sequences";
    my $chromatogram_id = "";
    
    if ($_debug) {
	print STDERR "sequence output article name, $sequences_output_article_name.\n";
	print STDERR "trace input article name, $traces_input_article_name.\n";
    }
    
    my $quality_output_article_name   = "base_quality_data";
    
    # Aqui escribimos las variables que necesitamos para la funcion.
    
    my $trim_alt;
    my $trim_cutoff;
    
    # Variables that will be passed to Phred_call
    my %chromatograms;
    my %parameters;
    
    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);
    
    # Get the parameters
    
    ($trim_alt) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "trim_alt");
    if (not defined $trim_alt) {
      $trim_alt = 0;
    }
    elsif (! ($trim_alt =~ /on|off/i)) {
      my $note = "'trim_alt' parameter, '$trim_alt', not accepted, should be ['On', 'Off']\n";
      print STDERR "$note\n";
      my $code = 222;
      my $moby_exception = INB::Exceptions::MobyException->new (
                                                              refElement => "trim_alt",
                                                              code       => $code,
                                                              type       => 'error',
                                                              queryID    => $queryID,
                                                              message    => "$note",
                                                               );
      push (@$moby_exceptions, $moby_exception);
      
      # Return an empty moby data object, as well as an exception telling what nothing got returned
		
      $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $quality_output_article_name);
      return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    # Map it into a boolean
    $trim_alt =~ /^on$/i  and $trim_alt = 1;
    $trim_alt =~ /^off$/i and $trim_alt = 0;
    
    ($trim_cutoff) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "trim_cutoff");
    if (not defined $trim_cutoff) {
      $trim_cutoff = 0.05;
    }

    # Check it is a float 
    if (! ($trim_cutoff =~ /\d\.*\d*/)) {
      print STDERR "trim_cutoff parameter, $trim_cutoff, is not valid, should be a float!!\n";
    }
    
    # check that if trim_cutoff is set up, trim_alt is settup too (as a warning)
    
    $parameters{trim_alt}    = $trim_alt;
    $parameters{trim_cutoff} = $trim_cutoff;
    
    # Tratamos a cada uno de los articulos
    foreach my $article (@articles) {
	
	# El articulo es una tupla que contiene el nombre de este
	# y su texto xml.

	my ($articleName, $DOM) = @{$article}; # get the named article
	
	# make it more 'interoperable' by testing also pattern matching!!!
	# i think it should be up to the clients to adjust the article names
	# taverna right now doesn't set up the article name properly...
	# ... no checking...
	
	if (($articleName =~ /trace/i) || isSimpleArticle ($DOM) || isCollectionArticle ($DOM)) {
	    
	    if ($_debug) {
		print STDERR "parsing the trace article(s)...\n";
	    }
	    
	    if (isSimpleArticle ($DOM) && ($input_type eq "collection")) {
		
		# not allowed
		
		my $note = "Received a simple input article instead of a collection";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => $traces_input_article_name,
									  code       => $code,
									  type       => 'error',
									  queryID    => $queryID,
									  message    => "$note",
									  );
		push (@$moby_exceptions, $moby_exception);
		
		# Return an empty moby data object, as well as an exception telling what nothing got returned
		
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $quality_output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
	    }
	    elsif (isCollectionArticle ($DOM) && ($input_type eq "simple")) {
		
		# not allowed
		
		my $note = "Received a collection input article instead of a simple";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => $traces_input_article_name,
									  code       => $code,
									  type       => 'error',
									  queryID    => $queryID,
									  message    => "$note",
									  );
		push (@$moby_exceptions, $moby_exception);
		
		# Return an empty moby data object, as well as an exception telling what nothing got returned
		
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $quality_output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
	    }
	    
	    # Validate the type of the simples in the collection - should all be Chromatogram objects
	    my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "Chromatogram_Encoded");
	    if (!$rightType) {
		my $note = "Expecting a Chromatogram_Encoded object, and receiving a $inputDataType object";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => $traces_input_article_name,
									  code       => $code,
									  type       => 'error',
									  queryID    => $queryID,
									  message    => "$note",
									  );
		push (@$moby_exceptions, $moby_exception);
		
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $quality_output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
	    }
	    
	    if ($input_type eq "simple") {
		($chromatogram_id) = getSimpleArticleIDs ([$DOM]);
		
		if (! defined $chromatogram_id) {
		    my $note = "Chromatogram_Encoded object is missing an indentifier. An identifier is required";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => $traces_input_article_name,
									      code       => $code,
									      type       => 'error',
									      queryID    => $queryID,
									      message    => "$note",
									      );
		    push (@$moby_exceptions, $moby_exception);
		    
		    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $quality_output_article_name);
		    return ($MOBY_RESPONSE, $moby_exceptions);
		}
		
		my $chromatogram_data = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "String");
		# Get the Object Datatype
		my $chromatogram_type = INB::GRIB::Utils::CommonUtilsSubs->getArticleDataType ($DOM);
		
		if (! defined $chromatogram_type) {
		    print STDERR "Error, the chromatogram type is not defined!\n";
		}
		
		if ($_debug) {
		    print STDERR "chromatogram type before, $chromatogram_type\n";
		}
		    
		$chromatogram_type =~ s/_Encoded//;
		$chromatogram_type =~ s/moby://;
		
		if ($_debug) {
		    print STDERR "chromatogram type after, $chromatogram_type\n";
		}
		
		my $chromatogram = {
		                    type    => $chromatogram_type,
				    rawdata => $chromatogram_data
				   };
		
		$chromatograms{$chromatogram_id} = $chromatogram;
	    }
	    else {
		my @chromatogram_articles_DOM = getCollectedSimples ($DOM);
		foreach my $chromatogram_article_DOM (@chromatogram_articles_DOM) {
		    
		    if ($_debug) {
			print STDERR "going through the collection of chromatograms...\n";
		    }
		    
		    ($chromatogram_id) = getSimpleArticleIDs ([$chromatogram_article_DOM]);
		    
		    if (! defined $chromatogram_id) {
			my $note = "Chromatogram_Encoded object is missing an indentifier. An identifier is required";
			print STDERR "$note\n";
			my $code = "201";
			my $moby_exception = INB::Exceptions::MobyException->new (
										  refElement => $traces_input_article_name,
										  code       => $code,
										  type       => 'error',
										  queryID    => $queryID,
										  message    => "$note",
										  );
			push (@$moby_exceptions, $moby_exception);
			
			$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $quality_output_article_name);
			return ($MOBY_RESPONSE, $moby_exceptions);
		    }
		    
		    my $chromatogram_data = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($chromatogram_article_DOM, "String");
		    # Get the Object Datatype
		    my $chromatogram_type = INB::GRIB::Utils::CommonUtilsSubs->getArticleDataType ($chromatogram_article_DOM);
		    
		    if (! defined $chromatogram_type) {
			print STDERR "Error, the chromatogram type is not defined!\n";
		    }
		    
		    if ($_debug) {
			print STDERR "chromatogram type before, $chromatogram_type\n";
		    }
		    
		    $chromatogram_type =~ s/_Encoded//;
		    $chromatogram_type =~ s/moby://;
		    
		    if ($_debug) {
			print STDERR "chromatogram type after, $chromatogram_type\n";
		    }
		    
		    if (defined $chromatograms{$chromatogram_id}) {
			my $note = "Chromatogram_Encoded object identifiers are not unique. The following one, $chromatogram_id, is duplicated. Make sure all identifiers are unique";
			print STDERR "$note\n";
			my $code = "201";
			my $moby_exception = INB::Exceptions::MobyException->new (
										  refElement => $traces_input_article_name,
										  code       => $code,
										  type       => 'error',
										  queryID    => $queryID,
										  message    => "$note",
										  );
			push (@$moby_exceptions, $moby_exception);
			
			$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $quality_output_article_name);
			return ($MOBY_RESPONSE, $moby_exceptions);
		    }
		    
		    my $chromatogram = {
					type    => $chromatogram_type,
					rawdata => $chromatogram_data
					};
		    $chromatograms{$chromatogram_id} = $chromatogram;
		}
	    }
	} # End parsing chromatograms article
    } # Next article
    
    # Check that we have parsed properly the chromatograms
    
    if (keys (%chromatograms) == 0) {
	my $note = "can't parse any chromatogram data...\n";
	print STDERR "$note\n";
	my $code = "201";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => $traces_input_article_name,
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $quality_output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    my ($fasta_sequences, $fasta_quality, $moby_exceptions_tmp) = Phred_call (chromatograms  => \%chromatograms, parameters => \%parameters, queryID => $queryID, debug => $_debug);
    push (@$moby_exceptions, @$moby_exceptions_tmp);
    
    if (((defined $fasta_sequences) && ($fasta_sequences ne "")) && ((defined $fasta_quality) && ($fasta_quality ne ""))) {
	
	if ($input_type eq "collection") {
	    my $sequences_fasta_moby_object = "<$sequences_moby_output_format namespace='$namespace' id='Default'>\n<String id='' namespace='' articleName='content'><![CDATA[$fasta_sequences]]></String>\n</$sequences_moby_output_format>\n";
	    
	    my $quality_fasta_moby_object = "<$quality_moby_output_format namespace='$namespace' id='Default'>\n<String id='' namespace='' articleName='content'><![CDATA[$fasta_quality]]></String>\n</$quality_moby_output_format>\n";
	    
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_DOUBLE_SIMPLE_RESPONSE ($sequences_fasta_moby_object, $sequences_output_article_name, $quality_fasta_moby_object, $quality_output_article_name, $queryID);
	}
	else {
	    my $sequence_fasta_moby_object = "<$sequences_moby_output_format namespace='$namespace' id='$chromatogram_id'>\n<String id='' namespace='' articleName='content'><![CDATA[$fasta_sequences]]></String>\n</$sequences_moby_output_format>\n";
	    
	    my $quality_fasta_moby_object = "<$quality_moby_output_format namespace='$namespace' id='$chromatogram_id'>\n<String id='' namespace='' articleName='content'><![CDATA[$fasta_quality]]></String>\n</$quality_moby_output_format>\n";
	    
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_DOUBLE_SIMPLE_RESPONSE ($sequence_fasta_moby_object, $sequences_output_article_name, $quality_fasta_moby_object, $quality_output_article_name, $queryID);
	}
    }
    else {
	
	if ($_debug) {
	    print STDERR "no sequence returned by phred !\n";
	}
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_DOUBLE_SIMPLE_RESPONSE ($queryID, $sequences_output_article_name, $quality_output_article_name);
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

########################################################################

=head2 runPhred

 Title   : runPhred
 Usage   : Esta función está pensada para llamarla desde un cliente SOAP.
	 : No obstante, se recomienda probarla en la misma máquina, antes
	 : de instalar el servicio. Para ello, podemos llamarla de la
	 : siguiente forma:
	 :
	 : my $result = Phred("call", $in);
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

sub runPhred {
	
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
    
    my $_sequences_moby_output_format = "FASTA_NA";
    my $_quality_moby_output_format   = "FASTA_Base_Quality";
    my $moby_logger = get_logger ("MobyServices");
    my $serviceName = "runPhred";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput(@queries){

	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    # print STDERR "query text: $query_str\n";
	    my $DOM_substr = substr ($query_str, 0, 200);
	    print STDERR "DOM first lines, $DOM_substr\n";
	}

	# En este punto es importante recordar que el objeto $query
	# es un XML::DOM::Node, y que si queremos trabajar con
	# el mensaje de texto debemos llamar a: $query->toString()
	my ($query_response, $moby_exceptions_tmp) = _do_query_Phred ($queryInput, $_sequences_moby_output_format, $_quality_moby_output_format, "simple");
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

=head2 runPhredCollection

 Title   : runPhredCollection
 Usage   : Esta función está pensada para llamarla desde un cliente SOAP.
	 : No obstante, se recomienda probarla en la misma máquina, antes
	 : de instalar el servicio. Para ello, podemos llamarla de la
	 : siguiente forma:
	 :
	 : my $result = Phred("call", $in);
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
 Returns : Return a collection of  objects, one object for each input sequence object

=cut

sub runPhredCollection {

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
    
    my $_sequences_moby_output_format = "FASTA_NA_multi";
    my $_quality_moby_output_format   = "FASTA_Base_Quality_multi";
    my $moby_logger = get_logger ("MobyServices");
    my $serviceName = "runPhredCollection";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){

	# En este punto es importante recordar que el objeto $query
	# es un XML::DOM::Node, y que si queremos trabajar con
	# el mensaje de texto debemos llamar a: $query->toString()

	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    # print STDERR "query text: $query_str\n";
	    my $DOM_substr = substr ($query_str, 0, 100);
	    print STDERR "DOM first lines, $DOM_substr\n";
	}

	my ($query_response, $moby_exceptions_tmp) = _do_query_Phred ($queryInput, $_sequences_moby_output_format, $_quality_moby_output_format, "collection");
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
