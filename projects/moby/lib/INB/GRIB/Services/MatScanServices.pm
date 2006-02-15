# $Id: MatScanServices.pm,v 1.20 2006-02-15 16:49:18 gmaster Exp $
#
# This file is an instance of a template written
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#
# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::MatScanServices  - Package for parser the Moby message to call the MatScan promoter regions analysis program .

=head1 SYNOPSIS

Con este package podremos parsear las llamadas xml de BioMoby para
poder llamar al servicio runMatScanGFF. Una vez que tengamos la salida de llamando
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

package INB::GRIB::Services::MatScanServices;

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
  &runMatScanGFF
  &runMatScanGFFCollection
  &runMatScanGFFCollectionVsInputMatrix
);

our $VERSION = '1.0';

my $_debug = 0;

# Preloaded methods go here.

###############################################################################

=head2 _do_query_MatScan

 Title   : _do_query_MatScan
	 :
	 : private function (NOT EXPORTED)
	 :
 Usage   : my $query_response = _do_query_MatScan($query);
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

sub _do_query_MatScan {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby
    my $queryInput_DOM = shift @_;
    # $_format is the type of output that returns MatScan (e.g. GFF)
    my $_format        = shift @_;
    # $_input_type tells if the service specifies a collection of input objects or just a simple input object
    # We need this because the two services, runMatScanGFF and MatScanGFFCollection call both this method, so they share the same code
    # But the output type will have to be different
    my $_input_type     = shift @_;

    if ($_debug) {
      print STDERR "MatScan Sequences Vs predefined matrix...\n";
    }

    if (($_input_type ne "simple") && ($_input_type ne "collection")) {
	# Don't know the type, don't know what to return !
	print STDERR "_input_type unknown (should be setup as being simple or collection)\n";
	exit 0;
    }

    my $MOBY_RESPONSE   = "";     # set empty response
    my $moby_exceptions = [];
    my $output_article_name = "matscan_predictions";
    
    # Aqui escribimos las variables que necesitamos para la funcion.
    my $matrix;
    my $matrix_mode;
    my $threshold;
    my $strands;

    # Variables that will be passed to MatScan_call
    my %sequences;
    my %parameters;

    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);

    # Get the parameters

    ($matrix) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "matrix");
    if (not defined $matrix) {
	# Default is to use Transfac Matrices
	$matrix = "Transfac";
    }

    ($matrix_mode) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "matrix mode");
    if (not defined $matrix_mode) {
	# Default is to use 'log-likelihood' mode
	$matrix_mode = "log-likelihood";
    }

    ($threshold) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "threshold");
    if (not defined $threshold) {
	# Default is 0.85
	$threshold = "0.85";
    }

    ($strands) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "strands");
    if (not defined $strands) {
	# Default is running MatScan on both strands
	$strands = "Both";
    }

    # Add the parsed parameters in a hash table

    $parameters{threshold}   = $threshold;
    $parameters{strands}     = $strands;
    $parameters{matrix}      = $matrix;
    $parameters{matrix_mode} = $matrix_mode;

    if ($_debug) {
	print STDERR "matrix, $matrix\n";
	print STDERR "matrix mode, $matrix_mode\n";
	print STDERR "threshold, $threshold\n";
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
		
		if ($_debug) {
		    print STDERR "sequences tag is a simple article...\n";
		}
		
		if ($_input_type eq "collection") {
		    
		    # not allowed
		    
		    my $note = "Received a collection input article instead of a simple";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => "upstream_sequences",
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
		my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "NucleotideSequence");
		if (!$rightType) {
		    my $note = "Expecting a NucleotideSequence object, and receiving a $inputDataType object";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => "upstream_sequences",
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
		    print STDERR "sequences is a collection article...\n";
		    print STDERR "Collection DOM: " . $DOM->toString() . "\n";
		}
		
		if ($_input_type eq "simple") {
		    # is it allowed ? - check this out
		    # at the moment - disallow this !!
		    # Anyway the client is supposed to deal with this before the execution of the service
		    
		    my $note = "Received a simple input article instead of a collection";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => "upstream_sequences",
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
		
		# Validate the type of the simples in the collection - should all be NucleotideSequence objects
		my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "NucleotideSequence");
		if (!$rightType) {
		    my $note = "Expecting a NucleotideSequence object, and receiving a $inputDataType object";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => "upstream_sequences",
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
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	if ($_input_type eq "simple") {
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
	}
	else {
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_SIMPLE_RESPONSE ($queryID, $output_article_name);
	}
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    # Una vez recogido todos los parametros necesarios, llamamos a
    # la funcion que nos devuelve el report.
    
    if ($_input_type eq "simple") {
	
	if ($_debug) {
	    print STDERR "making a simple response\n";
	}
	
	my ($report, $moby_exceptions_tmp) = MatScan_call (sequences  => \%sequences, format => $_format, queryID => $queryID, parameters => \%parameters, debug => $_debug);
	push (@$moby_exceptions, @$moby_exceptions_tmp);
	
	# Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX
	# nos queda encapsularla en un Objeto bioMoby. Esta operacio
	# la podriamos realizar en una funcion a parte si fuese compleja.
	
	if (defined $report) {
	    my ($sequenceIdentifier) = keys (%sequences);
	    
	    my $input = <<PRT;
<moby:$_format namespace='' id='$sequenceIdentifier'>
<String namespace='' id='' articleName='content'>
<![CDATA[
$report
]]>
</String>
</moby:$_format>
PRT

            $MOBY_RESPONSE = simpleResponse($input, $output_article_name, $queryID);   
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
    elsif ($_input_type eq "collection") {

	if ($_debug) {
	    print STDERR "making a collection response\n";
	}

	# The input was a collection of Sequences, so we have to return a collection of GFF objects

	my ($report, $moby_exceptions_tmp) = MatScan_call (sequences  => \%sequences, format => $_format, parameters => \%parameters, debug => $_debug);
	push (@$moby_exceptions, @$moby_exceptions_tmp);

	if (defined $report) {
	    my $output_objects = INB::GRIB::Utils::CommonUtilsSubs->parseSingleGFFIntoCollectionGFF ($report, $_format, "");
	    $MOBY_RESPONSE = collectionResponse($output_objects, $output_article_name, $queryID);
	}
	else {
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


sub _do_query_MatScanVsInputMatrix {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby
    my $queryInput_DOM = shift @_;
    # $_format is the type of output that returns MatScan (e.g. GFF)
    my $_format        = shift @_;
    # $_input_type tells if the service specifies a collection of input objects or just a simple input object
    # We need this because the two services, runMatScanGFF and MatScanGFFCollection call both this method, so they share the same code
    # But the output type will have to be different
    my $_input_type     = shift @_;

    if ($_debug) {
	print STDERR "MatScan Sequences Vs Input Matrix mode...\n";
    }

    if (($_input_type ne "simple") && ($_input_type ne "collection")) {
	# Don't know the type, don't know what to return !
	print STDERR "_input_type unknown (should be setup as being simple or collection)\n";
	exit 1;
    }

    my $MOBY_RESPONSE   = "";     # set empty response
    my $moby_exceptions = [];
    my $output_article_name = "matscan_predictions";
    
    # Aqui escribimos las variables que necesitamos para la funcion.
    my $threshold;
    my $strands;
    my $matrix_mode;

    # Variables that will be passed to MatScan_call
    my %sequences;
    my $matrix;
    my %parameters;

    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);

    # Get the parameters

    ($matrix_mode) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "matrix mode");
    if (not defined $matrix_mode) {
	# Default is to use 'log-likelihood' mode
	$matrix_mode = "log-likelihood";
    }

    ($threshold) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "threshold");
    if (not defined $threshold) {
	# Default is 0.85
	$threshold = "0.85";
    }

    ($strands) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "strands");
    if (not defined $strands) {
	# Default is running MatScan on both strands
	$strands = "Both";
    }

    # Add the parsed parameters in a hash table

    $parameters{threshold}   = $threshold;
    $parameters{strands}     = $strands;
    $parameters{matrix_mode} = $matrix_mode;

    if ($_debug) {
	print STDERR "matrix mode, $matrix_mode\n";
    }

    # Tratamos a cada uno de los articulos
    foreach my $article (@articles) {

	# El articulo es una tupla que contiene el nombre de este
	# y su texto xml.

	my ($articleName, $DOM) = @{$article}; # get the named article
	
	# make it more 'interoperable' by testing also pattern matching!!!
	# i think it should be up to the clients to adjust the article names
	
	if (($articleName eq "upstream_sequences") || ($articleName =~ /sequences/i)) {
	    
	    if (isSimpleArticle ($DOM)) {

		# not allowed
		
		print STDERR "parsing the sequences input...\n";

		my $note = "Received a simple input article instead of a collection";
		print STDERR "$note\n";

		print STDERR "DOM: " . $DOM->toString() . "\n";

		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => "upstream_sequences",
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
	    elsif (isCollectionArticle ($DOM)) {

		if ($_debug) {
		    print STDERR "sequences is a collection article...\n";
		    print STDERR "Collection DOM: " . $DOM->toString() . "\n";
		}
		
		# Validate the type of the simples in the collection - should all be NucleotideSequence objects
		my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "NucleotideSequence");
		if (!$rightType) {
		    my $note = "Expecting a NucleotideSequence object, and receiving a $inputDataType object";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => "upstream_sequences",
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
	elsif (($articleName eq "matrices") || ($articleName =~ /matrices/i)) {
	    
	    $matrix = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "text-formatted");
	    
	    if ($_debug) {
		print STDERR "parsed matrix, $matrix\n";
	    }

	    if (not defined $matrix) {
		print STDERR "Error, can't parse the input matrix...\n";
	    }
	}

    } # Next article

    # Check that we have parsed properly the sequences

    if ((keys (%sequences)) == 0) {
	my $note = "can't parse any sequences...\n";
	print STDERR "$note\n";
	my $code = "201";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => 'upstream_sequences',
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    if ((! defined $matrix) || ($matrix eq "")) {
	my $note = "can't parse any matrix...\n";
	print STDERR "$note\n";
	my $code = "201";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => 'matrix',
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
    
    # The input was a collection of Sequences, so we have to return a collection of GFF objects
    
    my ($report, $moby_exceptions_tmp) = MatScan_call (sequences  => \%sequences, matrix => $matrix, format => $_format, parameters => \%parameters, debug => $_debug);
    push (@$moby_exceptions, @$moby_exceptions_tmp);

    if (defined $report) {
	my $output_objects = INB::GRIB::Utils::CommonUtilsSubs->parseSingleGFFIntoCollectionGFF ($report, $_format, "");
	$MOBY_RESPONSE = collectionResponse($output_objects, $output_article_name, $queryID);
    }
    else {
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


=head2 runMatScanGFF

 Title   : runMatScanGFF
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

sub runMatScanGFF {
	
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

    my $_moby_output_format   = "GFF";
    my $moby_logger = get_logger ("MobyServices");
    my $serviceName = "runMatScanGFF";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput(@queries){

	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}

	# En este punto es importante recordar que el objeto $query
	# es un XML::DOM::Node, y que si queremos trabajar con
	# el mensaje de texto debemos llamar a: $query->toString()
	my ($query_response, $moby_exceptions_tmp) = _do_query_MatScan ($queryInput, $_moby_output_format, "simple");
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

=head2 runMatScanGFFCollection

 Title   : runMatScanGFFCollection
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
 Returns : Return a collection of GFF objects, one object for each input sequence object

=cut

sub runMatScanGFFCollection {

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
    my $serviceName = "runMatScanGFFCollection";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){

	# En este punto es importante recordar que el objeto $query
	# es un XML::DOM::Node, y que si queremos trabajar con
	# el mensaje de texto debemos llamar a: $query->toString()

	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}

	my ($query_response, $moby_exceptions_tmp) = _do_query_MatScan ($queryInput, $_format, "collection");
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


sub runMatScanGFFCollectionVsInputMatrix {

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
    my $serviceName = "runMatScanGFFCollectionVsInputMatrix";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){

	# En este punto es importante recordar que el objeto $query
	# es un XML::DOM::Node, y que si queremos trabajar con
	# el mensaje de texto debemos llamar a: $query->toString()

	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}

	my ($query_response, $moby_exceptions_tmp) = _do_query_MatScanVsInputMatrix ($queryInput, $_format, "collection");
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
