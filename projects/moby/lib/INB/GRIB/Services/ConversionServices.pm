# $Id: ConversionServices.pm,v 1.2 2006-06-16 10:22:40 gmaster Exp $
#
# This file is an instance of a template written 
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#
# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::ConversionServices  - Package for conversion services

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

package INB::GRIB::Services::ConversionServices;

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
  &fromGenericSequenceToFASTA
  &fromGenericSequenceCollectionToFASTA
  &fromFASTAToDNASequence
  &fromFASTAToDNASequenceCollection
  &fromFASTAToAminoAcidSequence
  &fromFASTAToAminoAcidSequenceCollection
  &fromFASTAToGenericSequenceCollection
);

our $VERSION = '1.0';

my $_debug = 0;

# Preloaded methods go here.

###############################################################################


=head2 _do_query_fromGenericSequencestoFASTA

 Title   : _do_query_fromGenericSequencestoFASTA
         : 
         : private function (NOT EXPORTED)
         : 
 Usage   : my $query_response = _do_query_GeneID($query);
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

sub _do_query_fromGenericSequencestoFASTA {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby 
    my ($queryInput_DOM, $inputType) = @_;
    
    my $MOBY_RESPONSE   = "";     # set empty response
    my $moby_exceptions = [];
    my $output_article_name;
    if ($inputType eq "simple") {
	$output_article_name = "sequence";
    }
    elsif ($inputType eq "collection") {
	$output_article_name = "sequences";
    }
    
    my $seqobjs = [];
    
    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);
    my $namespace = "";
    
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
	
	if (($articleName =~ /sequence/i) || (isSimpleArticle ($DOM) || (isCollectionArticle ($DOM)))) { 
	    
	    if (isSimpleArticle ($DOM)) {
		
		if ($_debug) {
		    print STDERR "\"sequence\" tag is a simple article...\n";
		}
		
		if ($inputType eq "collection") {
		    my $note = "Received a simple input article instead of a collection input";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => $output_article_name,
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

		# Get the namespace
		
		($namespace) = INB::GRIB::Utils::CommonUtilsSubs->getNamespace ($DOM);
		
		if (!defined $namespace) {
		    print STDERR "namespace not defined\n";
		    $namespace = "";
		}
		
		# Validate the object type (GenericSequence)
		
		my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "GenericSequence");
		if ((! defined $rightType) || !$rightType) {
		    my $note = "Expecting a GenericSequence object, and receiving a $inputDataType object";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => $output_article_name,
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
		
                my $seqobj = INB::GRIB::Utils::CommonUtilsSubs->parseMobySequenceObjectFromDOMintoBioperlObject ($DOM);
		push (@$seqobjs, $seqobj);
		
	    }
	    elsif (isCollectionArticle ($DOM)) {
		
		if ($_debug) {
		    print STDERR "sequences is a collection article...\n";
		    # print STDERR "Collection DOM: " . $DOM->toString() . "\n";
		}
		
		if ($inputType eq "simple") {
		    my $note = "Received a collection input article instead of a simple input";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => $output_article_name,
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
		
		# Validate the object type (GenericSequence)
		
		my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "GenericSequence");
		if ((! defined $rightType) || !$rightType) {
		    my $note = "Expecting a GenericSequence object, and receiving a $inputDataType object";
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
		    
		    # Empty response
		    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_RESPONSE ($queryID, $output_article_name);
		    return ($MOBY_RESPONSE, $moby_exceptions);
		}
		
		# Get the namespace
		
		my @namespaces = INB::GRIB::Utils::CommonUtilsSubs->getNamespace ($DOM);
		if (@namespaces > 1) {
		    my $note = "Collection contains sequence objects belonging to more than one namespace. Returning a FASTA object with an empty namespace\n";
		    print STDERR "$note\n";
		    my $code = "201";
		    my $moby_exception = INB::Exceptions::MobyException->new (
									      refElement => 'sequences',
									      code       => $code,
									      type       => 'warning',
									      queryID    => $queryID,
									      message    => "$note",
									      );
		    push (@$moby_exceptions, $moby_exception);
		}
		else {
		    ($namespace) = @namespaces;
		}
		
		my @sequence_articles_DOM = getCollectedSimples ($DOM);
		
		foreach my $sequence_article_DOM (@sequence_articles_DOM) {
		    my $seqobj = INB::GRIB::Utils::CommonUtilsSubs->parseMobySequenceObjectFromDOMintoBioperlObject ($sequence_article_DOM);
		    push (@$seqobjs, $seqobj);
		}
		
	    }
	    else {
		print STDERR "It is not a simple or collection article...\n";
		print STDERR "DOM: " . $DOM->toString() . "\n";
	    }
	} # End parsing sequences article tag
	
    } # Next article
    
    # Check that we have parsed properly the sequences and the predictions
    
    if (@$seqobjs == 0) {
	my $note = "can't parse any sequences...\n";
	print STDERR "$note\n";
	my $code = "201";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => $output_article_name,
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_SIMPLE_RESPONSE ($queryID, $output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    # Una vez recogido todos los parametros necesarios, llamamos a 
    # la funcion que nos devuelve el report.
    
    my $fasta_sequences = convertSequencesIntoFASTA ($seqobjs);
    
    # Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX 
    # nos queda encapsularla en un Objeto bioMoby. Esta operacio 
    # la podriamos realizar en una funcion a parte si fuese compleja.  
    
    if (not defined $fasta_sequences) {
	my $note = "sequences conversion into FASTA format failed\n";
	print STDERR "$note\n";
	my $code = "201";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => $output_article_name,
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_SIMPLE_RESPONSE ($queryID, $output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    else {
	my $moby_output_format  = "FASTA";	
	my $sequenceIdentifier = "";
	# Specify the sequence identifier only in case there is only one sequence in the FASTA output
	if (@$seqobjs == 1) {
	    my $seqobj = $seqobjs->[0];
	    $sequenceIdentifier = $seqobjs->[0]->display_id;
	    
	    if (! defined $sequenceIdentifier) {
		print STDERR "sequence identifier is not defined!\n";
		print STDERR "might be a problem with the fasta sequences...\n";
		print STDERR "\n$fasta_sequences\n\n";
	    }
	    
	}
	
	my $fasta_object = <<PRT;
<moby:$moby_output_format namespace='$namespace' id='$sequenceIdentifier'>
<String namespace='' id='' articleName='content'>
<![CDATA[
$fasta_sequences
]]>
</String>
</moby:$moby_output_format>
PRT
 
        $MOBY_RESPONSE = simpleResponse($fasta_object, $output_article_name, $queryID);
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

=head2 _do_query_fromFASTAtoMobySequence

 Title   : _do_query_fromFASTAtoMobySequence
         : 
         : private function (NOT EXPORTED)
         : 
 Usage   : my $query_response = _do_query_GeneID($query);
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

sub _do_query_fromFASTAtoMobySequence {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby 
    my $queryInput_DOM          = shift @_;
    my $moby_input_object_type  = shift @_;
    my $moby_output_object_type = shift @_;
    my $input_article_name      = shift @_;
    my $output_article_name     = shift @_;
    
    if ($_debug) {
	print STDERR "moby_output_object_type, $moby_output_object_type\n";
    }
    
    my $MOBY_RESPONSE   = "";     # set empty response
    my $moby_exceptions = [];
    
    my $fasta_seqs;
    
    my $queryID   = getInputID ($queryInput_DOM);
    my @articles  = getArticles($queryInput_DOM);
    my $namespace = "";
    
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
	# Sino sabemos que es el input articulo porque es un simple articulo

	if (isCollectionArticle ($DOM)) {
	    # not allowed
	    
	    my $note = "Received a collection input article instead of a simple";
	    print STDERR "$note\n";
	    my $code = "201";
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      refElement => $input_article_name,
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    push (@$moby_exceptions, $moby_exception);
	    
	    # Return an empty moby data object, as well as an exception telling why nothing got returned
	    
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_SIMPLE_RESPONSE ($queryID, $output_article_name);
	    return ($MOBY_RESPONSE, $moby_exceptions);
	}
	
	# It's not very nice but taverna doesn't set up easily article name for input data so we let the users not setting up the article name of the input (which should be 'sequence')
	
	if (((defined $articleName) && ($articleName =~ /sequence$/i)) || isSimpleArticle ($DOM)) {
	    
	    # Get the namespace
	    
	    my @namespaces = INB::GRIB::Utils::CommonUtilsSubs->getNamespace ($DOM);
	    if (@namespaces > 1) {
		
		# Will never happen as it is a simple...
		
		my $note = "Collection contains sequence objects belonging to more than one namespace. Returning a sequence object with an empty namespace\n";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => $input_article_name,
									  code       => $code,
									  type       => 'warning',
									  queryID    => $queryID,
									  message    => "$note",
									  );
		push (@$moby_exceptions, $moby_exception);
	    }
	    else {
		($namespace) = @namespaces;
	    }
	    
	    if ($_debug) {
		print STDERR "\"sequences\" tag is a simple article...\n";
		print STDERR "stringified DOM, " . $DOM->toString () . "\n";
	    }
	    
	    # Validate the object type (FASTA)
	    
	    my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, $moby_input_object_type);
	    if ((! defined $rightType) || !$rightType) {
		my $note = "Expecting a $moby_input_object_type object, and receiving a $inputDataType object";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => $input_article_name,
									  code       => $code,
									  type       => 'error',
									  queryID    => $queryID,
									  message    => "$note",
									  );
		push (@$moby_exceptions, $moby_exception);
		
		# Empty response
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_SIMPLE_RESPONSE ($queryID, $output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
	    }
	    
	    # with this method, if the FASTA object is not compliant with latest specs, it won't be parsed
	    # ($fasta_seqs) = getNodeContentWithArticle($DOM, "String", "content");	
	    # Let's be backward compatible for the time being...
	    $fasta_seqs = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, $moby_input_object_type);
	    if ((not defined $fasta_seqs) || ($fasta_seqs eq "")) {
		$fasta_seqs = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "String");
	    }
	}
    } # Next article
    
    # Parse the FASTA sequences
    
    if ($_debug) {
	print STDERR "fasta sequences,\n$fasta_seqs.\n";
    }
    
    my $seqobjs   = INB::GRIB::Utils::CommonUtilsSubs->createSequenceObjectsFromFASTA ($fasta_seqs, $moby_output_object_type, $namespace);
    
    # Check that we have parsed properly the sequences
    
    if ((! defined $seqobjs) || (@$seqobjs != 1)) {
	my $note = "Can't parse any sequences, Please check the syntax of your fasta sequence.\n";
	print STDERR "$note\n";
	my $code = "201";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => $input_article_name,
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_SIMPLE_RESPONSE ($queryID, $output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    my $seqobj     = $seqobjs->[0];
    $MOBY_RESPONSE = simpleResponse($seqobj, $output_article_name, $queryID);

    return ($MOBY_RESPONSE, $moby_exceptions);
}


=head2 _do_query_fromFASTAtoMobySequences

 Title   : _do_query_fromFASTAtoMobySequences
         : 
         : private function (NOT EXPORTED)
         : 
 Usage   : my $query_response = _do_query_GeneID($query);
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
    
sub _do_query_fromFASTAtoMobySequences {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby 
    my $queryInput_DOM          = shift @_;
    my $moby_input_object_type  = shift @_;
    my $moby_output_object_type = shift @_;
    my $input_article_name      = shift @_;
    my $output_article_name     = shift @_;
    
    if ($_debug) {
	print STDERR "moby_output_object_type, $moby_output_object_type\n";
    }
    
    my $MOBY_RESPONSE   = "";     # set empty response
    my $moby_exceptions = [];
    
    my $fasta_seqs;
    
    my $queryID   = getInputID ($queryInput_DOM);
    my @articles  = getArticles($queryInput_DOM);
    my $namespace = "";
    
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
	# Sino sabemos que es el input articulo porque es un simple articulo

	if (isCollectionArticle ($DOM)) {
	    # not allowed
	    
	    my $note = "Received a collection input article instead of a simple";
	    print STDERR "$note\n";
	    my $code = "201";
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      refElement => $input_article_name,
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
	
	# It's not very nice but taverna doesn't set up easily article name for input data so we let the users not setting up the article name of the input (which should be 'input_article_name')

	if (((defined $articleName) && ($articleName =~ /sequence/i)) || isSimpleArticle ($DOM)) {

	    my @namespaces = INB::GRIB::Utils::CommonUtilsSubs->getNamespace ($DOM);
	    if (@namespaces > 1) {
		
		# Will never happen as it is a simple...
		
		my $note = "Collection contains sequence objects belonging to more than one namespace. Returning sequence objects with an empty namespace\n";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => $input_article_name,
									  code       => $code,
									  type       => 'warning',
									  queryID    => $queryID,
									  message    => "$note",
									  );
		push (@$moby_exceptions, $moby_exception);
	    }
	    else {
		($namespace) = @namespaces;
	    }
	    
	    if ($_debug) {
		print STDERR "\"sequences\" tag is a simple article...\n";
		print STDERR "stringified DOM, " . $DOM->toString () . "\n";
	    }
	    
	    # Validate the object type
	    
	    my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, $moby_input_object_type);
	    if ((! defined $rightType) || !$rightType) {
		my $note = "Expecting a $moby_input_object_type object, and receiving a $inputDataType object";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => $input_article_name,
									  code       => $code,
									  type       => 'error',
									  queryID    => $queryID,
									  message    => "$note",
									  );
		push (@$moby_exceptions, $moby_exception);
		
		# Empty response
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
	    }
	    
	    # with this method, if the FASTA object is not compliant with latest specs, it won't be parsed
	    # ($fasta_seqs) = getNodeContentWithArticle($DOM, "String", "content");	
	    # Let's be backward compatible for the time being...
	    $fasta_seqs = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, $moby_input_object_type);
	    if ((not defined $fasta_seqs) || ($fasta_seqs eq "")) {
		$fasta_seqs = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "String");
	    }
	}
    } # Next article
    
    # Parse the FASTA sequences
    
    if ($_debug) {
	print STDERR "fasta sequences,\n$fasta_seqs.\n";
    }
    
    my $seqobjs = INB::GRIB::Utils::CommonUtilsSubs->createSequenceObjectsFromFASTA ($fasta_seqs, $moby_output_object_type, $namespace);
    
    # Check that we have parsed properly the sequences
    
    if (@$seqobjs == 0) {
	my $note = "Can't parse any sequences, Please check the syntax of your fasta sequences.\n";
	print STDERR "$note\n";
	my $code = "201";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => $input_article_name,
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }

    $MOBY_RESPONSE = collectionResponse($seqobjs, $output_article_name, $queryID);

    return ($MOBY_RESPONSE, $moby_exceptions);
}


#################################
# Interface service methods
#################################


=head2 fromGenericSequenceToFASTA

 Title   : fromGenericSequenceToFASTA
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

sub fromGenericSequenceToFASTA {
    
    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_;        # get the incoming MOBY query XML

    if ($_debug) {
	print STDERR "processing Moby fromGenericSequenceToFASTA query...\n";
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
    my $moby_logger     = get_logger ("MobyServices");
    my $serviceName     = "fromGenericSequenceToFASTA";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){
	
	# En este punto es importante recordar que el objeto $query 
	# es un XML::DOM::Node, y que si queremos trabajar con 
	# el mensaje de texto debemos llamar a: $query->toString() 
	
	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}
	
	my ($query_response, $moby_exceptions_tmp) = _do_query_fromGenericSequencestoFASTA ($queryInput, "simple");
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

=head2 fromGenericSequenceCollectionToFASTA

 Title   : fromGenericSequenceCollectionToFASTA
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

sub fromGenericSequenceCollectionToFASTA {
    
    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_;        # get the incoming MOBY query XML

    if ($_debug) {
	print STDERR "processing Moby fromGenericSequenceCollectionToFASTA query...\n";
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
    my $moby_logger     = get_logger ("MobyServices");
    my $serviceName     = "fromGenericSequenceCollectionToFASTA";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){
	
	# En este punto es importante recordar que el objeto $query 
	# es un XML::DOM::Node, y que si queremos trabajar con 
	# el mensaje de texto debemos llamar a: $query->toString() 
	
	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}
	
	my ($query_response, $moby_exceptions_tmp) = _do_query_fromGenericSequencestoFASTA ($queryInput, "collection");
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


=head2 fromFASTAToDNASequence

 Title   : fromDNASequenceToFASTA
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

sub fromFASTAToDNASequence {
    
    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_;        # get the incoming MOBY query XML

    if ($_debug) {
	print STDERR "processing Moby fromFASTAToDNASequence query...\n";
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
    my $moby_logger     = get_logger ("MobyServices");
    my $serviceName     = "fromFASTAToDNASequence";

    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){
	
	# En este punto es importante recordar que el objeto $query 
	# es un XML::DOM::Node, y que si queremos trabajar con 
	# el mensaje de texto debemos llamar a: $query->toString() 
	
	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}
	
	my ($query_response, $moby_exceptions_tmp) = _do_query_fromFASTAtoMobySequence ($queryInput, "FASTA_NA", "DNASequence", "sequence", "sequence");
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

=head2 fromFASTAToDNASequenceCollection

 Title   : fromDNASequenceCollectionToFASTA
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

sub fromFASTAToDNASequenceCollection {
    
    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_;        # get the incoming MOBY query XML

    if ($_debug) {
	print STDERR "processing Moby fromFASTAToDNASequenceCollection query...\n";
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
    my $moby_logger     = get_logger ("MobyServices");
    my $serviceName     = "fromFASTAToDNASequenceCollection";

    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){
	
	# En este punto es importante recordar que el objeto $query 
	# es un XML::DOM::Node, y que si queremos trabajar con 
	# el mensaje de texto debemos llamar a: $query->toString() 
	
	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}
	
	my ($query_response, $moby_exceptions_tmp) = _do_query_fromFASTAtoMobySequences ($queryInput, "FASTA_NA_multi", "DNASequence", "sequences", "sequences");
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

=head2 fromFASTAToAminoAcidSequence

 Title   : fromAminoAcidSequenceToFASTA
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

sub fromFASTAToAminoAcidSequence {
    
    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_;        # get the incoming MOBY query XML

    if ($_debug) {
	print STDERR "processing Moby fromFASTAToAminoAcidSequence query...\n";
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
    my $moby_logger     = get_logger ("MobyServices");
    my $serviceName     = "fromFASTAToAminoAcidSequence";

    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){
	
	# En este punto es importante recordar que el objeto $query 
	# es un XML::DOM::Node, y que si queremos trabajar con 
	# el mensaje de texto debemos llamar a: $query->toString() 
	
	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}
	
	my ($query_response, $moby_exceptions_tmp) = _do_query_fromFASTAtoMobySequence ($queryInput, "FASTA_AA", "AminoAcidSequence", "sequence", "sequence");
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

=head2 fromFASTAToAminoAcidSequenceCollection

 Title   : fromAminoAcidSequenceCollectionToFASTA
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

sub fromFASTAToAminoAcidSequenceCollection {
    
    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_;        # get the incoming MOBY query XML

    if ($_debug) {
	print STDERR "processing Moby fromFASTAToAminoAcidSequenceCollection query...\n";
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
    my $moby_logger     = get_logger ("MobyServices");
    my $serviceName     = "fromFASTAToAminoAcidSequenceCollection";

    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){
	
	# En este punto es importante recordar que el objeto $query 
	# es un XML::DOM::Node, y que si queremos trabajar con 
	# el mensaje de texto debemos llamar a: $query->toString() 
	
	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}
	
	my ($query_response, $moby_exceptions_tmp) = _do_query_fromFASTAtoMobySequences ($queryInput, "FASTA_AA_multi", "AminoAcidSequence", "sequences", "sequences");
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



=head2 fromFASTAtoGenericSequenceCollection

 Title   : fromGenericSequenceCollectiontoFASTA
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

sub fromFASTAToGenericSequenceCollection {
    
    # El parametro $message es un texto xml con la peticion.
    my ($caller, $message) = @_;        # get the incoming MOBY query XML

    if ($_debug) {
	print STDERR "processing Moby fromFASTAToGenericSequenceCollection query...\n";
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
    my $moby_logger     = get_logger ("MobyServices");
    my $serviceName     = "fromFASTAToGenericSequenceCollection";
    
    # Para cada query ejecutaremos el _execute_query.
    foreach my $queryInput (@queries){
	
	# En este punto es importante recordar que el objeto $query 
	# es un XML::DOM::Node, y que si queremos trabajar con 
	# el mensaje de texto debemos llamar a: $query->toString() 
	
	if ($_debug) {
	    my $query_str = $queryInput->toString();
	    print STDERR "query text: $query_str\n";
	}
	
	my ($query_response, $moby_exceptions_tmp) = _do_query_fromFASTAtoMobySequences ($queryInput, "GenericSequence");
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

