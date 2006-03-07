# $Id: MemeServices.pm,v 1.23 2006-03-07 14:28:29 gmaster Exp $
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
    my $output_article_name = "meme_predictions";
    
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
    elsif (! (($motif_distribution eq "zero or one") || ($motif_distribution eq "one") || ($motif_distribution eq "any number of repetitions"))) {
	my $note = "motif distribution parameter, '$motif_distribution', not accepted, should be ['zero or one','one','any number of repetitions']";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "motif distribution",
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

    ($maximum_number_motifs) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "maximum number of motifs");
    if (not defined $maximum_number_motifs) {
	$maximum_number_motifs = 8;
    }
    elsif (($maximum_number_motifs > 12) || ($maximum_number_motifs < 0)) {
	my $note = "maximum number of motifs is out of range, [0,12]";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "maximum number motifs",
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
    
    ($minimum_number_sites) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "minimum sites for each motif");
    if (not defined $minimum_number_sites) {
	# no default
    }
    elsif (($minimum_number_sites < 2) ||  ($minimum_number_sites > 300)) {
	my $note = "minimum number of sites is out of range [2,300]";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "minimum number sites",
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
    
    ($maximum_number_sites) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "maximum sites for each motif");
    if (not defined $maximum_number_sites) {
	# no default
    }
    elsif (($maximum_number_sites < 2) ||  ($maximum_number_sites > 300)) {
	my $note = "maximum number of sites is out of range [2,300]";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "maximum number sites",
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
    
    if (defined $minimum_number_sites && defined $maximum_number_sites && $minimum_number_sites > $maximum_number_sites) {
	my $note = "Range incompatibility, you gave a minimum number of sites greater than the maximum one";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "minimum number sites",
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

    ($minimum_motif_width) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "minimum optimum width");
    if (not defined $minimum_motif_width) {
	$minimum_motif_width = 6;
    }
    elsif (($minimum_motif_width < 2) || ($minimum_motif_width > 300)) {
	my $note = "minimum motif width is out of range [2,300]";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "minimum motif width",
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

    ($maximum_motif_width) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "maximum optimum width");
    if (not defined $maximum_motif_width) {
	$maximum_motif_width = 15;
    }
    elsif (($maximum_motif_width < 2) || ($maximum_motif_width > 300)) {
	my $note = "maximum motif width is out of range [2,300]";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "maximum motif width",
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
    
    if ($minimum_motif_width > $maximum_motif_width) {
	my $note = "Range incompatibility, you gave a minimum number of sites greater than the maximum one";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "minimum motif width",
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
    
    ($e_value_cutoff) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "motif E-value cutoff");
    if (not defined $e_value_cutoff) {
	$e_value_cutoff = "1";
    }

    ($background_order) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "background markov model training (value is the model order)");
    if (not defined $background_order) {
	$background_order = 1;
    }
    elsif (($background_order =~ /\d+/) && (($background_order > 3) || ($background_order < 1))) {
	my $note = "background order, $background_order, is out of range [1,3]";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "background markov model training (value is the model order)",
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
    elsif ((! $background_order =~ /\d+/) && (lc ($background_order) ne "none")) {
	my $note = "background order, $background_order, doesn't match any expected value ['None', 1, 2, 3]";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "background markov model training (value is the model order)",
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

	if (isSimpleArticle ($DOM)) {
	    my $note = "Received a simple input article instead of a collection";
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
	    
	    # Return an empty moby data object, as well as an exception telling what nothing got returned
	    
	    $MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_RESPONSE ($queryID, $output_article_name);
	    return ($MOBY_RESPONSE, $moby_exceptions);
	}

	if (($articleName eq "sequences") || (isCollectionArticle ($DOM)) || ($articleName =~ /sequences/i)) {
	    
	    if ($_debug) {
		print STDERR "sequences is a collection article...\n";
		print STDERR "Collection DOM: " . $DOM->toString() . "\n";
	    }
	    
	    # Validate the type of the simples in the collection - should all be NucleotideSequence objects
	    my ($rightType, $inputDataType) = INB::GRIB::Utils::CommonUtilsSubs->validateDataType ($DOM, "GenericSequence");
	    if (!$rightType) {
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
		
		# Simple Response doesn't fit !! (the simple article is not empty as it should be!), so we need to create the string from scratch !
		$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_RESPONSE ($queryID, $output_article_name);
		return ($MOBY_RESPONSE, $moby_exceptions);
	    }
	    
	    my @sequence_articles_DOM = getCollectedSimples ($DOM);
	    
	    foreach my $sequence_article_DOM (@sequence_articles_DOM) {
		%sequences = INB::GRIB::Utils::CommonUtilsSubs->parseMobySequenceObjectFromDOM ($sequence_article_DOM, \%sequences);
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
	
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_RESPONSE ($queryID, $output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }

    # Una vez recogido todos los parametros necesarios, llamamos a
    # la funcion que nos devuelve el report.
    
    if ($_debug) {
	print STDERR "making a simple response\n";
    }

    my ($report, $moby_exceptions_tmp) = MEME_call (sequences  => \%sequences, format => $_output_format, debug => $_debug, queryID => $queryID, parameters => \%parameters);
    push (@$moby_exceptions, @$moby_exceptions_tmp);
    
    # Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX
    # nos queda encapsularla en un Objeto bioMoby. Esta operacio
    # la podriamos realizar en una funcion a parte si fuese compleja.
    
    if (defined $report) {
	my $input = <<PRT;
<moby:$_output_format namespace='' id=''>
<String namespace='' id='' articleName='content'>
<![CDATA[
$report
]]>
</String>
</moby:$_output_format>
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

sub _do_query_MemeMotifMatrices {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby
    my $queryInput_DOM = shift @_;
    my $_output_format = shift @_;

    my $MOBY_RESPONSE   = "";     # set empty response
    my $moby_exceptions = [];
    my $output_article_name = "meme_matrices";
    
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
    elsif (! ((lc ($matrix_mode) eq "raw format") || (lc ($matrix_mode) eq "log-likelihood"))) {
	my $note = "matrix mode parameter, '$matrix_mode', not accepted, should be ['raw format','log-likelihood']";
	print STDERR "$note\n";
	my $code = "222";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  refElement => "matrix mode",
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
    

    # Match the matrix object name and the type of its element depending on matrix mode, ie
    # raw format     => MatrixInteger
    # log-likelihood => MatrixFloat

    my $matrix_object_name;
    my $matrix_element_type;
    if ($matrix_mode eq "log-likelihood") {
	$matrix_object_name  = "MatrixFloat";
	$matrix_element_type = "Float";
    }
    elsif ($matrix_mode eq "raw format") {
	$matrix_object_name  = "MatrixInteger";
	$matrix_element_type = "Integer";
    }
    
    # Add the parsed parameters in a hash table
    
    $parameters{matrix_mode} = $matrix_mode;

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
	    
	    if (isCollectionArticle ($DOM)) {
		my $note = "Received a collection input article instead of a simple one";
		print STDERR "$note\n";
		my $code = "201";
		my $moby_exception = INB::Exceptions::MobyException->new (
									  refElement => "meme_predictions",
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

	    # Validation ...
	    
	    $meme_predictions = INB::GRIB::Utils::CommonUtilsSubs->getTextContentFromXML ($DOM, "MEME_Text");

	    if ($_debug) {
		print STDERR "parsed meme_predictions,\n$meme_predictions.\n";
	    }

	}
	if (isCollectionArticle ($DOM)) {
	    my $note = "Received a collection input article instead of a simple one";
	    print STDERR "$note\n";
	    my $code = "201";
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      refElement => "meme_predictions",
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
    my $namespace = "";
    
    if (not defined $matrices_aref) {
	# Return an emtpy message !
	$MOBY_RESPONSE = INB::GRIB::Utils::CommonUtilsSubs->MOBY_EMPTY_COLLECTION_RESPONSE ($queryID, $output_article_name);
	return ($MOBY_RESPONSE, $moby_exceptions);
    }
    
    my $meme_matrix_objects = [];    
    foreach my $matrix (@$matrices_aref) {
	my $moby_matrix = INB::GRIB::Utils::CommonUtilsSubs->convert_tabularPositionWeightMatrix_into_MobyMatrix ($matrix, $matrix_object_name, $matrix_element_type);
	push (@$meme_matrix_objects, $moby_matrix);
    }

    # Bien!!! ya tenemos el objeto de salida del servicio , solo nos queda
    # volver a encapsularlo en un objeto biomoby de respuesta. Pero
    # en este caso disponemos de una funcion que lo realiza. Si tuvieramos
    # una respuesta compleja (de verdad, esta era simple ;) llamariamos
    # a collection response.
    # IMPORTANTE: el identificador de la respuesta ($queryID) debe ser
    # el mismo que el de la query.
    
    $MOBY_RESPONSE = collectionResponse($meme_matrix_objects, $output_article_name, $queryID);
    return ($MOBY_RESPONSE, $moby_exceptions);
    
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
 Returns : Return a MEME_Text object with MEME results

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
    
    my $_format = "MEME_Text";
    my $moby_logger = get_logger ("MobyServices");
    my $serviceName = "runMemeText";
    
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
    
    my $response = INB::GRIB::Utils::CommonUtilsSubs->setMobyResponse ($MOBY_RESPONSE, $moby_exceptions, $moby_logger, $serviceName);
    
    return $response;

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
 Returns : Return a text_html object with MEME results

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
    # The moby output format for this service is text_html
    #
    my $_moby_output_format   = "text_html";
    my $moby_logger = get_logger ("MobyServices");
    my $serviceName = "runMemeHTML";
    
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
    
    my $response = INB::GRIB::Utils::CommonUtilsSubs->setMobyResponse ($MOBY_RESPONSE, $moby_exceptions, $moby_logger, $serviceName);
    
    return $response;
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
    
    my $_format = "MEME_Text";
    my $moby_logger = get_logger ("MobyServices");
    my $serviceName = "parseMotifMatricesfromMEME";
    
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
    
    my $response = INB::GRIB::Utils::CommonUtilsSubs->setMobyResponse ($MOBY_RESPONSE, $moby_exceptions, $moby_logger, $serviceName);
    
    return $response;

}


1;

__END__
