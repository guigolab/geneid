# $Id: Factory.pm,v 1.103 2006-07-28 14:56:15 gmaster Exp $
#
# INBPerl module for INB::GRIB::geneid::Factory
#
# This file is an instance of a template written
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#

# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::Factory - Package for calling geneid service.

=head1 SYNOPSIS

  # En este package tenemos las funciones necesarias para hacer las
  # llamadas al programa Geneid  La forma de utilizarlo en esta
  # versión es la siguiente:

 @params = ('nucleotide' => "ATCAGAGAGCGAGCAGATGCAGAGTAGAGCGAGCGAGCGCT");
 $report = INB::GRIB::Services::Factory::factory_call(@params);

  # Esta llamada nos devuelve una variable que contiene el texto con la
  # salida del programa geneid.

=head1 DESCRIPTION

Este package sirve para hacer la llamadas al programa geneid.
Se han utilizado las libreriras de geneid, por tanto serán
necesarias instalarlas previamente.

Es requisito necesario para ejecutar este módulo, que las variables de
entorno DIR y DATADIR estén descritas en el sistema. Dentro de
la variable $DIR se encuentra el path absoluto a la aplicacion geneid,
mientras que en $DATADIR pondremos la localización al directorio que
contiene las bases de datos.

En siguientes versiones, en este modulo se podran añadir las diferentes
funciones para llamar a los programas que se requieran.

=head1 AUTHOR

Arnaud Kerhornou, akerhornou@imim.es

=head1 COPYRIGHT

Copyright (c) 2004, Roman Roset Mayals and INB - Nodo Computacional UPC/CIRI.
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

package INB::GRIB::Services::Factory;

use strict;
use warnings;
use Data::Dumper;

use Carp;
# Just for Running GeneID - CGI call
use CGI;
use LWP::UserAgent;
use HTTP::Request;
use HTTP::Request::Common;

# Create temporary data files and directories
use File::Temp qw/tempfile/;
use File::Temp qw/tempdir/;

# MIME encoding/decoding
use MIME::Base64;

# Bioperl

use Bio::SeqIO;
use Bio::Seq;

# Moby Exceptions module
use INB::Exceptions::MobyException;

# Report Unix Error codes
use Errno qw (EINTR EIO :POSIX);

# Module to delete directories
use File::Path;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw() ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

# Aqui pondremos las funciones que exportaremos a otros módulos. Ejemplo
#
#  @EXPORT = qw( &func1 &func2);
#
our @EXPORT = qw(
  &GeneID_call_CGI
  &GeneID_call
  &SGP2_call
  &GOstat_call
  &TranslateGeneIDPredictions_call
  &PromoterExtraction_call
  &MatScan_call
  &Clover_call
  &MetaAlignment_call
  &MultiMetaAlignment_call
  &generateScoreMatrix_call
  &MEME_call
  &meme2matrix_call
  &RepeatMasker_call
  &Dust_call
  &Dust_FASTA_call
  &CrossMatchToScreenVector_call
  &Phrap_call
  &Phred_call
  &SequenceFilteringByLength_call
  &KMeans_call
  &GFF2JPEG_call
);

our $VERSION = '1.00';


###############################################################################
sub _exists_env {
	my $var_name = shift;
	croak "Environment variable \$$var_name not defined.\n"
	if (not exists $ENV{$var_name});
}

BEGIN {
	## in este array escribiremos las variables de entorno que le hacen
	## falta a la aplicacion para que funcione correctamente.
	## Ejemplo:  my @needed_vars = qw(DIR DATADIR);
	my @needed_vars = qw();
	_exists_env($_) foreach @needed_vars;
}

###############################################################################

=head2 GeneID_call_CGI

 Title   : GeneID_call_CGI
 Usage   : $report = GeneID_call_CGI (@params);
	 :
	 : ## where @params are,
	 : @params = ('arg1'  => "WKRPPEICENPRFIIGGANRTDIAAIACLTLNERL",
	 :            'arg2'  => "query 1", ## optional
	 :            'arg3'  => "nr");     ## optional (default: nr)
 Returns : Devuelve un string que contiene el resultado de la ejecución.

=cut

sub GeneID_call_CGI {
	my %args = @_;

	# relleno los parametros por defecto GeneID_call

	my $sequences          = $args{sequences}  || undef;
	my $format             = $args{format}     || "";
	my $parameters         = $args{parameters} || undef;
	
	# Get the parameters
	
	my $profile = $parameters->{profile};
	my $strands = $parameters->{strands};
	my $mode    = $parameters->{engine};
	
	# print STDERR "GeneID parameters (profile, strands, format, engine): $profile, $strands, $format, $engine.\n";
	
	# Output definition
	
	my $results         = "";
	my $moby_exceptions = [];
	
	# Parse the sequences hash
	
	my @seqIds = keys (%$sequences);
	
	foreach my $sequenceIdentifier (@seqIds) {
	    
	    my $nucleotides = $sequences->{$sequenceIdentifier};

	    # print STDERR "Submitting sequence, $sequenceIdentifier, of length, " . length ($nucleotides) . ", to GeneID CGI...\n";

	    # Llama a GeneID usando el interface web
	    my $diagplot_url = "http://genome.imim.es/cgi-bin/GeneID_cgi/geneid_2002/geneid_2002.cgi";
	    my $agent_diag   = LWP::UserAgent->new(timeout => 30000);
	    my $request_diag = POST($diagplot_url,
				    "Content-type" => "form-data",
				    "Content" => [
						  "seq"     => ">$sequenceIdentifier\n".$nucleotides,
						  "profile" => $profile,
						  "engines" => lc ($mode),
						  "strands" => lc ($strands),
						  "format"  => lc ($format)
						  ]);
	    my $result_diag = $agent_diag->request($request_diag);
	    $results       .= $result_diag->content;
	}
	
	return ($results, $moby_exceptions);
    }


=head2 GeneID_call

 Title   : GeneID_call
 Usage   : $report = GeneID_call (@params);
	 :
	 : ## where @params are,
	 : @params = ('arg1'  => "WKRPPEICENPRFIIGGANRTDIAAIACLTLNERL",
	 :            'arg2'  => "query 1", ## optional
	 :            'arg3'  => "nr");     ## optional (default: nr)
 Returns : Devuelve un string que contiene el resultado de la ejecución.

=cut

sub GeneID_call {
    my %args = @_;

    # output specs
    my $geneid_output   = "";
    my $peptides_output = "";
    my $moby_exceptions = [];

    # relleno los parametros por defecto GeneID_call

    my $sequences          = $args{sequences}  || undef;
    my $format             = $args{format}     || "";
    my $parameters         = $args{parameters} || undef;
    my $queryID            = $args{queryID}    || "";

    # Get the parameters

    my $profile     = $parameters->{profile};
    my $strands     = $parameters->{strands};
    my $exons_ref   = $parameters->{exons};
    my $signals_ref = $parameters->{signals};
    my $engine      = $parameters->{engine};
    
    # Llama a GeneID en local
    my $_geneid_dir  = "/home/ug/gmaster/projects/geneid";
    my $_geneid_bin  = "bin/geneid_latest_withGFF3output";
    my $_geneid_args = "";
    
    # Check that the binary is in place
    if (! -f "$_geneid_dir/$_geneid_bin") {
	my $note = "Geneid binary not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, undef, [$moby_exception]);
    }

    if (($format eq "GFF") || ($format eq "GFF2")) {
	$_geneid_args .= "-G";
    }
    
    if ($format eq "GFF3") {
	$_geneid_args .= "-3";
    }
    
    if ($strands eq "Both") {
	# Default anyway
    }
    elsif ($strands eq "Forward") {
	$_geneid_args .= "W";
    }
    elsif ($strands eq "Reverse") {
	$_geneid_args .= "C";
    }
    
    SWITCH: {
	if ($profile eq "Homo sapiens (suitable for mammals)")      { $_geneid_args .= "P $_geneid_dir/params/human.param"; last SWITCH; }
	if ($profile eq "Tetraodon nigroviridis (pupper fish)")     { $_geneid_args .= "P $_geneid_dir/params/tetraodon.param"; last SWITCH; }
	if ($profile eq "Drosophila melanogaster (fruit fly)")      { $_geneid_args .= "P $_geneid_dir/params/dros.param"; last SWITCH; }
	if ($profile eq "Caenorhabditis elegans (worm)") { $_geneid_args .= "P $_geneid_dir/params/celegans.param"; last SWITCH; }
	if ($profile eq "Triticum aestivum (wheat)")     { $_geneid_args .= "P $_geneid_dir/params/wheat.param"; last SWITCH; }
	if ($profile eq "Arabidopsis thaliana (weed)")   { $_geneid_args .= "P $_geneid_dir/params/arabidopsis.param"; last SWITCH; }
	if ($profile eq "Oryza sativa (rice)")           { $_geneid_args .= "P $_geneid_dir/params/rice.param"; last SWITCH; }
	if ($profile eq "Plasmodium falciparum (malaria parasite)") { $_geneid_args .= "P $_geneid_dir/params/plasmodium.param"; last SWITCH; }
	if ($profile eq "Dictyostelium discoideum (slime mold)")    { $_geneid_args .= "P $_geneid_dir/params/dictyostelium.param"; last SWITCH; }
	if ($profile eq "Aspergillus nidulans")          { $_geneid_args .= "P $_geneid_dir/params/aspergillus.param"; last SWITCH; }
	if ($profile eq "Neurospora crassa")             { $_geneid_args .= "P $_geneid_dir/params/neurospora.param"; last SWITCH; }
	if ($profile eq "Cryptococcus neomorfans")       { $_geneid_args .= "P $_geneid_dir/params/cneomorfans.param"; last SWITCH; }
	if ($profile eq "Coprinus cinereus")             { $_geneid_args .= "P $_geneid_dir/params/cinereus.param"; last SWITCH; }
	if ($profile eq "Apis mellifera (honey bee)")    { $_geneid_args .= "P $_geneid_dir/params/bee.param"; last SWITCH; }
	if ($profile eq "Chaetomium globosum")           { $_geneid_args .= "P $_geneid_dir/params/chaetomium.param"; last SWITCH; }
	if ($profile eq "Schistosoma japonica")          { $_geneid_args .= "P $_geneid_dir/params/sjaponica.param"; last SWITCH; }
	if ($profile eq "Stagnospora nodorum")           { $_geneid_args .= "P $_geneid_dir/params/snodorum.param"; last SWITCH; }
	if ($profile eq "Solanaceae")                    { $_geneid_args .= "P $_geneid_dir/params/solanaceae.param"; last SWITCH; }
	if ($profile eq "Sclerotinia sclerotiorum")      { $_geneid_args .= "P $_geneid_dir/params/ssclerotiorum.param"; last SWITCH; }
	if ($profile eq "Coccidioides immitis")          { $_geneid_args .= "P $_geneid_dir/params/cimmitis.param"; last SWITCH; }
	if ($profile eq "Histoplasma capsulatum")        { $_geneid_args .= "P $_geneid_dir/params/hcapsulatum.param"; last SWITCH; }
	# Default is Human
	$_geneid_args .= "P $_geneid_dir/params/human.param";
    }

    foreach my $exon (@$exons_ref) {
      SWITCH: {
	  if ($exon eq "None")                { last SWITCH; }
	  if ($exon eq "First exons")         { $_geneid_args .= " -f"; last SWITCH; }
	  if ($exon eq "Internal exons")      { $_geneid_args .= " -i"; last SWITCH; }
	  if ($exon eq "Terminal exons")      { $_geneid_args .= " -t"; last SWITCH; }
	  if ($exon eq "All exons")           { $_geneid_args .= " -x"; last SWITCH; }
	  if ($exon eq "Single genes")        { $_geneid_args .= " -s"; last SWITCH; }
	  if ($exon eq "Open reading frames") { $_geneid_args .= " -z"; last SWITCH; }
	  # Default is to leave blank the GeneID parameters line for not reporting Potential Exons Features
      }
    }

    foreach my $signal (@$signals_ref) {
      SWITCH: {
	  if ($signal eq "None")                  { last SWITCH; }
	  if ($signal eq "Donor splice sites")    { $_geneid_args .= " -d";  last SWITCH; }
	  if ($signal eq "Acceptor splice sites") { $_geneid_args .= " -a";  last SWITCH; }
	  if ($signal eq "All splice sites")      { $_geneid_args .= " -ad"; last SWITCH; }
	  if ($signal eq "Start codons")          { $_geneid_args .= " -b";  last SWITCH; }
	  if ($signal eq "Stop codons")           { $_geneid_args .= " -e";  last SWITCH; }
	  if ($signal eq "All codons")            { $_geneid_args .= " -be"; last SWITCH; }
	  if ($signal eq "All")                   { $_geneid_args .= " -adbe"; last SWITCH; }
	  # No signals reported by default
      }
    }
    
    SWITCH: {
	if ($engine eq "Normal")             { last SWITCH; }
	if ($engine eq "Exon Mode")          { $_geneid_args .= " -o"; last SWITCH; }
	# Default is Normal - nothing to add then
    }
    
    # Generate a temporary file locally with the sequence(s) in FASTA format
    # locally, ie not on a NFS mounted directory, for speed sake
    
    my ($seq_fh, $seqfile);
    eval {
	($seq_fh, $seqfile) = tempfile("/tmp/GENEID.XXXXXX", UNLINK => 0);
    };
    if ($@) {
	my $note = "Internal System Error. Can not open geneid input temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, undef, [$moby_exception]);
    }

    # Bioperl sequence factory

    my $sout = Bio::SeqIO->new (
				-fh     => $seq_fh,
				-format => 'fasta'
				);

    my @seqIds = keys (%$sequences);
    my $seqId = $seqIds[0];
    foreach my $sequenceIdentifier (@seqIds) {
	my $nucleotides = $sequences->{$sequenceIdentifier};
	# bioperl object
	my $seqobj = Bio::Seq->new (
				    -display_id => $sequenceIdentifier,
				    -seq        => $nucleotides
				    );
	$sout->write_seq ($seqobj);
    }
    close $seq_fh;
    
    # Test empty file
    if (-z $seqfile) {
	my $note = "Internal System Error. Empty geneid input sequence file...\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, undef, [$moby_exception]);
    }
    
    # print STDERR "Running GeneID, with this command:\n";
    # print STDERR "$_geneid_dir\/$_geneid_bin $_geneid_args $seqfile \n";
    
    $geneid_output = qx/$_geneid_dir\/$_geneid_bin $_geneid_args $seqfile/;
    
    # Comment this line if you want to keep the file...
    unlink $seqfile;
    
    if ((! defined $geneid_output) || (length $geneid_output < 1)) {
	my $note = "Internal System Error. Geneid has failed!\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	return (undef, undef, $moby_exceptions);
    }

    # Just for speed sake, we check the format but would work also with data in GFF2 format
    # runGeneIDGFF doesn't specify the peptides output so run this code just when it is runGeneIDGFF3 service
    
    if ($format eq "GFF3") {
	
	# Translate the gene predictions
	
	my %translation_parameters;
	$translation_parameters{translation_table} = "Standard (1)";
	my %geneid_predictions;
	$geneid_predictions{$seqId} = $geneid_output;
	
	my $moby_exceptions_tmp;
	($peptides_output, $moby_exceptions_tmp) = TranslateGeneIDPredictions_call (sequences => $sequences, predictions => \%geneid_predictions, parameters => \%translation_parameters);
	push (@$moby_exceptions, @$moby_exceptions_tmp);
	
	
	if (defined $peptides_output && (length $peptides_output > 0)) {
	    
	    # Cleaning the mRNA identifier in the headers
	    $peptides_output =~ s/Parent=//g;
	    $peptides_output =~ s/;ExonType=\w+//g;
	    
	    return ($geneid_output, $peptides_output, $moby_exceptions);
	}
	else {
	    my $note = "Internal System Error. Geneid has failed!\n";
	    print STDERR "$note\n";
	    my $code = 701;
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    push (@$moby_exceptions, $moby_exception);
	    return (undef, undef, $moby_exceptions);
	}
    }
    else {
	return ($geneid_output, undef, $moby_exceptions);
    }
}

=head2 SGP2_call

 Title   : SGP2_call
 Usage   : $report = SGP2_call (@params);
	 :
	 : ## where @params are,
	 : @params = ('arg1'  => "WKRPPEICENPRFIIGGANRTDIAAIACLTLNERL",
	 :            'arg2'  => "query 1", ## optional
	 :            'arg3'  => "nr");     ## optional (default: nr)
 Returns : Devuelve un string que contiene el resultado de la ejecución.

=cut

sub SGP2_call {
	my %args = @_;

	# output specs declaration
	my $sgp2_output     = "";
	my $moby_exceptions = [];

	# relleno los parametros por defecto SGP2_call (nucleotide  => $nucleotide, seqIdentifier => $sequenceIdentifier);
	my $sequences          = $args{sequences}      || undef;
	my $tblastx_output     = $args{tblastx_output} || undef;
	my $format             = $args{format}         || "";
	my $parameters         = $args{parameters}     || undef;
	my $queryID            = $args{queryID}        || "";
	my $debug              = $args{debug}          || 0;
	# Get the parameters
	
	my $profile = $parameters->{profile};
	# print STDERR "GeneID parameters (profile, format): $profile, $format.\n";
	
	# Llama a SGP2 localmente
	my $_sgp2_dir  = "/home/ug/gmaster/projects/sgp2/";
	$_sgp2_dir     = $ENV{SGP2};
	my $_sgp2_bin  = "bin/sgp2";
	my $_sgp2_args = "";
	
	# Check that the binary is in place
	if (! -f "$_sgp2_dir/$_sgp2_bin") {
	    my $note = "Internal System Error. SGP2 binary not found";
	    print STDERR "$note\n";
	    my $code = 701;
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    return (undef, [$moby_exception]);
	}
	
	if ($format eq "GFF") {
	    $_sgp2_args .= "-g G";
	}
	
        SWITCH: {
	    if ($profile eq "Human Vs Mouse")      { $_sgp2_args .= " -P $_sgp2_dir/params/human3iso.sgp.Hs-Mm.param"; last SWITCH; }
	    if ($profile eq "Human Vs Chicken")     { $_sgp2_args .= " -P $_sgp2_dir/params/human3iso.sgp.Gg-Hs.param"; last SWITCH; }
	    # Default is Human Vs Mouse
	    $_sgp2_args .= " -P $_sgp2_dir/params/human3iso.sgp.Hs-Mm.param";
	}
	
	# Generate a temporary file locally with the sequence in FASTA format
	# locally, ie not on a NFS mounted directory, for speed sake

	my ($seq_fh, $seqfile);
	eval {
	    ($seq_fh, $seqfile) = tempfile("/tmp/SGP2_Sequence.XXXXXX", UNLINK => 0);
	};
	if ($@) {
	    my $note = "Internal System Error. Can not open SGP2 sequence input temporary file!\n";
	    my $code = 701;
	    print STDERR "$note\n";
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    return (undef, [$moby_exception]);
	}
	
	# Bioperl sequence factory
	
	my $sout = Bio::SeqIO->new (
				    -fh     => $seq_fh,
				    -format => 'fasta'
				    );

	my @seqIds = keys (%$sequences);
	foreach my $sequenceIdentifier (@seqIds) {
	    my $nucleotides = $sequences->{$sequenceIdentifier};
	    
	    # bioperl sequence object
	    
	    my $seqobj = Bio::Seq->new (
					-display_id => $sequenceIdentifier,
					-seq        => $nucleotides
					);
	    $sout->write_seq ($seqobj);
	}
	close $seq_fh;

	# TBLASTX Output File

	# Generate a temporary file locally with the TBLASTX Output
	# locally, ie not on a NFS mounted directory, for speed sake

	my ($blast_fh, $tblastx_output_file);
	eval {
	    ($blast_fh, $tblastx_output_file) = tempfile("/tmp/SGP2_TBLASTX.XXXXXX", UNLINK => 0);
	    print $blast_fh "$tblastx_output";
	    close $blast_fh;
	};
	if ($@) {
	    my $note = "Internal System Error. Can not open SGP2 tblastx input temporary file!\n";
	    my $code = 701;
	    print STDERR "$note\n";
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    return (undef, [$moby_exception]);
	}

	# Test empty files

	if (-z $seqfile) {
	    my $note = "Internal System Error. Empty SGP2 input sequence file...\n";
	    print STDERR "$note\n";
	    my $code = 701;
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    return (undef, [$moby_exception]);
	}

	if (-z $tblastx_output_file) {
	    my $note = "Internal System Error. Empty SGP2 input tblastx file...\n";
	    print STDERR "$note\n";
	    my $code = 701;
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    return (undef, [$moby_exception]);
	}
	
	if ($debug) {
	    print STDERR "Running SGP2, with this command:\n";
	    print STDERR "$_sgp2_dir\/$_sgp2_bin $_sgp2_args -1 $seqfile -t $tblastx_output_file\n";
	}
	
	$sgp2_output = qx/$_sgp2_dir\/$_sgp2_bin $_sgp2_args -1 $seqfile -t $tblastx_output_file/;
	
	if (! $debug) {
	    unlink $seqfile;
	    unlink $tblastx_output_file;
	}
	
	if ((defined $sgp2_output) && (length ($sgp2_output) > 1)) {
	    return ($sgp2_output, $moby_exceptions);
	}
	else {
	    my $note = "Internal System Error. SGP2 has failed!\n";
	    print STDERR "$note\n";
	    
	    # print STDERR "Running SGP2, with this command:\n";
	    print STDERR "$_sgp2_dir\/$_sgp2_bin $_sgp2_args -1 $seqfile -t $tblastx_output_file\n";
	    
	    my $code = 701;
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    push (@$moby_exceptions, $moby_exception);
	    return (undef, $moby_exceptions);
	}
}


=head2 GOstat_call

 Title   : GOstat_call
 Usage   : $report = GOstat_call (@params);
	 :
	 : ## where @params are,
	 : @params = ('arg1'  => "WKRPPEICENPRFIIGGANRTDIAAIACLTLNERL",
	 :            'arg2'  => "query 1", ## optional
	 :            'arg3'  => "nr");     ## optional (default: nr)
 Returns : Devuelve un string que contiene el resultado de la ejecución.

=cut

sub GOstat_call {
    my %args = @_;

    # output specs declaration
    my $gostat_output   = "";
    my $moby_exceptions = [];

    my $regulated_genes    = $args{regulated_genes} || undef;
    my $reference_genes    = $args{reference_genes} || undef;
    my $format             = $args{format}          || "";
    my $parameters         = $args{parameters}      || undef;
    my $queryID            = $args{queryID}         || "";

    # Llama a GOstat localmente
    my $_gostat_dir  = "/home/ug/gmaster/projects/gostat";
    my $_gostat_bin  = "gostat.pl";
    my $_gostat_args = "";

    # Make two temporary files for both input lists of genes

    my ($regulated_genes_fh, $regulated_genes_file) = tempfile("/tmp/REGULATED_GENES.XXXXXX", UNLINK => 0);
    close ($regulated_genes_fh);
    my ($reference_genes_fh, $reference_genes_file) = tempfile("/tmp/REFERENCE_GENES.XXXXXX", UNLINK => 0);
    close ($reference_genes_fh);

    open (FILE, ">$regulated_genes_file") or die "can't open temp file, $regulated_genes_file!\n";
    print FILE (join ("\n", @$regulated_genes) . "\n");
    close FILE;

    open (FILE, ">$reference_genes_file") or die "can't open temp file, $reference_genes_file!\n";
    print FILE (join ("\n", @$reference_genes) . "\n");
    close FILE;

    $gostat_output = qx/$_gostat_dir\/$_gostat_bin $_gostat_args --reg $regulated_genes_file --ref $reference_genes_file/;

    unlink $regulated_genes_file;
    unlink $reference_genes_file;

    if (defined $gostat_output) {
	return $gostat_output;
    }
    else {
	# What else better to return ??
	return undef;
    }

}

sub TranslateGeneIDPredictions_call {
    my %args = @_;

    # output specs declaration
    my $translateGeneID_output = "";
    my $moby_exceptions        = [];

    my $sequences          = $args{sequences}    || undef;
    my $geneid_predictions = $args{predictions}  || undef;
    my $parameters         = $args{parameters}   || undef;
    my $queryID            = $args{queryID}      || "";

    my $translation_table  = $parameters->{translation_table};
    my $translation_code;
    SWITCH: {
	if ($translation_table eq "Standard (1)")   { $translation_code = 1; last SWITCH; }
	if ($translation_table eq "bacterial (11)") { $translation_code = 11; last SWITCH; }
	# Default is Standard
	$translation_code = 1;
    }

    # Llama a GOstat localmente
    my $_translateGeneID_dir  = "/home/ug/gmaster/projects/GFF_Translations";
    my $_translateGeneID_bin  = "GFF_Features_Translation.pl";
    my $_translateGeneID_args = "-t $translation_code";

    # Check that the binary is in place
    if (! -f "$_translateGeneID_dir/$_translateGeneID_bin") {
	my $note = "Internal System Error. GFF_Translations binary not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }

    # Make two temporary files for both input lists of genes
    
    my ($seq_fh, $seqfile);
    eval {
	($seq_fh, $seqfile) = tempfile("/tmp/SEQS.TRANSLATION.XXXXXX", UNLINK => 0);
    };
    my ($feature_fh, $featurefile);
    eval {
	($feature_fh, $featurefile) = tempfile("/tmp/FEATURES.TRANSLATION.XXXXXX", UNLINK => 0);
    };

    # Bioperl sequence factory
    my $sout = Bio::SeqIO->new (
				-fh     => $seq_fh,
				-format => 'fasta'
				);

    my @seqIds = keys (%$sequences);
    foreach my $sequenceIdentifier (@seqIds) {
	# Sequence
	my $nucleotides = $sequences->{$sequenceIdentifier};
	
	# bioperl sequence object
	my $seqobj = Bio::Seq->new (
				    -display_id => $sequenceIdentifier,
				    -seq        => $nucleotides
				    );
	$sout->write_seq ($seqobj);
	
	# GeneID Predictions
	my $geneid_prediction = $geneid_predictions->{$sequenceIdentifier};
	print $feature_fh "$geneid_prediction";
	
    }
    close $seq_fh;
    close $feature_fh;

    # Test empty file
    if (-z $seqfile) {
	my $note = "Internal System Error when calling translateGeneIDPredictions sequence. Empty sequence file...\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }

    # print STDERR "Running the following command, $_translateGeneID_dir\/$_translateGeneID_bin $_translateGeneID_args -s $seqfile -f $featurefile...\n";
    
    $translateGeneID_output = qx/$_translateGeneID_dir\/$_translateGeneID_bin $_translateGeneID_args -s $seqfile -f $featurefile/;
    
    if ((defined $translateGeneID_output) && ($translateGeneID_output ne "")) {
	# Comment these two lines if you want to keep those files...
	unlink $seqfile;
	unlink $featurefile;
	
	return ($translateGeneID_output, $moby_exceptions);
    }
    else {
	my $note = "Internal System Error. The translation of GeneID GFF exon features has failed!!\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	
	print STDERR "Ran the following command, $_translateGeneID_dir\/$_translateGeneID_bin $_translateGeneID_args -s $seqfile -f $featurefile...\n";
	# Let's keep the files for debugging purposes
	# unlink $seqfile;
	# unlink $featurefile;

	return (undef, $moby_exceptions);
    }

}

sub PromoterExtraction_call {
    my %args = @_;

    # output specs declaration
    my $promoterExtraction_output = "";
    my $moby_exceptions           = [];

    # relleno los parametros por defecto PromoterExtraction_call
    
    my $genes_ref  = $args{genes}      || undef;
    my $parameters = $args{parameters} || undef;
    my $queryID    = $args{queryID}    || "";
    my $debug      = $args{debug}      || 0;
    
    # Get the parameters
    
    my $organism          = $parameters->{organism};
    my $dbrelease         = $parameters->{dbrelease};
    my $upstream_length   = $parameters->{upstream_length};
    my $downstream_length = $parameters->{downstream_length};
    my $intergenic_only   = $parameters->{intergenic_only};
    my $orthologous_mode  = $parameters->{orthologous_mode};
    
    # Llama a GeneID en local
    my $_promExtraction_dir  = "/home/ug/gmaster/projects/promoter_extraction";
    my $_promExtraction_bin  = "promoter_extraction.pl";
    my $_promExtraction_args = "";
    
    # Check that the binary is in place
    if (! -f "$_promExtraction_dir/$_promExtraction_bin") {
	my $note = "Internal System Error. promoter_extraction.pl script not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    if (lc ($intergenic_only) eq "true") {
	$_promExtraction_args .= "-i";
    }
    if (lc ($orthologous_mode) eq "true") {
	$_promExtraction_args .= " -o";
    }
    
    $_promExtraction_args .= " -u $upstream_length";
    $_promExtraction_args .= " -d $downstream_length";
    
    SWITCH: {
	if (lc ($organism) eq "homo sapiens")             { $_promExtraction_args .= " -s homo_sapiens"; last SWITCH; }
	if (lc ($organism) eq "pan troglodytes")          { $_promExtraction_args .= " -s pan_troglodytes"; last SWITCH; }
	if (lc ($organism) eq "mus musculus")             { $_promExtraction_args .= " -s mus_musculus"; last SWITCH; }
	if (lc ($organism) eq "rattus norvegicus")        { $_promExtraction_args .= " -s rattus_norvegicus"; last SWITCH; }
	if (lc ($organism) eq "canis familiaris")         { $_promExtraction_args .= " -s canis_familiaris"; last SWITCH; }
	if (lc ($organism) eq "bos taurus")               { $_promExtraction_args .= " -s bos_taurus"; last SWITCH; }
	if (lc ($organism) eq "monodelphis domestica")               { $_promExtraction_args .= " -s monodelphis-domestica"; last SWITCH; }
	if (lc ($organism) eq "gallus gallus")            { $_promExtraction_args .= " -s gallus_gallus"; last SWITCH; }
	if (lc ($organism) eq "xenopus tropicalis")       { $_promExtraction_args .= " -s xenopus_tropicalis"; last SWITCH; }
	if (lc ($organism) eq "danio rerio")              { $_promExtraction_args .= " -s danio_rerio"; last SWITCH; }
	if (lc ($organism) eq "takifugu rubripes")        { $_promExtraction_args .= " -s takifugu_rubripes"; last SWITCH; }
	if (lc ($organism) eq "tetraodon nigroviridis")   { $_promExtraction_args .= " -s tetraodon_nigroviridis"; last SWITCH; }
	if (lc ($organism) eq "ciona intestinalis")       { $_promExtraction_args .= " -s ciona_intestinalis"; last SWITCH; }
	if (lc ($organism) eq "drosophila melanogaster")  { $_promExtraction_args .= " -s drosophila_melanogaster"; last SWITCH; }
	if (lc ($organism) eq "anopheles gambiae")        { $_promExtraction_args .= " -s anopheles_gambiae"; last SWITCH; }
	if (lc ($organism) eq "apis mellifera")           { $_promExtraction_args .= " -s apis_mellifera"; last SWITCH; }
	if (lc ($organism) eq "caenorhabditis elegans")   { $_promExtraction_args .= " -s caenorhabditis_elegans"; last SWITCH; }
	if (lc ($organism) eq "saccharomyces cerevisiae") { $_promExtraction_args .= " -s saccharomyces_cerevisiae"; last SWITCH; }
	# Default is Human
	$_promExtraction_args .= " -s homo_sapiens";
    }
    
    $_promExtraction_args .= " -r $dbrelease";
    
    # Make a temporary file for the input list of genes
    
    my ($genes_list_fh, $genes_list_file);
    eval {
	($genes_list_fh, $genes_list_file) = tempfile("/tmp/PROM_EXTRACTION_GENES.XXXXXX", UNLINK => 1);
	print $genes_list_fh (join ("\n", @$genes_ref) . "\n");
    };
    if ($@) {
	close $genes_list_fh;
	my $note = "Internal System Error. Can not open gene identifier list input temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    close $genes_list_fh;
    
    if ($debug) {
	print STDERR "running command,\n";
	print STDERR "$_promExtraction_dir\/$_promExtraction_bin $_promExtraction_args -f $genes_list_file\n";
    }
    
    $promoterExtraction_output = qx/$_promExtraction_dir\/$_promExtraction_bin $_promExtraction_args -f $genes_list_file/;
    
    if (!$debug) {
	unlink $genes_list_file;
    }
    
    if (defined $promoterExtraction_output) {
	return ($promoterExtraction_output, $moby_exceptions);
    }
    else {
	my $note = "Internal System Error. prmoter extraction script execution has failed!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
}


=head2 MatScan_call

 Title   : MatScan_call
 Usage   : $report = MatScan_call (@params);
	 :
	 : ## where @params are,
	 : @params = ('arg1'  => "WKRPPEICENPRFIIGGANRTDIAAIACLTLNERL",
	 :            'arg2'  => "query 1", ## optional
	 :            'arg3'  => "nr");     ## optional (default: nr)
 Returns : Devuelve un string que contiene el resultado de la ejecución.

=cut

sub MatScan_call {
    my %args = @_;

    # output specs declaration
    my $matscan_output   = "";
    my $moby_exceptions  = [];

    # relleno los parametros por defecto MatScan_call

    my $sequences      = $args{sequences}  || undef;
    my $matrices_input = $args{matrices}   || undef;
    my $format         = $args{format}     || "";
    my $parameters     = $args{parameters} || undef;
    my $debug          = $args{debug}      || 0;
    my $queryID        = $args{queryID}    || "";

    # Get the parameters

    my $threshold          = $parameters->{threshold};
    my $strands            = $parameters->{strands};
    my $database_parameter = $parameters->{motif_database};
    my $matrix_mode        = $parameters->{matrix_mode};
    
    if ($debug) {
	print STDERR "threshold, $threshold\n";
	print STDERR "motif database parameter, $database_parameter\n";
	print STDERR "matrix mode, $matrix_mode\n";
    }
    
    # Llama a MatScan en local
    my $_matscan_dir  = "/home/ug/gmaster/projects/Meta/";
    my $_matscan_bin  = "bin/matscan";
    my $_matscan_args = "-T $threshold";
    # Check that the binary is in place
    if (! -f "$_matscan_dir/$_matscan_bin") {
	my $note = "Internal System Error. Matscan binary not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }

    my ($matrix_fh, $matrix_file);

    if ($strands eq "Both") {
	# Default anyway
    }
    elsif ($strands eq "Forward") {
	$_matscan_args .= " -W";
    }
    elsif ($strands eq "Reverse") {
	$_matscan_args .= " -C";
    }

    if (defined $database_parameter) {

	if ($debug) {
	    print STDERR "matrix as a parameter...\n";
	}

	if ($matrix_mode eq "raw format") {

	    if ($debug) {
		print STDERR "raw mode\n";
	    }

	  SWITCH: {
	      if (lc ($database_parameter) eq "transfac") { $matrix_file = "$_matscan_dir/matrices/Transfac_raw_format.matrices"; last SWITCH; }
	      if (lc ($database_parameter) eq "jaspar")   { $matrix_file = "$_matscan_dir/matrices/Jaspar_raw_format.matrices"; last SWITCH; }
	      # Default is Transfac
	      $matrix_file = "$_matscan_dir/matrices/Transfac_raw_format.matrices";
	  }
	}
	elsif ($matrix_mode eq "log-likelihood") {
	    if ($debug) {
		print STDERR "log-likelihood mode\n";
	    }

	  SWITCH: {
	      if (lc ($database_parameter) eq "transfac") { $_matscan_args .= " -sm"; $matrix_file = "$_matscan_dir/matrices/Transfac_likelihood.matrices"; last SWITCH; }
	      if (lc ($database_parameter) eq "jaspar")   { $_matscan_args .= " -s"; $matrix_file = "$_matscan_dir/matrices/Jaspar_likelihood.matrices"; last SWITCH; }
	      # Default is Transfac
	      $_matscan_args .= " -sm";
	      $matrix_file = "$_matscan_dir/matrices/Transfac_likelihood.matrices";
	  }
	}
	else {
	    # should be validated before...
	    print STDERR "don't know anything about matrix mode, $matrix_mode!\n";
	    exit 0;
	}
    }
    elsif (defined $matrices_input) {
	if ($debug) {
	    print STDERR "matrices as an input...\n";
	}

	# Make a temporary file with the matrix input
	$_matscan_args .= " -s";
	eval {
	    ($matrix_fh, $matrix_file) = tempfile("/tmp/MATSCAN_MATRIX.XXXXXX", UNLINK => 0);
	    print $matrix_fh "$matrices_input";
	    close $matrix_fh;
	};
	if ($@) {
	    my $note = "Internal System Error. Can not open MatScan matrix input temporary file!\n";
	    my $code = 701;
	    print STDERR "$note\n";
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    return (undef, [$moby_exception]);
	}

    }
    else {
	# To change !
	print STDERR "matrices_input neither database_parameter are defined!!\n";
	exit 0;
    }

    if ((not defined $matrix_file) || (-z $matrix_file)) {
	# could well be possible if the input matrix set is empty, ie if MEME didn't predict any !!
	# But i guess in that case, no need to go up there, validation will be done before MatScan_call
	print STDERR "Internal System Error. No defined matrix file!\n";
	return (undef, []);
    }

    # Generate a temporary file locally with the sequence(s) in FASTA format
    # locally, ie not on a NFS mounted directory, for speed sake

    my ($seq_fh, $seqfile);
    eval {
	($seq_fh, $seqfile) = tempfile("/tmp/MATSCAN_SEQS.XXXXXX", UNLINK => 0);
    };
    if ($@) {
	my $note = "Internal System Error. Can not open MatScan sequence input temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }

    # Bioperl sequence factory

    my $sout = Bio::SeqIO->new (
				-fh     => $seq_fh,
				-format => 'fasta'
				);

    my @seqIds = keys (%$sequences);
    foreach my $sequenceIdentifier (@seqIds) {
	my $nucleotides = $sequences->{$sequenceIdentifier};

	# bioperl sequence object

	my $seqobj = Bio::Seq->new (
				    -display_id => $sequenceIdentifier,
				    -seq        => $nucleotides
				    );
	$sout->write_seq ($seqobj);
    }
    close $seq_fh;

    # Test empty file
    if (-z $seqfile) {
	my $note = "Internal System Error. Empty MatScan input sequence file...\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }

    if ($debug) {
	print STDERR "Running Matscan, with this command:\n";
	print STDERR "$_matscan_dir\/$_matscan_bin $_matscan_args $seqfile $matrix_file\n";
    }
    
    $matscan_output = qx/$_matscan_dir\/$_matscan_bin $_matscan_args $seqfile $matrix_file | grep MatScan | sort +3n/;
    
    unlink $seqfile unless $debug;
    if ((not $debug) && (defined $matrices_input)) {
	# Only remove the matrix input file when this matrix collection is given by the user!!
	# If it is transfac or jaspar, don't remove it !!! - as they are stored locally !
	unlink $matrix_file;
    }

    if (defined $matscan_output) {
	return ($matscan_output, $moby_exceptions);
    }
    else {
	my $note = "Internal System Error. MatScan has failed!\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	return (undef, $moby_exceptions);
    }
}

=head2 Clover_call

 Title   : Clover_call
 Usage   : $report = Clover_call (@params);
	 :
	 : ## where @params are,
	 : @params = ('arg1'  => "WKRPPEICENPRFIIGGANRTDIAAIACLTLNERL",
	 :            'arg2'  => "query 1", ## optional
	 :            'arg3'  => "nr");     ## optional (default: nr)
 Returns : Devuelve un string que contiene el resultado de la ejecución.

=cut

sub Clover_call {
    my %args = @_;

    # output specs declaration
    my $clover_output_gff   = "";
    my $moby_exceptions  = [];

    # relleno los parametros por defecto Clover_call

    my $sequences      = $args{sequences}  || undef;
    my $matrices_input = $args{matrices}   || undef;
    my $format         = $args{format}     || "";
    my $parameters     = $args{parameters} || undef;
    my $debug          = $args{debug}      || 0;
    my $queryID        = $args{queryID}    || "";

    # Get the parameters

    my $pvalue_threshold     = $parameters->{pvalue_threshold};
    my $score_threshold      = $parameters->{score_threshold};
    my $database_parameter   = $parameters->{motif_database};
    my $matrix_mode          = $parameters->{matrix_mode};
    my $background_sequences = $parameters->{background_sequences};
    
    if ($debug) {
	print STDERR "pvalue threshold, $pvalue_threshold\n";
	print STDERR "score threshold, $score_threshold\n";
	print STDERR "motif database parameter, $database_parameter\n";
	print STDERR "matrix mode, $matrix_mode\n";
	print STDERR "background sequences, $background_sequences\n";
    }
    
    # Llama a Clover en local
    my $_clover_dir  = "/home/ug/gmaster/projects/Clover";
    my $_clover_bin  = "bin/clover";
    my $_clover_args = "-t $pvalue_threshold -s $score_threshold";
    
    my $_clover_bg_dir  = $_clover_dir . "/backgrounds";
    my $clover_bg_file = "";
    
    # Check that the binary is in place
    if (! -f "$_clover_dir/$_clover_bin") {
	my $note = "Internal System Error. Clover binary not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    # Clover_to_GFF perl script
    
    my $_clover_to_gff_bin = "/home/ug/gmaster/projects/clover_to_GFF/clover_to_GFF.pl";
    
    if (! -f "$_clover_to_gff_bin") {
	my $note = "Internal System Error. Clover Output to GFF conversion script not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    my ($matrix_fh, $matrix_file);

    if (defined $database_parameter) {

	if ($debug) {
	    print STDERR "matrix as a parameter...\n";
	}

	if ($matrix_mode eq "raw format") {

	    if ($debug) {
		print STDERR "raw mode\n";
	    }

	  SWITCH: {
	      if (lc ($database_parameter) eq "transfac") { $matrix_file = "$_clover_dir/matrices/Transfac-latest_raw_format.matrices.fa"; last SWITCH; }
	      if (lc ($database_parameter) eq "jaspar")   { $matrix_file = "$_clover_dir/matrices/Jaspar-latest_raw_format.matrices.fa"; last SWITCH; }
	      # Default is Transfac
	      $matrix_file = "$_clover_dir/matrices/Transfac_raw_format.matrices.fa";
	  }
	}
	elsif ($matrix_mode eq "log-likelihood") {
	    if ($debug) {
		print STDERR "log-likelihood mode\n";
	    }

	  SWITCH: {
	      if (lc ($database_parameter) eq "transfac") { $matrix_file = "$_clover_dir/matrices/Transfac-latest_likelihood.matrices.fa"; last SWITCH; }
	      if (lc ($database_parameter) eq "jaspar")   { $matrix_file = "$_clover_dir/matrices/Jaspar-latest_likelihood.matrices.fa"; last SWITCH; }
	      # Default is Transfac
	      $matrix_file = "$_clover_dir/matrices/Transfac-latest_likelihood.matrices.fa";
	  }
	}
	else {
	    # should be validated before...
	    print STDERR "don't know anything about matrix mode, $matrix_mode!\n";
	    exit 0;
	}
    }
    elsif (defined $matrices_input) {
	if ($debug) {
	    print STDERR "matrices as an input...\n";
	}

	# Make a temporary file with the matrix input
	$_clover_args .= " -s";
	eval {
	    ($matrix_fh, $matrix_file) = tempfile("/tmp/CLOVER_MATRIX.XXXXXX", UNLINK => 0);
	    print $matrix_fh "$matrices_input";
	    close $matrix_fh;
	};
	if ($@) {
	    my $note = "Internal System Error. Can not open Clover matrix input temporary file!\n";
	    my $code = 701;
	    print STDERR "$note\n";
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    return (undef, [$moby_exception]);
	}

    }
    else {
	# To change !
	print STDERR "matrices_input neither database_parameter are defined!!\n";
	exit 0;
    }

    if ((not defined $matrix_file) || (-z $matrix_file)) {
	# could well be possible if the input matrix set is empty, ie if MEME didn't predict any !!
	# But i guess in that case, no need to go up there, validation will be done before Clover_call
	print STDERR "Internal System Error. No defined matrix file!\n";
	return (undef, []);
    }
    
    # Background sequences parameter
    
    if (lc ($background_sequences) eq "human") {
	$clover_bg_file = $_clover_bg_dir . "/hs_chr20.mfa";
    }
    
    # Generate a temporary file locally with the sequence(s) in FASTA format
    # locally, ie not on a NFS mounted directory, for speed sake
    
    my ($seq_fh, $seqfile);
    eval {
	($seq_fh, $seqfile) = tempfile("/tmp/CLOVER_SEQS.XXXXXX", UNLINK => 0);
    };
    if ($@) {
	my $note = "Internal System Error. Can not open Clover sequence input temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }

    # Bioperl sequence factory

    my $sout = Bio::SeqIO->new (
				-fh     => $seq_fh,
				-format => 'fasta'
				);

    my @seqIds = keys (%$sequences);
    foreach my $sequenceIdentifier (@seqIds) {
	my $nucleotides = $sequences->{$sequenceIdentifier};

	# bioperl sequence object

	my $seqobj = Bio::Seq->new (
				    -display_id => $sequenceIdentifier,
				    -seq        => $nucleotides
				    );
	$sout->write_seq ($seqobj);
    }
    close $seq_fh;
    
    # Test empty file
    if (-z $seqfile) {
	my $note = "Internal System Error. Empty Clover input sequence file...\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    my ($clover_output_fh, $clover_output_filename) = tempfile("/tmp/CLOVER_OUTPUT.XXXXXX", UNLINK => 0);;
    close $clover_output_fh;
    
    if ($debug) {
	print STDERR "Running Clover, with this command:\n";
	print STDERR "$_clover_dir\/$_clover_bin $_clover_args $matrix_file $seqfile $clover_bg_file > $clover_output_filename\n";
    }
    
    # qx/$_clover_dir\/$_clover_bin $_clover_args $matrix_file $seqfile $clover_bg_file > $clover_output_filename/;
    
    unlink $seqfile unless $debug;
    
    ##
    # Test in taverna
    $clover_output_filename = "/tmp/CLOVER_OUTPUT.TEST";
    ##
    
    if (! -z $clover_output_filename) {
	
	# Convert it into GFF format
	$clover_output_gff = qx/$_clover_to_gff_bin $clover_output_filename/;
	
	if ($debug) {
	    # print STDERR "clover results in GFF format,\n$clover_output_gff\n";
	}
	
	unlink $clover_output_filename unless $debug;
	
	if ((! defined $clover_output_gff) || ($clover_output_gff eq "")) {
	    my $note = "Internal System Error. Converting Clover output into GFF format has failed!\n";
	    print STDERR "$note\n";
	    my $code = 701;
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    push (@$moby_exceptions, $moby_exception);
	    return (undef, $moby_exceptions);
	}
	
	return ($clover_output_gff, $moby_exceptions);
    }
    else {
	my $note = "Internal System Error. Clover has failed!\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	return (undef, $moby_exceptions);
    }
}


sub MetaAlignment_call {
    my %args = @_;

    # method output
    my $meta_output     = "";
    my $moby_exceptions = [];

    # relleno los parametros por defecto MetaAlignment_call

    my $map1 = $args{map1};
    my $map2 = $args{map2};
    my $queryID    = $args{queryID};
    my $parameters = $args{parameters} || undef;

    # Get the parameters

    my $alpha_penalty  = $parameters->{alpha_penalty};
    my $lambda_penalty = $parameters->{lambda_penalty};
    my $mu_penalty     = $parameters->{mu_penalty};
    my $debug = 0;
    
    my $output_format  = $parameters->{output_format};

    # Llama a Meta-alignment en local
    my $_meta_alignment_dir  = "/home/ug/gmaster/projects/Meta";
    my $_meta_alignment_bin  = "bin/meta";
    my $_meta_alignment_args = "-a $alpha_penalty -l $lambda_penalty -m $mu_penalty";

    # Check that the binary is in place
    if (! -f "$_meta_alignment_dir/$_meta_alignment_bin") {
	my $note = "Internal System Error. Meta-alignment binary not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }

    if ($output_format eq "GFF") {
	$_meta_alignment_args .= " -g";
    }

    # Create the temp map files

    my ($map1_fh, $map1_file);
    eval {
	($map1_fh, $map1_file) = tempfile("/tmp/META_MAP1.XXXXXX", UNLINK => 0);
	close ($map1_fh);
	open (FILE, ">$map1_file");
    };
    if ($@) {
	my $note = "Internal System Error. Can not open meta-alignment input (map1) temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    print FILE $map1;
    close FILE;

    # Create the temp map files

    my ($map2_fh, $map2_file);
    eval {
	($map2_fh, $map2_file) = tempfile("/tmp/META_MAP2.XXXXXX", UNLINK => 0);
	close ($map2_fh);
	open (FILE, ">$map2_file");
    };
    if ($@) {
	my $note = "Internal System Error. Can not open meta-alignment input (map2) temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    print FILE $map2;
    close FILE;
    
    # Sorting
    
    my $map1_sorted = qx/cat $map1_file | sort +3n/;
    
    open (FILE, ">$map1_file") or die "can't open temp file, $map1_file!\n";
    print FILE $map1_sorted;
    close FILE;

    my $map2_sorted = qx/cat $map2_file | sort +3n/;

    open (FILE, ">$map2_file") or die "can't open temp file, $map2_file!\n";
    print FILE $map2_sorted;
    close FILE;

    # Run meta

    # print STDERR "Running Meta-alignment, with this command:\n";
    # print STDERR "$_meta_alignment_dir\/$_meta_alignment_bin $_meta_alignment_args $map1_file $map2_file\n";

    my ($stdout_fh, $stdout_file);
    eval {
	($stdout_fh, $stdout_file) = tempfile("/tmp/META_OUTPUT.XXXXXX", UNLINK => 0);
	close $stdout_fh;
    };
    if ($@) {
	my $note = "Internal System Error. Can not open a temporary file to store meta-alignment output!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }

    if ($debug) {
	print STDERR "Executing meta-alignment,\n";
	print STDERR "$_meta_alignment_dir\/$_meta_alignment_bin $_meta_alignment_args $map1_file $map2_file > $stdout_file\n";
    }
    
    my @args = ("$_meta_alignment_dir\/$_meta_alignment_bin $_meta_alignment_args $map1_file $map2_file > $stdout_file");
    
    my $failed = system (@args);
    if ($failed > 0) {
	if (($! != ENOTTY) || ($! ne "Inappropriate ioctl for device")) {
	    # This is not an error, just mean that stdout is a terminal !!
	    print STDERR "Error, '$!'\n";
	}
	else {
	    my $note = "Internal System Error.  Meta-alignment system call died (with error code, $?).\n";
	    print STDERR "$note\n";
	    my $code = 701;
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    push (@$moby_exceptions, $moby_exception);
	    return (undef, $moby_exceptions);
	}
    }
    else {
	$meta_output = qx/cat $stdout_file/;
    }

    if (!$debug) {
	unlink $stdout_file;
	unlink $map1_file;
	unlink $map2_file;
    }
    
    if (defined $meta_output) {
	return ($meta_output, $moby_exceptions);
    }
    else {
	my $note = "Internal System Error. Meta-alignment has failed!\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	return (undef, $moby_exceptions);
    }
}


sub MultiMetaAlignment_call {
    my %args = @_;
    
    # method output
    my $mmeta_output     = "";
    my $moby_exceptions = [];
    
    # relleno los parametros por defecto MultiMetaAlignment_call
    
    my $maps       = $args{maps};
    my $queryID    = $args{queryID};
    my $parameters = $args{parameters} || undef;
    my $debug      = $args{debug} || 0;
    
    # Get the parameters
    
    my $alpha_penalty  = $parameters->{alpha_penalty};
    my $lambda_penalty = $parameters->{lambda_penalty};
    my $mu_penalty     = $parameters->{mu_penalty};
    my $gap_penalty    = $parameters->{gap_penalty};
    my $non_colinear_penalty = $parameters->{non_colinear_penalty};
    
    my $output_format  = $parameters->{output_format};

    # Llama a Multiple-Meta-alignment en local
    my $_mmeta_alignment_dir  = "/home/ug/gmaster/projects/MultiMeta";
    my $_mmeta_alignment_bin  = "bin/mmeta";
    my $_mmeta_alignment_args = "-a $alpha_penalty -l $lambda_penalty -m $mu_penalty -p $gap_penalty -c $non_colinear_penalty";
    
    # Check that the binary is in place
    if (! -f "$_mmeta_alignment_dir/$_mmeta_alignment_bin") {
	my $note = "Internal System Error. Multiple-Meta-alignment binary not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    if ($output_format eq "GFF") {
	$_mmeta_alignment_args .= " -g";
    }
    
    # Create the temp map files
    
    my ($maps_fh, $maps_file);
    eval {
	($maps_fh, $maps_file) = tempfile("/tmp/MMETA_MAPS.XXXXXX", UNLINK => 0);
	close ($maps_fh);
    };
    if ($@) {
	my $note = "Internal System Error. Can not open mmeta-alignment input (maps) temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    foreach my $map (@$maps) {
	
	if ($debug) {
	    print STDERR "concatenating map\n";
	}
	
	my $result = qx/echo "$map" | sort +3n >> $maps_file/;

	if ($debug) {
	    print STDERR "result, $result\n";
	}
	
    }
    
    # Run mmeta

    if ($debug) {
	print STDERR "Running Multiple-Meta-alignment, with this command:\n";
	print STDERR "$_mmeta_alignment_dir\/$_mmeta_alignment_bin $_mmeta_alignment_args $maps_file\n";
    }
    
    my ($stdout_fh, $stdout_file);
    eval {
	($stdout_fh, $stdout_file) = tempfile("/tmp/MMETA_OUTPUT.XXXXXX", UNLINK => 0);
	close $stdout_fh;
    };
    if ($@) {
	my $note = "Internal System Error. Can not open a temporary file to store mmeta-alignment output!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }

    if ($debug) {
	print STDERR "Executing mmeta-alignment,\n";
	print STDERR "$_mmeta_alignment_dir\/$_mmeta_alignment_bin $_mmeta_alignment_args $maps_file > $stdout_file\n";
    }
    
    my @args = ("$_mmeta_alignment_dir\/$_mmeta_alignment_bin $_mmeta_alignment_args $maps_file > $stdout_file");
    
    my $failed = system (@args);
    if ($failed > 0) {
	if (($! != ENOTTY) || ($! ne "Inappropriate ioctl for device")) {
	    # This is not an error, just mean that stdout is a terminal !!
	    print STDERR "Error, '$!'\n";
	}
	else {
	    my $note = "Internal System Error.  Multiple-Meta-alignment system call died (with error code, $?).\n";
	    print STDERR "$note\n";
	    my $code = 701;
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    push (@$moby_exceptions, $moby_exception);
	    return (undef, $moby_exceptions);
	}
    }
    else {
	$mmeta_output = qx/cat $stdout_file/;
    }

    if (!$debug) {
	unlink $stdout_file;
	unlink $maps_file;
    }
    
    if (defined $mmeta_output) {
	return ($mmeta_output, $moby_exceptions);
    }
    else {
	my $note = "Internal System Error. Multiple-Meta-alignment has failed!\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	return (undef, $moby_exceptions);
    }
}


sub generateScoreMatrix_call {
  my %args = @_;

  # output specs declaration
  my $matrix_output   = "";
  my $moby_exceptions = [];
  
  # relleno los parametros por defecto generateScoreMatrix_call
  
  my $inputdata_arrayref = $args{similarity_results};
  my $parameters         = $args{parameters} || undef;
  my $queryID            = $args{queryID}    || "";
  my $debug              = $args{debug}       || 0;
  
  # Llama a generateScoreMatrix script en local
  my $_application_dir  = "/home/ug/gmaster/projects/generateScoreMatrices";
  my $_application_bin  = "generateScoreMatrix.pl";
  my $_application_args    = "";
  
  # Check that the binary is in place
  if (! -f "$_application_dir/$_application_bin") {
      my $note = "Internal System Error. generateScoreMatrix script not found";
      print STDERR "$note\n";
      my $code = 701;
      my $moby_exception = INB::Exceptions::MobyException->new (
								code       => $code,
								type       => 'error',
								queryID    => $queryID,
								message    => "$note",
								);
      return (undef, [$moby_exception]);
  }
  
  # Generate a temporary file to store meta-alignment data
  
  my ($meta_fh, $meta_file);
  eval {
      ($meta_fh, $meta_file) = tempfile("/tmp/META_OUTPUT.XXXXXX", UNLINK => 0);
      print $meta_fh "@$inputdata_arrayref\n";
      close $meta_fh;
  };
  if ($@) {
      my $note = "Internal System Error. Can not open generateScoreMatrix meta input temporary file!\n";
      my $code = 701;
      print STDERR "$note\n";
      my $moby_exception = INB::Exceptions::MobyException->new (
								code       => $code,
								type       => 'error',
								queryID    => $queryID,
								message    => "$note",
								);
      return (undef, [$moby_exception]);
  }
  
  # Check the file is empty ?
  # I don't think so, all validation of the input data should be made before calling Factory.pm method
  
  # Run generateScoreMatrix

  if ($debug) {
      print STDERR "got " . @$inputdata_arrayref . " array references (meta alignments)\n";
      print STDERR "Running generateScoreMatrix.pl, with this command:\n";
      print STDERR "cat $meta_file | $_application_dir\/$_application_bin $_application_args\n";
  }
  
  $matrix_output = qx/cat $meta_file | $_application_dir\/$_application_bin $_application_args/;
  
  if ($debug) {
      print STDERR "matrix ouput, $matrix_output\n";
  }
  
  if (!$debug) {
      unlink $meta_file;
  }
  
  if (defined $matrix_output) {
      return ($matrix_output, $moby_exceptions);
  }
  else {
      print STDERR "Internal System Error. Matrix_output not defined!!\n";
      return (undef, $moby_exceptions);
  }

}

sub MEME_call {
    my %args = @_;

    # output specs declaration
    my $meme_output   = "";
    my $moby_exceptions = [];

    my $sequences  = $args{sequences};
    my $format     = $args{format};
    my $parameters = $args{parameters};
    my $_debug     = $args{debug};
    my $queryID    = $args{queryID}    || "";
    
    my $motif_distribution    = $parameters->{motif_distribution};
    my $maximum_number_motifs = $parameters->{maximum_number_motifs};
    my $minimum_number_sites  = $parameters->{minimum_number_sites};
    my $maximum_number_sites  = $parameters->{maximum_number_sites};
    my $minimum_motif_width   = $parameters->{minimum_motif_width};
    my $maximum_motif_width   = $parameters->{maximum_motif_width};
    my $e_value_cutoff        = $parameters->{e_value_cutoff};
    my $background_order      = $parameters->{background_order};

    # Llama a Meme en local
    my $_meme_dir   = "/usr/local/molbio/Install/meme-3.5.2";
    my $_meme_bin  = "bin/meme";
    my $_meme_args = "-nostatus -time 160 -maxiter 20 ";
    
    # Check that the binary is in place
    if (! -f "$_meme_dir/$_meme_bin") {
	my $note = "Internal System Error. MEME binary not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    # Setting up the MEME parameters

    # Default is HTML

    if ($_debug) {
	print STDERR "format, $format\n";
    }

    if ($format =~ /meme.text/i) {

	if ($_debug) {
	    print STDERR "MEME report will be text formatted\n";
	}

	$_meme_args .= "-text ";
    }

    if (not defined $motif_distribution) {
	$_meme_args .= "-mod zoops";
    }
    else {
      SWITCH: {
	  if ($motif_distribution eq "zero or one")               { $_meme_args .= "-mod zoops"; last SWITCH; }
	  if ($motif_distribution eq "one")                       { $_meme_args .= "-mod oops"; last SWITCH; }
	  if ($motif_distribution eq "any number of repetitions") { $_meme_args .= "-mod anr"; last SWITCH; }
	  # Default is "zero or one motif per sequence"
	  $_meme_args .= "-mod zoops";
      }
    }

    if (defined $maximum_number_motifs) {
	$_meme_args .= " -nmotifs $maximum_number_motifs";
    }
    else {
	$_meme_args .= " -nmotifs 3";
    }

    if (defined $minimum_number_sites) {
	$_meme_args .= " -minsites $minimum_number_sites";
    }
    else {
	# No default, ie let meme figure out what is the default !!
	# Will make it up on the fly, according a hardcoded protocol of their own.
    }

    if (defined $maximum_number_sites) {
	$_meme_args .= " -maxsites $maximum_number_sites";
    }
    else {
	# No default, ie let meme figure out what is the default !!
	# Will make it up on the fly, according a hardcoded protocol of their own.
    }

    if (defined $minimum_motif_width) {
	$_meme_args .= " -minw $minimum_motif_width";
    }
    else {
	$_meme_args .= " -minw 6";
    }

    if (defined $maximum_motif_width) {
	$_meme_args .= " -maxw $maximum_motif_width";
    }
    else {
	$_meme_args .= " -maxw 50";
    }

    if (defined $e_value_cutoff) {
	$_meme_args .= " -evt $e_value_cutoff";
    }
    else {
	$_meme_args .= " -evt 1";
    }

    # Generate a temporary file locally with the sequence(s) in FASTA format
    # locally, ie not on a NFS mounted directory, for speed sake
    
    my ($seq_fh, $seqfile);
    eval {
	($seq_fh, $seqfile) = tempfile("/tmp/MEME.XXXXXX", UNLINK => 0);
    };
    if ($@) {
	my $note = "Internal System Error. Can not open MEME input sequences temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    # Get the alphabet
    my $alphabet;
    
    # Bioperl sequence factory
    
    my $sout = Bio::SeqIO->new (
				-fh     => $seq_fh,
				-format => 'fasta'
				);
    
    my @seqIds = keys (%$sequences);
    foreach my $sequenceIdentifier (@seqIds) {
	my $nucleotides = $sequences->{$sequenceIdentifier};
	
	# bioperl sequence object
	
	my $seqobj = Bio::Seq->new (
				    -display_id => $sequenceIdentifier,
				    -seq        => $nucleotides
				    );
	
	if (not defined $alphabet) {
	    $alphabet = $seqobj->alphabet();
	    if ($_debug) {
		print STDERR "alphabet: $alphabet\n";
	    }
	}
	
	$sout->write_seq ($seqobj);
    }
    close $seq_fh;

    # Test empty file
    if (-z $seqfile) {
	print STDERR "Internal System Error when calling MEME software. Empty input sequences file...\n";
    }

    # MEME default is protein alphabet so if it is dna sequences, specify it, as the MEME model will be then different

    if (defined $alphabet && (($alphabet eq "dna") || ($alphabet eq "rna"))) {
	$_meme_args .= " -dna -revcomp";
    }

    my ($bfile, $bfh);
    # Train a background model if required
    if (defined $background_order && (lc ($background_order) eq "none")) {
	if ($_debug) {
	    print STDERR "no training...\n";
	}
    }
    else {
	if ($_debug) {
	    print STDERR "background training with order, $background_order\n";
	}

	# Do a background training
	my $training_bin = "bin/fasta-get-markov";

	my $training_args = "-m $background_order";
	if ($alphabet eq "protein") {
	    $training_args .= " -p";
	}

	eval {
	    ($bfh, $bfile) = tempfile("/tmp/MEME_BG.XXXXXX", UNLINK => 0);
	    close $bfh;
	    my @args = ("$_meme_dir/$training_bin $training_args < $seqfile > $bfile");
	    system (@args);
	};
	if ($@) {
	    my $note = "Internal System Error. Can not open MEME background temporary file!\n";
	    my $code = 701;
	    print STDERR "$note\n";
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    return (undef, [$moby_exception]);
	}
	
	if (-z $bfile) {
	    print STDERR "made an empty background model file while preparing meme execution!!\n";
	}
	
	$_meme_args .= " -bfile $bfile";
    }

    if ($_debug) {
	print STDERR "Running Meme, with this command:\n";
	print STDERR "$_meme_dir\/$_meme_bin $seqfile $_meme_args\n";
    }

    $meme_output = qx/source $_meme_dir\/etc\/meme.sh; $_meme_dir\/$_meme_bin $seqfile $_meme_args/;

    # Comment this line if you want to keep the file...
    unlink $seqfile unless $_debug;
    if (defined $bfile && (not $_debug)) {
	unlink $bfile;
    }

    if (defined $meme_output) {
	return ($meme_output, $moby_exceptions);
    }
    else {
	my $note = "Internal System Error. MEME has failed!\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	return (undef, $moby_exceptions);
    }
}


sub meme2matrix_call {
    my %args = @_;

    # output specs declaration
    my @PWMs = ();
    my $moby_exceptions = [];

    # relleno los parametros por defecto meme2matrix_call

    my $meme_predictions   = $args{meme_predictions} || undef;
    my $parameters         = $args{parameters}       || undef;
    my $_debug             = $args{debug};
    my $queryID            = $args{queryID}          || "";

    # Get the parameters

    my $matrix_mode = $parameters->{matrix_mode};

    # Llama a Meme2matrix en local
    my $_meme2matrix_dir  = "/home/ug/gmaster/projects/meme2matrix";
    my $_meme2matrix_bin  = "meme2matrix.pl";
    my $_meme2matrix_args = "";
    
    # Check that the binary is in place
    if (! -f "$_meme2matrix_dir/$_meme2matrix_bin") {
	my $note = "Internal System Error. meme2matrix script not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    if ($matrix_mode eq "raw format") {
	$_meme2matrix_args = "score";
    }
    elsif ($matrix_mode eq "log-likelihood") {
	$_meme2matrix_args = "probability";
    }
    
    # Create the temp map files

    my ($meme2matrix_fh, $meme2matrix_file);
    eval {
	($meme2matrix_fh, $meme2matrix_file) = tempfile("/tmp/MEME2MATRIX.XXXXXX", UNLINK => 0);
	print $meme2matrix_fh $meme_predictions;
	close $meme2matrix_fh;
    };
    if ($@) {
	my $note = "Internal System Error. Can not open meme2matrix input temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    if ($_debug) {
	print STDERR "Running meme2matrix, with this command:\n";
	print STDERR "$_meme2matrix_dir\/$_meme2matrix_bin $meme2matrix_file $_meme2matrix_args\n";
    }

    my $matrices_output = qx/$_meme2matrix_dir\/$_meme2matrix_bin $meme2matrix_file $_meme2matrix_args/;

    # Comment this line if you want to keep the file...
    unlink $meme2matrix_file;
    
   # Return an array of matrices from the output (which is a string)
    # Because we want to return a collection of text-formatted objects (each object containing a matrix)
    # So better to return an array rather than just a string

    if (! defined $matrices_output) {
	my $note = "Internal System Error. The parsing of MEME data has failed!\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	return (undef, $moby_exceptions);
    }
    
    # To split the matrices, check which position they start, ie chech th position of 'MEME'
    my @positions = ();
    while ($matrices_output =~ m/MEME/g) {
	push (@positions, pos ($matrices_output) - 4);
    }
    
    if ($_debug) {
	print STDERR "positions: @positions\n";
    }
    
    if (@positions == 0) {
	# It's not an error, it is a warning !!
	# it can be normal not find find motifs !!
	my $note = "meme2matrix didn't find any motif matrix in meme data file!\n";
	print STDERR "$note\n";
	my $code = 700;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'information',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	return ([], $moby_exceptions);
    }
    
    my $position_1 = $positions[0];
    my $index = 1;
    while (my $position_2 = $positions[$index]) {
	my $matrix = substr ($matrices_output, $position_1, $position_2 - $position_1);
	push (@PWMs, $matrix);
	$position_1 = $position_2;
	$index++;
    }
    my $end = length ($matrices_output);
    my $matrix = substr ($matrices_output, $position_1, $end - $position_1);
    push (@PWMs, $matrix);

    if ($_debug) {
	print STDERR "nb matrices: " . @PWMs . ".\n";
	print STDERR "Dumping matrices array,\n" . Dumper (@PWMs) . "\n";
    }

    if (@PWMs > 0) {
	return (\@PWMs, $moby_exceptions);
    }
    else {
	# It's not an error, it is a warning !!
	# it can be normal not find find motifs !!
	my $note = "meme2matrix didn't find any motif matrix in meme data file!\n";
	print STDERR "$note\n";
	my $code = 700;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'information',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	return ([], $moby_exceptions);
    }

}

sub RepeatMasker_call {
    my %args = @_;
    
    # output specs declaration
    my @masked_seqs_fasta = "";
    my $moby_exceptions   = [];
    
    # relleno los parametros por defecto RepeatMasker_call
    
    my $sequences          = $args{sequences}  || undef;
    my $parameters         = $args{parameters} || undef;
    my $debug              = $args{debug};
    my $queryID            = $args{queryID}    || "";
    
    # parameters

    my $species = $parameters->{species};
    my $engine  = $parameters->{engine};
    
    # Llama a RepeatMasker en local
    my $_repeatmasker_dir  = "/usr/local/molbio/Install/repeatmasker-3.1.5";
    my $_repeatmasker_bin  = "RepeatMasker";
    my $_repeatmasker_args = "";
    
    if ((defined $species) && (lc ($species) ne "none")) {
	$_repeatmasker_args .= "-species $species ";
    }

    if (defined $engine) {
	$_repeatmasker_args .= "-e $engine";
    }
    else {
	# Default is cross-match
	$_repeatmasker_args .= "-e crossmatch";
    }
    
    # Check that the binary is in place
    if (! -f "$_repeatmasker_dir/$_repeatmasker_bin") {
	my $note = "Internal System Error. repeatmasker script not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    # Create the temp sequences files

    my ($seq_fh, $repeatmasker_file);
    eval {
	($seq_fh, $repeatmasker_file) = tempfile("/tmp/REPEATMASKER.XXXXXX", UNLINK => 0);
    };
    if ($@) {
	my $note = "Internal System Error. Can not open repeatmasker input temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    # Bioperl sequence factory
    my $sout = Bio::SeqIO->new (
				-fh     => $seq_fh,
				-format => 'fasta'
				);
    
    my @seqIds = keys (%$sequences);
    foreach my $sequenceIdentifier (@seqIds) {
	my $nucleotides = $sequences->{$sequenceIdentifier};
	
	# bioperl sequence object
	my $seqobj = Bio::Seq->new (
				    -display_id => $sequenceIdentifier,
				    -seq        => $nucleotides
				    );
	$sout->write_seq ($seqobj);
    }
    close $seq_fh;
    
    # Test empty file
    if (-z $repeatmasker_file) {
	my $note = "Internal System Error. Empty repeatmasker input sequence file...\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    if ($debug) {
	print STDERR "Running repeatmasker, with this command:\n";
	print STDERR "$_repeatmasker_dir\/$_repeatmasker_bin $_repeatmasker_args $repeatmasker_file\n";
    }
    
    my $result = qx/$_repeatmasker_dir\/$_repeatmasker_bin $_repeatmasker_args $repeatmasker_file >& \/dev\/null/;
    
    my $masked_sequences_file = $repeatmasker_file . ".masked";
    my $masked_sequences      = qx/cat $masked_sequences_file/;
    
    # Comment this line if you want to keep the file...
    if (! $debug) {
	unlink $repeatmasker_file;
	unlink $repeatmasker_file . ".log";
	unlink $repeatmasker_file . ".masked";
	unlink $repeatmasker_file . ".cat";
	unlink $repeatmasker_file . ".cat.all";
	unlink $repeatmasker_file . ".tbl";
	unlink $repeatmasker_file . ".out";
    }
    
    if (! defined $masked_sequences) {
	my $note = "Internal System Error. The parsing of the masked sequence data has failed!\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	return (undef, $moby_exceptions);
    }
    else {
	return ($masked_sequences, $moby_exceptions);
    }
    
    
}

sub Dust_call {
    my %args = @_;
    
    # output specs declaration
    my @masked_seqs_fasta = "";
    my $moby_exceptions   = [];
    
    # relleno los parametros por defecto Duct_call
    
    my $sequences          = $args{sequences}  || undef;
    my $parameters         = $args{parameters} || undef;
    my $debug              = $args{debug};
    my $queryID            = $args{queryID}    || "";
    
    # No parameters
    
    # Llama a Dust en local
    my $_dust_dir  = "/home/ug/gmaster/projects/bin";
    my $_dust_bin  = "dust";
    
    # Check that the binary is in place
    if (! -f "$_dust_dir/$_dust_bin") {
	my $note = "Internal System Error. dust binary not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    # Create the temp sequences file

    my ($seq_fh, $dust_file);
    eval {
	($seq_fh, $dust_file) = tempfile("/tmp/DUST.XXXXXX", UNLINK => 0);
    };
    if ($@) {
	my $note = "Internal System Error. Can not open dust input temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    # Bioperl sequence factory
    my $sout = Bio::SeqIO->new (
				-fh     => $seq_fh,
				-format => 'fasta'
				);
    
    my @seqIds = keys (%$sequences);
    foreach my $sequenceIdentifier (@seqIds) {
	my $nucleotides = $sequences->{$sequenceIdentifier};
	
	# bioperl sequence object
	my $seqobj = Bio::Seq->new (
				    -display_id => $sequenceIdentifier,
				    -seq        => $nucleotides
				    );
	$sout->write_seq ($seqobj);
    }
    close $seq_fh;

    # Test empty file
    if (-z $dust_file) {
	my $note = "Internal System Error. Empty dust input sequence file...\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    if ($debug) {
	print STDERR "Running dust, with this command:\n";
	print STDERR "$_dust_dir\/$_dust_bin $dust_file\n";
    }

    my $masked_sequences = qx/$_dust_dir\/$_dust_bin $dust_file/;
    
    # Comment this line if you want to keep the file...
    unlink $dust_file;
    
    if (! defined $masked_sequences) {
	my $note = "Internal System Error. The parsing of masked sequence data has failed!\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	return (undef, $moby_exceptions);
    }
    else {
	return ($masked_sequences, $moby_exceptions);
    }
    
}

sub Dust_FASTA_call {
    my %args = @_;
    
    # output specs declaration
    my @masked_seqs_fasta = "";
    my $moby_exceptions   = [];
    
    # relleno los parametros por defecto Duct_call
    
    my $sequences          = $args{sequences}  || undef;
    my $parameters         = $args{parameters} || undef;
    my $debug              = $args{debug};
    my $queryID            = $args{queryID}    || "";
    
    # No parameters
    
    # Llama a Dust en local
    my $_dust_dir  = "/home/ug/gmaster/projects/bin";
    my $_dust_bin  = "dust";
    
    # Check that the binary is in place
    if (! -f "$_dust_dir/$_dust_bin") {
	my $note = "Internal System Error. dust binary not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    # Create the temp sequences file

    my ($seq_fh, $dust_file);
    eval {
	($seq_fh, $dust_file) = tempfile("/tmp/DUST.XXXXXX", UNLINK => 0);
    };
    if ($@) {
	my $note = "Internal System Error. Can not open dust input temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    print $seq_fh $sequences;
    close $seq_fh;

    # Test empty file
    if (-z $dust_file) {
	my $note = "Internal System Error. Empty dust input sequence file...\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    if ($debug) {
	print STDERR "Running dust, with this command:\n";
	print STDERR "$_dust_dir\/$_dust_bin $dust_file\n";
    }
    
    my $masked_sequences = qx/$_dust_dir\/$_dust_bin $dust_file/;
    
    # Comment this line if you want to keep the file...
    unlink $dust_file;
    
    if (! defined $masked_sequences) {
	my $note = "Internal System Error. The parsing of masked sequence data has failed!\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	return (undef, $moby_exceptions);
    }
    else {
	return ($masked_sequences, $moby_exceptions);
    }
    
}

sub CrossMatchToScreenVector_call {
    my %args = @_;
    
    # output specs declaration
    my @screened_seqs_fasta = "";
    my $moby_exceptions     = [];
    
    # relleno los parametros por defecto CrossMatchToScreenVector_call
    
    my $fasta_sequences_str = $args{sequences}  || undef;
    my $parameters          = $args{parameters} || undef;
    my $debug               = $args{debug};
    my $queryID             = $args{queryID}    || "";
    
    # parameters
    
    my $minmatch = $parameters->{minmatch};
    my $minscore = $parameters->{minscore};
    
    # Llama a Cross_Match en local
    my $_cross_match_dir    = "/usr/local/molbio/Install/cross_match";
    my $_cross_match_bin    = "cross_match";
    my $_cross_match_args   = "-screen ";
    my $_cross_match_vectors = "/home/ug/gmaster/projects/assembly/data_libraries/UniVec_Core.fa";

    if (defined $minmatch) {
	$_cross_match_args .= "-minmatch $minmatch ";
    }
    if (defined $minscore) {
	$_cross_match_args .= "-minscore $minscore ";
    }
    
    # Check that the binary is in place
    if (! -f "$_cross_match_dir/$_cross_match_bin") {
	my $note = "Internal System Error. cross_match binary not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    # Create the temp sequences file

    my ($seq_fh, $cross_match_file);
    eval {
	($seq_fh, $cross_match_file) = tempfile("/tmp/CROSS_MATCH.XXXXXX", UNLINK => 0);
    };
    if ($@) {
	my $note = "Internal System Error. Can not open cross_match input temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    print $seq_fh $fasta_sequences_str;
    close $seq_fh;

    # Test empty file
    if (-z $cross_match_file) {
	my $note = "Internal System Error. Empty cross_match input sequence file...\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    if ($debug) {
	print STDERR "Running cross_match, with this command:\n";
	print STDERR "$_cross_match_dir\/$_cross_match_bin $cross_match_file $_cross_match_vectors $_cross_match_args\n";
    }

    my $result = qx/$_cross_match_dir\/$_cross_match_bin $cross_match_file $_cross_match_vectors $_cross_match_args >& \/dev\/null/;
    my $screened_sequences_file = $cross_match_file . ".screen";
    my $screened_sequences = qx/cat $screened_sequences_file/;
    
    # Comment this line if you want to keep the file...
    if (! $debug) {
	unlink $cross_match_file;
	unlink $screened_sequences_file;
	unlink $cross_match_file . ".log"
    }
    
    if (! defined $screened_sequences) {
	my $note = "Internal System Error. The parsing of the screened sequence data has failed!\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	return (undef, $moby_exceptions);
    }
    else {
	return ($screened_sequences, $moby_exceptions);
    }   
    
}

sub Phrap_call {
    my %args = @_;
    
    # output specs declaration
    my $assembled_sequences;
    my $ace_data;
    my $moby_exceptions      = [];
    
    # relleno los parametros por defecto Phrap_call
    
    my $sequences          = $args{sequences}    || undef;
    my $fasta_quality_data = $args{quality_data} || undef;
    my $parameters         = $args{parameters}   || undef;
    my $debug              = $args{debug};
    my $queryID            = $args{queryID}      || "";
    
    # parameters
    
    my $node_seg   = $parameters->{node_seg};
    my $node_space = $parameters->{node_space};
    
    # Llama a Phrap en local
    my $_phrap_dir    = "/usr/local/molbio/Install/cross_match";
    my $_phrap_bin    = "phrap";
    # $_phrap_bin    = "phrap.bigmalloc";
    # $_phrap_bin    = "phrap.manyreads";
    my $_phrap_args   = "-new_ace";
    
    my %ace_parsing_options;
    if (($_phrap_bin eq "phrap") || ($_phrap_bin eq "phrap.manyreads") || ($_phrap_bin eq "phrap.bigmalloc")) {
	%ace_parsing_options = (
				'ace_name_suffix' => '.ace',
				'contig_pattern'  => 'CO\s(\S+)\s+.*',
				'AF_test'         => '^AF',
				'AF_pattern'      => 'AF\s+(\S+)\s+(C|U)\s+-?\d+.*',
				'member_sequence' => 'RD\s+(\S+).*',
				);
    }
    else {
	# Old formatting
	%ace_parsing_options = (
				'ace_name_suffix'  => '.ace',
				'contig_pattern'   => 'DNA\s(\S+)',
				'AF_test'          => '^Assembled_from',
				'AF_pattern'       => 'Assembled_from\*?\s+(\S+)\s+-?\d+\s+-?\d+',
				'member_sequence'  => 'DNA\s(\S+)',
				);
    }
    
    if (defined $node_seg) {
	$_phrap_args .= " -node_seg $node_seg";
    }
    if (defined $node_space) {
	$_phrap_args .= " -node_space $node_space";
    }
    
    # Check that the binary is in place
    if (! -f "$_phrap_dir/$_phrap_bin") {
	my $note = "Internal System Error. phrap binary not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, undef, [$moby_exception]);
    }
    
    # Create the temp sequences file

    my ($seq_fh, $sequences_phrap_file);
    eval {
	($seq_fh, $sequences_phrap_file) = tempfile("/tmp/PHRAP.XXXXXX", UNLINK => 0);
    };
    if ($@) {
	my $note = "Internal System Error. Can not open phrap input temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, undef, [$moby_exception]);
    }
    
    print $seq_fh "$sequences";
    close $seq_fh;
    
    if (-z $sequences_phrap_file) {
	    my $note = "Internal System Error. Empty phrap input sequences data file...\n";
	    print STDERR "$note\n";
	    my $code = 701;
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    return (undef, undef, [$moby_exception]);
	}
    
    # Parse the quality data if any
    
    my $phrap_quality_file;
    if (defined $fasta_quality_data) {
	
	$phrap_quality_file = $sequences_phrap_file . ".qual";
	
	if ($debug) {
	    print STDERR "quality data hash defined!\n";
	}
	
	open FILE, ">$phrap_quality_file" or die "can't open temporary file for phrap quality data, $phrap_quality_file!\n";
	print FILE "$fasta_quality_data";
	close FILE;
	
	if (-z $phrap_quality_file) {
	    my $note = "Internal System Error. Empty phrap input quality data file...\n";
	    print STDERR "$note\n";
	    my $code = 701;
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    return (undef, undef, [$moby_exception]);
	}
	
    }
           
    # Free memory
    undef $sequences;
    undef $fasta_quality_data;
    
    if ($debug) {
	print STDERR "Running phrap, with this command:\n";
	print STDERR "$_phrap_dir\/$_phrap_bin $sequences_phrap_file $_phrap_args\n";
    }
    
    my $result;
    if ($debug) {
	$result = qx/$_phrap_dir\/$_phrap_bin $sequences_phrap_file $_phrap_args/;
    }
    else {
	$result = qx/$_phrap_dir\/$_phrap_bin $sequences_phrap_file $_phrap_args >& \/dev\/null/;
    }
    
    if ($debug) {
	print STDERR "Phrap done\n";
    }
    
    if (! -f "$sequences_phrap_file.singlets") {
	my $note = "Phrap failed";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, undef, [$moby_exception]);
    }
    
    # Process the phrap results to get a unique FASTA file
    
    my $singlets_seqs = qx/cat $sequences_phrap_file".singlets"/;
    
    # Replace contigs of one member by singlets
    # Also gives some information on the cluster members
    # To do so, parse the ace file
    
    my $ace_filename  = $sequences_phrap_file . ".ace";
    my ($clusters_list, $extra_singlets_list) = _parse_ace_file ($ace_filename, $debug, %ace_parsing_options);
    
    my $pattern = join ('|', @$extra_singlets_list);
    
    if ($debug) {
	print STDERR "pattern, $pattern\n";
    }
    
    my $sequence_search_bin = "/home/ug/gmaster/projects/bin/sequence_search.pl";
    
    # Extract the singlet sequences from the sequences input file
    my $extra_singlets_seqs = qx/$sequence_search_bin -f $sequences_phrap_file -h -p \"$pattern\"/;
    
    # reprocess the contigs, using clusters_list data
    
    my $contigs_seqs = "";
    
    foreach my $cluster_info_href (@$clusters_list) {
	my $contig_href = $cluster_info_href->{cluster_info};
	my $header = $contig_href->{line_info};
	my $seq    = $contig_href->{sequence};
	
	$contigs_seqs .= ">$header\n";
	# FASTA Formatting of the contig sequence
	while ($seq =~ /(.{1,60})/g) {
	    $contigs_seqs .= "$1\n";
	}
    }
    
    $assembled_sequences = $contigs_seqs . $singlets_seqs . $extra_singlets_seqs;
    
    $ace_data = qx/cat $ace_filename/;
    
    # Comment this line if you want to keep the file...
    if (!$debug) {
	# The phrap input files
	unlink $sequences_phrap_file;
	if (defined $phrap_quality_file) {
	    unlink $phrap_quality_file;
	}
	# Also the phrap output files
	unlink $sequences_phrap_file . ".log";
	unlink $sequences_phrap_file . ".ace";
	unlink $sequences_phrap_file . ".contigs";
	unlink $sequences_phrap_file . ".contigs.qual";
	unlink $sequences_phrap_file . ".singlets";
	unlink $sequences_phrap_file . ".problems";
	unlink $sequences_phrap_file . ".problems.qual";
    }
    
    if (! defined $assembled_sequences) {
	my $note = "Internal System Error. The parsing of the assembled sequence data has failed!\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	return (undef, undef, $moby_exceptions);
    }
    else {
	return ($assembled_sequences, $ace_data, $moby_exceptions);
    }
    
}

sub Phred_call {
    my %args = @_;
    
    # output specs declaration
    my @seqs_fasta      = ();
    my @quality_fasta   = ();
    my $moby_exceptions = [];
    
    # relleno los parametros por defecto Phred_call
    
    my $chromatograms = $args{chromatograms} || undef;
    my $parameters    = $args{parameters}    || undef;
    my $debug         = $args{debug};
    my $queryID       = $args{queryID}       || "";
    
    # parameters
    
    my $trim_alt    = $parameters->{trim_alt};
    my $trim_cutoff = $parameters->{trim_cutoff};
    
    # Llama a Phred en local
    my $_phred_dir     = "/usr/local/molbio/Install/phred-0.020425.c";
    my $_phred_bin     = "phred";
    my $_phred_args    = "";
    my $_phred_config_file = "/home/ug/gmaster/projects/assembly/phredpar.dat";
    my $_phd2fasta_bin = "phd2fasta";
    
    if ($trim_alt) {
	$_phred_args .= "-trim_alt \"\" -trim_phd ";
    }
    
    if (defined $trim_cutoff) {
      $_phred_args .= "-trim_cutoff $trim_cutoff ";
    }
    
    # Check that the binary is in place
    if (! -f "$_phred_dir/$_phred_bin") {
	my $note = "Internal System Error. phred binary not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, undef, [$moby_exception]);
    }
    
    # Prepare the input chromatograms data
    
    # Set up a temp input dir and a temp output dir
    my $phred_input_dir = tempdir( "/tmp/PHRED.INPUT.XXXXXX" );
    my $phred_output_dir = tempdir( "/tmp/PHRED.OUTPUT.XXXXXX" );
    
    if ($debug) {
	print STDERR "phred_input_dir, $phred_input_dir, phred_output_dir, $phred_output_dir\n";
    }
    
    my @chromatogram_ids = keys (%$chromatograms);
    foreach my $chromatogram_id (@chromatogram_ids) {
	
	print STDERR "chromatogram identifier, $chromatogram_id.\n";
	
	my $chromatogram = $chromatograms->{$chromatogram_id};
	
	my $chromatogram_type     = $chromatogram->{type};
	my $chromatogram_data_b64 = $chromatogram->{rawdata};
	
	my $tempfile = "$phred_input_dir/$chromatogram_id" . "." . $chromatogram_type;
	open FILE, ">$tempfile";
	
	if ($debug) {
	    print STDERR "temp input file name, $tempfile\n";
	}
	
	# decode
	my $chromatogram_data = decode_base64 ($chromatogram_data_b64);
	# put the data in it
	print FILE "$chromatogram_data";
	close FILE;
    }
    
    # Set up a temporary output sequence FASTA file and quality FASTA file
    
    my ($seq_fasta_fh, $seq_fasta_file) = tempfile ("/tmp/PHRED_SEQS_XXXXX", SUFFIX => ".dna");
    close $seq_fasta_fh;
    my $qual_fasta_file = $seq_fasta_file . ".qual";
    
    if ($debug) {
	print STDERR "fasta sequence file, $seq_fasta_file\n";
	print STDERR "fasta quality file, $qual_fasta_file\n";
    }
    
    if ($debug) {
	print STDERR "Running phred, with this command:\n";
	print STDERR "$_phred_dir\/$_phred_bin -id $phred_input_dir -pd $phred_output_dir $_phred_args\n";
	print STDERR "then running:\n";
	print STDERR "$_phred_dir\/$_phd2fasta_bin -id $phred_output_dir -os $seq_fasta_file -oq $qual_fasta_file\n";
    }
    
    my $result = qx/export PHRED_PARAMETER_FILE=$_phred_config_file; $_phred_dir\/$_phred_bin -id $phred_input_dir -pd $phred_output_dir $_phred_args/;
    # then
    $result = qx/$_phred_dir\/$_phd2fasta_bin -id $phred_output_dir -os $seq_fasta_file -oq $qual_fasta_file/;
    
    # Process the results
    
    my $fasta_sequences = qx/cat $seq_fasta_file/;
    my $fasta_quality   = qx/cat $qual_fasta_file/;
    
    # Remove temprary files
    
    if (! $debug) {
	unlink $seq_fasta_file;
	unlink $qual_fasta_file;
	# doesn't work !!!????
	# rmtree (["$phred_input_dir", "$phred_output_dir"], 1, 1);
    }
    
    if ((! defined $fasta_sequences) || (length $fasta_sequences < 1)) {
	my $note = "Internal System Error. Phred base calling processing has failed!\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	return (undef, undef, $moby_exceptions);
    }
    
    return ($fasta_sequences, $fasta_quality, $moby_exceptions);
}

sub SequenceFilteringByLength_call {
    my %args = @_;
    
    # output specs declaration
    my $filtered_sequences;
    my $filtered_quality_data = "";
    my $moby_exceptions    = [];
    
    # relleno los parametros por defecto de SequenceFilteringByLength_call
    
    my $sequences          = $args{sequences}    || undef;
    my $fasta_quality_data = $args{quality_data} || undef;
    my $parameters         = $args{parameters}   || undef;
    my $debug              = $args{debug};
    my $queryID            = $args{queryID}      || "";
    
    # parameters
    
    my $trim_masked_regions = $parameters->{trim_masked_regions};
    my $length_cutoff       = $parameters->{length_cutoff};
    
    # Llama a sequence filtering script en local
    my $_sequence_filtering_dir    = "/home/ug/gmaster/projects/bin";
    my $_sequence_filtering_bin    = "filter_sequences_by_length.pl";
    my $_sequence_filtering_args   = "-c $length_cutoff";
    
    if ($trim_masked_regions) {
	$_sequence_filtering_args .= " -t";
    }
    
    # Check that the binary is in place
    if (! -f "$_sequence_filtering_dir/$_sequence_filtering_bin") {
	my $note = "Internal System Error. $_sequence_filtering_bin script not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, undef, [$moby_exception]);
    }
    
    # Create the temp sequences file

    my ($seq_fh, $sequences_file);
    eval {
	($seq_fh, $sequences_file) = tempfile("/tmp/SEQ_FILTERING.XXXXXX", UNLINK => 0);
    };
    if ($@) {
	my $note = "Internal System Error. Can not open sequences input temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, undef, [$moby_exception]);
    }
    
    print $seq_fh "$sequences";
    close $seq_fh;
    
    if (-z $sequences_file) {
	my $note = "Internal System Error. Empty input sequences data file...\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								      );
	return (undef, undef, [$moby_exception]);
    }
    
    # Parse the quality data if any
    
    my $quality_data_file;
    if (defined $fasta_quality_data) {
	
	$quality_data_file = $sequences_file . ".qual";
	
	if ($debug) {
	    print STDERR "quality data hash defined!\n";
	}
	
	open FILE, ">$quality_data_file" or die "can't open temporary file for quality quality data, $quality_data_file!\n";
	print FILE "$fasta_quality_data";
	close FILE;
	
	if (-z $quality_data_file) {
	    my $note = "Internal System Error. Empty input quality data file...\n";
	    print STDERR "$note\n";
	    my $code = 701;
	    my $moby_exception = INB::Exceptions::MobyException->new (
								      code       => $code,
								      type       => 'error',
								      queryID    => $queryID,
								      message    => "$note",
								      );
	    return (undef, undef, [$moby_exception]);
	}
    }
           
    if ($debug) {
	print STDERR "Running sequence filtering script, with this command:\n";
	print STDERR "$_sequence_filtering_dir\/$_sequence_filtering_bin -i $sequences_file $_sequence_filtering_args\n";
    }
    
    $filtered_sequences = qx/$_sequence_filtering_dir\/$_sequence_filtering_bin -i $sequences_file $_sequence_filtering_args/;
    
    if (defined $quality_data_file) {
	
	# Put in sync the quality data
	
	# Get also a list of removed sequences
	
	my $removed_sequences_list = qx/$_sequence_filtering_dir\/$_sequence_filtering_bin -i $sequences_file $_sequence_filtering_args -l/;
	
	$removed_sequences_list =~ s/\n/|/g;
	chop $removed_sequences_list;
	
	if ($debug) {
	    print STDERR "removed_sequences_list pattern, $removed_sequences_list\n";
	}
	
	if (length $removed_sequences_list > 0) {
	    my @quality_data_array = split ("\n", $fasta_quality_data);

	    if ($debug) {
		print STDERR "Quality array length, " . @quality_data_array  ."\n";
	    }
	    
	    # Filtering processing
	    
	    my $keep = 0;
	    foreach my $line (@quality_data_array) {
		chomp $line;
		if ($line =~ /^>/) {
		    
		    if ($debug) {
			print STDERR "header line\n";
		    }
		    
		    if ($line =~ /$removed_sequences_list/) {

			if ($debug) {
			    print STDERR "header match!\n";
			}
			
			$keep = 0;
		    }
		    else {
			$keep = 1;
			$filtered_quality_data .= "$line\n";
		    }
		}
		elsif ($keep) {
		    $filtered_quality_data .= "$line\n";
		}
	    }
	}
	else {
	    # Empty, ie no sequences have been removed !
	    $filtered_quality_data = $fasta_quality_data;
	}
    }
    
    if (!$debug) {
	unlink $sequences_file;
	if (defined $quality_data_file) {
	    unlink $quality_data_file;
	}
    }
    
    if (! defined $filtered_sequences || (length $filtered_sequences < 1)) {
	my $note = "Internal System Error. the sequence fitlering script execution has failed!\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	push (@$moby_exceptions, $moby_exception);
	return (undef, undef, $moby_exceptions);
    }
    
    if ($debug) {
	print STDERR "Sequence Filtering done\n";
    }
    
    return ($filtered_sequences, $filtered_quality_data, $moby_exceptions);
    
}

sub KMeans_call {
    my %args = @_;
    
    # output specs declaration
    my $gene_clusters_aref = [];
    my $moby_exceptions    = [];
    
    # relleno los parametros por defecto de SequenceFilteringByLength_call
    
    my $gene_matrix        = $args{gene_matrix} || undef;
    my $parameters         = $args{parameters}  || undef;
    my $debug              = $args{debug}       || 0;
    my $queryID            = $args{queryID}     || "";
    
    # parameters
    
    my $gene_centering   = $parameters->{gene_centering};
    my $iteration_number = $parameters->{iteration_number};
    my $cluster_number   = $parameters->{cluster_number};
    
    # Llama a cluster binary en local
    my $_cluster_dir    = "/home/ug/gmaster/projects/cluster";
    my $_cluster_bin    = "cluster";
    my $_cluster_args   = "-k $cluster_number -r $iteration_number";
    
    # Add the gene centering option
    
    if ($gene_centering =~ /none/i) {
	# nothing - it is default !
    }
    elsif ($gene_centering =~ /k-means/i) {
	$_cluster_args .= " -cg a";
    }
    elsif ($gene_centering =~ /k-medians/i) {
	$_cluster_args .= " -cg m";
    }
    else {
	print STDERR "problem with the 'gene centering' parameter, $gene_centering, is unknown!\n";
    }
    
    # Check that the binary is in place
    if (! -f "$_cluster_dir/$_cluster_bin") {
	my $note = "Internal System Error. $_cluster_bin script not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    # Create the temp sequences file

    my ($gene_matrix_fh, $gene_matrix_file);
    eval {
	($gene_matrix_fh, $gene_matrix_file) = tempfile("/tmp/GENE_MATRIX.XXXXXX", UNLINK => 0);
    };
    if ($@) {
	my $note = "Internal System Error. Can not open gene matrix input temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    # Be careful, i add an extra tabulation because the xml parsing didn't keep in place !!!
    print $gene_matrix_fh "$gene_matrix\t\n";
    close $gene_matrix_fh;
    
    if (-z $gene_matrix_file) {
	my $note = "Internal System Error. Empty input gene matrix data file...\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								      );
	return (undef, [$moby_exception]);
    }
    
    # output prefix
    
    my $_output_prefix  = $gene_matrix_file;
    my $output_filename = $_output_prefix . "_K_G" . $cluster_number . ".kgg";
    $_cluster_args     .= " -u $_output_prefix";
    
    if ($debug) {
	print STDERR "Running k-means clustering, with this command:\n";
	print STDERR "$_cluster_dir\/$_cluster_bin -f $gene_matrix_file $_cluster_args\n";
    }
    
    my $result = qx/$_cluster_dir\/$_cluster_bin -f $gene_matrix_file $_cluster_args/;
    chomp $result;
    
    if ($debug) {
	print STDERR "clustering result, $result\n";
    }
    
    if (! -f $output_filename) {
	my $note = "Internal System Error. K-means clustering has failed, here the error that has been given back by cluster software, '$result'\n";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	
	return (undef, [$moby_exception]);
    }
    else {
	$result    = qx/cat $output_filename/;
	
	# set up the array of clusters
	
	if ($debug) {
	    print STDERR "parsing k-means clustering output file...\n";
	}
	
	my @lines = split ('\n', $result);
	my $cluster_index;
	my $cluster = "";
	foreach my $line (@lines) {
	    if ($line =~ /^([^\t]+)\t(\d+)/) {
		my $gene_identifier    = $1;
		my $cluster_index_tmp  = $2;
		
		if (defined $cluster_index_tmp) {

		    if ($debug) {
			print STDERR "gene identifier, $gene_identifier\n";
			print STDERR "cluster index, $cluster_index_tmp\n";
		    }
		    
		    if (! defined $cluster_index || $cluster_index eq $cluster_index_tmp) {
			# Same cluster
			$cluster .= "$gene_identifier\n";
			$cluster_index = $cluster_index_tmp;
		    }
		    else {
			# new cluster
			push (@$gene_clusters_aref, $cluster);
			
			# Initialisation
			$cluster = "";
			
			$cluster .= "$gene_identifier\n";
			$cluster_index = $cluster_index_tmp;
		    }
		}
		else {
		    print STDERR "problem parsing cluster info line, $line\n";
		}
	    }
	}
	
	if (defined $cluster_index) {
	    push (@$gene_clusters_aref, $cluster);
	}
	
	if ($debug) {
	    print STDERR "parsing done!\n";
	}
	
	if (! $debug) {
	    unlink $gene_matrix_file;
	    unlink $output_filename;
	    unlink $gene_matrix_file . "_K_G" . $cluster_number . ".cdt";
	}
	
	return ($gene_clusters_aref, $moby_exceptions);
    }
}

sub GFF2JPEG_call {
    my %args = @_;
    
    # output specs declaration
    my $image = undef;
    my $moby_exceptions    = [];
    
    my $gff_maps           = $args{maps}       || undef;
    my $parameters         = $args{parameters} || undef;
    my $debug              = $args{debug}      || 0;
    my $queryID            = $args{queryID}    || "";
    
    my $title              = $args{title}      || "annotations maps";
    
    # Llama a cluster binary en local
    my $_gff2ps_dir    = "/home/ug/gmaster/projects/gff2ps";
    my $_gff2ps_bin    = "bin/gff2ps";
    my $_gff2ps_params = "params/matrices.gff2ps.param";
    my $_gff2ps_data   = "data/empty.gff";
    my $_gff2ps_args   = "-C $_gff2ps_dir/$_gff2ps_params -T \"$title\" -loOa";
    
    # Check that the binary is in place
    if (! -f "$_gff2ps_dir/$_gff2ps_bin") {
	my $note = "Internal System Error. $_gff2ps_bin script not found";
	print STDERR "$note\n";
	my $code = 701;
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    
    # Create the temp GFF maps input files
    
    my @gff_files    = ();
    my $gff2ps_input = "";
    foreach my $gff_map (@$gff_maps) {
    
      if (length ($gff2ps_input) != 0) {
        $gff2ps_input .= " $_gff2ps_dir/$_gff2ps_data ";
      }
    
      my ($gff_map_fh, $gff_map_file);
      eval {
	($gff_map_fh, $gff_map_file) = tempfile("/tmp/GFF_MAP.XXXXXX", UNLINK => 0);
      };
      if ($@) {
	my $note = "Internal System Error. Can not open GFF map input temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
      }
      print $gff_map_fh "$gff_map\n";
      close $gff_map_fh;
      push (@gff_files, $gff_map_file);
      $gff2ps_input .= "$gff_map_file";
    }
    
    # Setup temporary output files
    
    my ($ps_fh, $ps_file);
    eval {
	($ps_fh, $ps_file) = tempfile("/tmp/GFF2PS_PS.XXXXXX", UNLINK => 0);
    };
    if ($@) {
	my $note = "Internal System Error. Can not open gff2ps postcript output temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    close $ps_fh;
    # rename the file so it has the extension .ps !!"
    # otherwise the conversion into jpeg will not work !!
    qx/mv $ps_file $ps_file.ps/;
    $ps_file = $ps_file . ".ps";
    
    my ($jpeg_fh, $jpeg_file);
    eval {
	($jpeg_fh, $jpeg_file) = tempfile("/tmp/GFF2PS_JPEG.XXXXXX", UNLINK => 0);
    };
    if ($@) {
	my $note = "Internal System Error. Can not open convert jpeg output temporary file!\n";
	my $code = 701;
	print STDERR "$note\n";
	my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
	return (undef, [$moby_exception]);
    }
    close $jpeg_fh;
    # rename the file so it has the extension .jpg !!"
    qx/mv $jpeg_file $jpeg_file.jpg/;
    $jpeg_file = $jpeg_file . ".jpg";
    
    if ($debug) {
      print STDERR "executing gf2ps with the following command:\n";
      print STDERR "$_gff2ps_dir/$_gff2ps_bin $_gff2ps_args $gff2ps_input > $ps_file\n";
    }
    
    my $gff2ps_result = qx/$_gff2ps_dir\/$_gff2ps_bin $_gff2ps_args $gff2ps_input > $ps_file/;
    
    if ($debug) {
      print STDERR "gff2ps execution done\n\n";
      print STDERR "converting into JPEG format the postcript file, with the following command:\n";
      print STDERR "convert -rotate 90 $ps_file $jpeg_file\n";
    }
    
    # I expect convert tool to be in the PATH !
    my $convert_result = qx/convert -rotate 90 $ps_file $jpeg_file/;
    
    if ($debug) {
      print STDERR "conversion into JPEG format done\n\n";
    }
    
    # if not empty file
    if (! -z $jpeg_file) {
      $image = qx/cat $jpeg_file/;
    }
    
    if (!$debug) {
      foreach my $map_file (@gff_files) {
        unlink $map_file;
      }
      unlink $ps_file;
      unlink $jpeg_file;
    }
    
    if (! defined $image) {
      my $note = "Internal System Error. gff2ps has failed!\n";
      my $code = 701;
      print STDERR "$note\n";
      my $moby_exception = INB::Exceptions::MobyException->new (
								  code       => $code,
								  type       => 'error',
								  queryID    => $queryID,
								  message    => "$note",
								  );
      return (undef, [$moby_exception]);
    }
    
    return ($image, $moby_exceptions);
}

#####################################################################################################

sub _parse_ace_file {
    my ($ace_file_name, $debug, %ace_parsing_options) = @_;
    
    if ($debug) {
	print STDERR "parsing the ace file, $ace_file_name...\n";
    }
    
    # Open the ace file
    
    open (ACE, "< $ace_file_name") or die "can't open the ace file, $ace_file_name !!!\n";
    
    # parse it
    
    # this data structure stores everything !!!!
    # this array will be returned to main
    my @clusters_list = ();
    # CLusters of one member are returned as being singlets
    my @singlets = ();
    
    while (my $line = <ACE>) {
	
	chomp $line;
	if ($debug) {
	    print STDERR "processing line, $line...\n";
	}
	
	if ($line =~ /$ace_parsing_options{contig_pattern}/) {
	    my $contig_name = $1;
	    if ($debug) {
		print STDERR "new cluster, contig name : $contig_name\n";
	    }
	    # retrieve Contig sequence
	    my $sequence_contig = "";
	    if ($debug) {
		print STDERR "parsing contig sequence...\n";
	    }
	    my $new_line = <ACE>;
	    while (lc ($new_line) =~ /(a|c|t|g|n)\S+/) {
		chomp ($new_line);
		$sequence_contig = $sequence_contig . $new_line;
		$new_line = <ACE>;
	    }
	    if ($debug) {
		print STDERR "sequence_contig : $sequence_contig\n";
	    }
	    
	    $sequence_contig =~ s/\*//g;
	    
	    # then Base quality
	    while (not $new_line =~ /$ace_parsing_options{AF_test}/) {
		$new_line = <ACE>;
	    }
	    my @members = ();
	    # then information about the number of members within the current cluster
	    while ($new_line =~ /$ace_parsing_options{AF_test}/) {
		$new_line =~ /$ace_parsing_options{AF_pattern}/;
		my $member_name = $1;
		
		if (! defined $member_name) {
		    print STDERR "error, can't get the member name from line, $new_line!\n";
		}
		
		if (not _is_in ($member_name, \@members)) {
		    push (@members, $member_name);
		}
		
		if ($debug) {
		    print STDERR "found member $member_name, members list : @members, new_line : $new_line\n";
		}
		
		$new_line = <ACE>;
	    }
	    my $members_nb = (@members+0);
	    my $cpt = 0;
	    my @members_list = ();
	    
	    if ($debug) {
		print STDERR "found $members_nb members for the current cluster\n";
	    }
	    
	    if ($members_nb > 1) {
	    
		# Retrieve the name and the sequence of each member
		
		while ($cpt < $members_nb) {
		    if ($new_line =~ /$ace_parsing_options{member_sequence}/) {
			my $member_name = $1;
			if ($debug) {
			    print STDERR "member name : $member_name\n";
			}
			my $new_line = <ACE>;
			my $sequence = "";
			while (lc ($new_line) =~ /(a|c|t|g|n)\S+/) {
			    chomp ($new_line);
			    $sequence = $sequence . $new_line;
			    $new_line = <ACE>;
			}
			if ($debug) {
			    print STDERR "member sequence : $sequence\n";
			}
			
			# the sequence may contain * which don't fit with the fasta format, so remove them from the sequence
			# this * is for a shift in the alignment
			
			$sequence =~ s/\*//g;
			
			my %h_member = (line_info => "$member_name",
					sequence  => "$sequence",
					);
			$cpt ++;
			push (@members_list, \%h_member);
		    }
		    $new_line = <ACE>;
		}
		
		my $line_info = $contig_name . " $members_nb members: " . join ( ", ", @members);
		
		my %h_contig = (line_info => "$line_info",
				sequence  => "$sequence_contig",
				);
		
		# this data structure stores all the information about one contig and its members
		my %cluster_info = (cluster_info => \%h_contig,
				    members_info => \@members_list,
				    );
		push (@clusters_list, \%cluster_info);
		
		if ($debug) {
		    print STDERR "Contig done - current line : $new_line\n";
		}
	    }
	    else {
		# Only one member, so it is a singlet !!
		
		my $member_id = $members[0];
		
		if ($debug) {
		    print STDERR "singlet, $member_id\n";
		}
		
		push (@singlets, $member_id);
	    }
	}
    }
    close (ACE);
    
    if ($debug) {
	print STDERR "parsing $ace_file_name done\n";
	print STDERR @clusters_list . " clusters found\n";
	print STDERR @singlets . " singlets found\n";
	
    }
    
    return (\@clusters_list, \@singlets);
    
}

sub _is_in {
    my ($member_name, $members) = @_;
    
    foreach my $element (@$members) {
	if ($element eq $member_name) {
	    return 1;
	}
    }
    
    return 0;
}

1;

__END__
