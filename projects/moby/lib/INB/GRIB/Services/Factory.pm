# $Id: Factory.pm,v 1.68 2006-03-08 14:13:03 gmaster Exp $
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
  # versi�n es la siguiente:

 @params = ('nucleotide' => "ATCAGAGAGCGAGCAGATGCAGAGTAGAGCGAGCGAGCGCT");
 $report = INB::GRIB::Services::Factory::factory_call(@params);

  # Esta llamada nos devuelve una variable que contiene el texto con la
  # salida del programa geneid.

=head1 DESCRIPTION

Este package sirve para hacer la llamadas al programa geneid.
Se han utilizado las libreriras de geneid, por tanto ser�n
necesarias instalarlas previamente.

Es requisito necesario para ejecutar este m�dulo, que las variables de
entorno DIR y DATADIR est�n descritas en el sistema. Dentro de
la variable $DIR se encuentra el path absoluto a la aplicacion geneid,
mientras que en $DATADIR pondremos la localizaci�n al directorio que
contiene las bases de datos.

En siguientes versiones, en este modulo se podran a�adir las diferentes
funciones para llamar a los programas que se requieran.

=head1 AUTHOR

Francisco Camara, fcamara@imim.es

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

# Create temporary data files
use File::Temp qw/tempfile/;

# Bioperl

use Bio::SeqIO;
use Bio::Seq;

# Moby Exceptions module
use INB::Exceptions::MobyException;

# Report Unix Error codes
use Errno qw (EINTR EIO :POSIX);

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw() ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

# Aqui pondremos las funciones que exportaremos a otros m�dulos. Ejemplo
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
  &MetaAlignment_call
  &generateScoreMatrix_call
  &MEME_call
  &meme2matrix_call
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
 Returns : Devuelve un string que contiene el resultado de la ejecuci�n.

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
 Returns : Devuelve un string que contiene el resultado de la ejecuci�n.

=cut

sub GeneID_call {
    my %args = @_;

    # output specs
    my $geneid_output   = "";
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
    my $_geneid_bin  = "bin/geneid";
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
	return (undef, [$moby_exception]);
    }

    if ($format eq "GFF") {
	$_geneid_args .= "-G";
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
	return (undef, [$moby_exception]);
    }

    # print STDERR "Running GeneID, with this command:\n";
    # print STDERR "$_geneid_dir\/$_geneid_bin $_geneid_args $seqfile \n";

    $geneid_output = qx/$_geneid_dir\/$_geneid_bin $_geneid_args $seqfile/;

    # Comment this line if you want to keep the file...
    unlink $seqfile;

    if (defined $geneid_output) {
	return ($geneid_output, $moby_exceptions);
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
	return (undef, $moby_exceptions);
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
 Returns : Devuelve un string que contiene el resultado de la ejecuci�n.

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

	# No parameters yet !!

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
	    close($blast_fh);
	    qx/echo "$tblastx_output" > $tblastx_output_file/;
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

	# print STDERR "Running SGP2, with this command:\n";
	# print STDERR "$_sgp2_dir\/$_sgp2_bin $_sgp2_args -1 $seqfile -t $tblastx_output_file\n";

	$sgp2_output = qx/$_sgp2_dir\/$_sgp2_bin $_sgp2_args -1 $seqfile -t $tblastx_output_file/;

	# Comment these two lines if you want to keep the file...
	unlink $seqfile;
	unlink $tblastx_output_file;

	if (defined $sgp2_output) {
	    return ($sgp2_output, $moby_exceptions);
	}
	else {
	    my $note = "Internal System Error. SGP2 has failed!\n";
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


=head2 GOstat_call

 Title   : GOstat_call
 Usage   : $report = GOstat_call (@params);
	 :
	 : ## where @params are,
	 : @params = ('arg1'  => "WKRPPEICENPRFIIGGANRTDIAAIACLTLNERL",
	 :            'arg2'  => "query 1", ## optional
	 :            'arg3'  => "nr");     ## optional (default: nr)
 Returns : Devuelve un string que contiene el resultado de la ejecuci�n.

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

    # relleno los parametros por defecto GeneID_call

    my $genes_ref  = $args{genes}      || undef;
    my $parameters = $args{parameters} || undef;
    my $queryID    = $args{queryID}    || "";

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

    $_promExtraction_args .= " -r $dbrelease" || die "no ensembl release was given!\n";

    # Make a temporary file for the input list of genes

    my ($genes_list_fh, $genes_list_file) = tempfile("/tmp/PROM_EXTRACTION_GENES.XXXXXX", UNLINK => 1);
    close ($genes_list_fh);

    open (FILE, ">$genes_list_file") or die "can't open temp file, $genes_list_file!\n";
    print FILE (join ("\n", @$genes_ref) . "\n");
    close FILE;

    # print STDERR "running command,\n";
    # print STDERR "$_promExtraction_dir\/$_promExtraction_bin $_promExtraction_args -f $genes_list_file\n";

    $promoterExtraction_output = qx/$_promExtraction_dir\/$_promExtraction_bin $_promExtraction_args -f $genes_list_file/;

    unlink $genes_list_file;

    if (defined $promoterExtraction_output) {
	return $promoterExtraction_output;
    }
    else {
	return undef;
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
 Returns : Devuelve un string que contiene el resultado de la ejecuci�n.

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
    
    $matscan_output = qx/$_matscan_dir\/$_matscan_bin $_matscan_args $seqfile $matrix_file | grep MatScan/;
    
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

sub generateScoreMatrix_call {
  my %args = @_;

  # output specs declaration
  my $matrix_output   = "";
  my $moby_exceptions = [];
  
  # relleno los parametros por defecto generateScoreMatrix_call
  
  my $inputdata_arrayref = $args{similarity_results};
  my $parameters         = $args{parameters} || undef;
  my $queryID            = $args{queryID}    || "";
  my $debug = 0;
  
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
    my $_meme_dir   = "/usr/local/molbio/Install/meme-3.5.0";
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

    if ($format =~ /text/i) {
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
	my $note = "Internal System Error. the parsing of MEME data has failed!\n";
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

1;

__END__
