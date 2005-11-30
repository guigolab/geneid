# $Id: Factory.pm,v 1.49 2005-11-30 16:14:53 gmaster Exp $
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
	my $mode    = "normal";

	# print STDERR "GeneID parameters (profile, strands, format): $profile, $strands, $format.\n";

	my $results = "";

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
	
	return $results;
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

        # relleno los parametros por defecto GeneID_call

        my $sequences          = $args{sequences} || undef;
	my $format             = $args{format} || "";
	my $parameters         = $args{parameters} || undef;

	# Get the parameters

	my $profile     = $parameters->{profile};
	my $strands     = $parameters->{strands};
	my $exons_ref   = $parameters->{exons};
	my $signals_ref = $parameters->{signals};

        # Llama a GeneID en local
        my $_geneid_dir  = "/home/ug/gmaster/GeneID/geneid_2002";
        my $_geneid_bin  = "bin/geneid";
        my $_geneid_args = "";
	
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
	    if ($profile eq "Human")         { $_geneid_args .= "P $_geneid_dir/human.param"; last SWITCH; }
	    if ($profile eq "Tetraodon")     { $_geneid_args .= "P $_geneid_dir/tetraodon.param"; last SWITCH; }
	    if ($profile eq "Drosophila")    { $_geneid_args .= "P $_geneid_dir/dros.param"; last SWITCH; }
	    if ($profile eq "Celegans")      { $_geneid_args .= "P $_geneid_dir/celegans.param"; last SWITCH; }
	    if ($profile eq "Wheat")         { $_geneid_args .= "P $_geneid_dir/wheat.param"; last SWITCH; }
	    if ($profile eq "Arabidopsis")   { $_geneid_args .= "P $_geneid_dir/arabidopsis.param"; last SWITCH; }
	    if ($profile eq "Rice")          { $_geneid_args .= "P $_geneid_dir/rice.param"; last SWITCH; }
	    if ($profile eq "Plasmodium")    { $_geneid_args .= "P $_geneid_dir/plasmodium.param"; last SWITCH; }
	    if ($profile eq "Dictyostelium") { $_geneid_args .= "P $_geneid_dir/dictyostelium.param"; last SWITCH; }
	    if ($profile eq "Aspergillus")   { $_geneid_args .= "P $_geneid_dir/aspergillus.param"; last SWITCH; }
	    if ($profile eq "Neurospora")    { $_geneid_args .= "P $_geneid_dir/neurospora.param"; last SWITCH; }
	    if ($profile eq "Cryptococcus")  { $_geneid_args .= "P $_geneid_dir/cneomorfans.param"; last SWITCH; }
	    if ($profile eq "Coprinus")      { $_geneid_args .= "P $_geneid_dir/cinereus.param"; last SWITCH; }
	    # Default is Human
	    $_geneid_args .= "P $_geneid_dir/human.param";
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

	# Generate a temporary file locally with the sequence(s) in FASTA format
	# locally, ie not on a NFS mounted directory, for speed sake

	my ($seq_fh, $seqfile) = tempfile("/tmp/GENEID.XXXXXX", UNLINK => 0);

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
	    print STDERR "Error, empty sequence file...\n";
	}

	# print STDERR "Running GeneID, with this command:\n";
	# print STDERR "$_geneid_dir\/$_geneid_bin $_geneid_args $seqfile \n";

        my $geneid_output = qx/$_geneid_dir\/$_geneid_bin $_geneid_args $seqfile/;
        
	# Comment this thine if you want to keep the file...
	unlink $seqfile;

        if (defined $geneid_output) {
		return $geneid_output;
	}	
	else {
		# What else better to return ??
		return undef;
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
        # relleno los parametros por defecto SGP2_call (nucleotide  => $nucleotide, seqIdentifier => $sequenceIdentifier);
        my $sequences          = $args{sequences}      || undef;
        my $tblastx_output     = $args{tblastx_output} || undef;
	my $format             = $args{format}         || "";
	my $parameters         = $args{parameters}     || undef;

	# No parameters yet !!

        # Llama a SGP2 localmente
        my $_sgp2_dir  = "/home/ug/gmaster/sgp2/sgp2_2003/";
	$_sgp2_dir     = $ENV{SGP2};
        my $_sgp2_bin  = "bin/sgp2";
        my $_sgp2_args = "";
	
        if ($format eq "GFF") {
	    $_sgp2_args .= "-g G";
        }
	
	# Generate a temporary file locally with the sequence in FASTA format
	# locally, ie not on a NFS mounted directory, for speed sake

	my ($fh1, $seqfile) = tempfile("/tmp/SGP2_Sequence.XXXXXX", UNLINK => 0);
	close($fh1);

	my @seqIds = keys (%$sequences);
	foreach my $sequenceIdentifier (@seqIds) {

	    my $nucleotides = $sequences->{$sequenceIdentifier};
	    
	    # bioperl object
	    
	    my $seqobj = Bio::Seq->new (
					-display_id => $sequenceIdentifier,
					-seq        => $nucleotides
					);
	    
	    # Bioperl sequence factory
	    
	    my $sout = Bio::SeqIO->new (
	     			    -file   => ">$seqfile",
	    			    -format => 'fasta'
	    			    );
	    $sout->write_seq ($seqobj);
	    
	}
	
	# TBLASTX Output File

	# Generate a temporary file locally with the TBLASTX Output
	# locally, ie not on a NFS mounted directory, for speed sake

	my ($fh2, $tblastx_output_file) = tempfile("/tmp/SGP2_TBLASTX.XXXXXX", UNLINK => 0);
	close($fh2);
	
	qx/echo "$tblastx_output" > $tblastx_output_file/;
	
	# Test empty files

	if (-z $seqfile) {
	    print STDERR "Error, empty sequence file...\n";
	    exit 0;
	}

	if (-z $tblastx_output_file) {
	    print STDERR "Error, empty tblastx output file...\n";
	    exit 0;
	}

	# print STDERR "Running SGP2, with this command:\n";
	# print STDERR "$_sgp2_dir\/$_sgp2_bin $_sgp2_args -1 $seqfile -t $tblastx_output_file\n";
	
        my $sgp2_output = qx/$_sgp2_dir\/$_sgp2_bin $_sgp2_args -1 $seqfile -t $tblastx_output_file/;
        
	unlink $seqfile;
	unlink $tblastx_output_file;
	
        if (defined $sgp2_output) {
	    return $sgp2_output;
	}	
	else {
	    # What else better to return ??
	    return undef;
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

    my $regulated_genes    = $args{regulated_genes} || undef;
    my $reference_genes    = $args{reference_genes} || undef;
    my $format             = $args{format}          || "";
    my $parameters         = $args{parameters}      || undef;
    
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
    
    my $gostat_output = qx/$_gostat_dir\/$_gostat_bin $_gostat_args --reg $regulated_genes_file --ref $reference_genes_file/;
        
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

    my $sequences          = $args{sequences}    || undef;
    my $geneid_predictions = $args{predictions}  || undef;
    my $parameters         = $args{parameters}   || undef;
    
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
    
    # Make two temporary files for both input lists of genes
    
    my ($seq_fh, $seqfile) = tempfile("/tmp/SEQS.TRANSLATION.XXXXXX", UNLINK => 0);
    # close ($seq_fh);

    my ($feature_fh, $featurefile) = tempfile("/tmp/FEATURES.TRANSLATION.XXXXXX", UNLINK => 0);
    close ($feature_fh);
    open (FILE, ">$featurefile") or die "can't open temp file, $featurefile!\n";

    # Bioperl sequence factory
    
    my $sout = Bio::SeqIO->new (
				-fh     => $seq_fh,
				-format => 'fasta'
				);
    
    my @seqIds = keys (%$sequences);
    
    foreach my $sequenceIdentifier (@seqIds) {
	
	# Sequence

	my $nucleotides = $sequences->{$sequenceIdentifier};
	
	# bioperl object
	
	my $seqobj = Bio::Seq->new (
				    -display_id => $sequenceIdentifier,
				    -seq        => $nucleotides
				    );
	
	$sout->write_seq ($seqobj);

	# GeneID Predictions

	my $geneid_prediction = $geneid_predictions->{$sequenceIdentifier};

	print FILE "$geneid_prediction";
	
    }
    
    close $seq_fh;
    close FILE;

    # Test empty file
    if (-z $seqfile) {
	print STDERR "Error, empty sequence file...\n";
    }

    # print STDERR "Running the following command, $_translateGeneID_dir\/$_translateGeneID_bin $_translateGeneID_args -s $seqfile -f $featurefile...\n";
    
    my $translateGeneID_output = qx/$_translateGeneID_dir\/$_translateGeneID_bin $_translateGeneID_args -s $seqfile -f $featurefile/;
        
    unlink $seqfile;
    unlink $featurefile;

    if (defined $translateGeneID_output) {
	return $translateGeneID_output;
    }
    else {
	# What else better to return ??
	return undef;
    }

}

sub PromoterExtraction_call {
    my %args = @_;

    # relleno los parametros por defecto GeneID_call
    
    my $genes_ref  = $args{genes}      || undef;
    my $parameters = $args{parameters} || undef;
    
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
    
    my $promoterExtraction_output = qx/$_promExtraction_dir\/$_promExtraction_bin $_promExtraction_args -f $genes_list_file/;
    
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
 Returns : Devuelve un string que contiene el resultado de la ejecución.

=cut

sub MatScan_call {
        my %args = @_;
	
        # relleno los parametros por defecto MatScan_call
	
        my $sequences    = $args{sequences}  || undef;
	my $matrix_input = $args{matrix}     || undef;
	my $format       = $args{format}     || "";
	my $parameters   = $args{parameters} || undef;
	my $debug        = $args{debug}      || 0;
	
	# Get the parameters
	
	my $threshold        = $parameters->{threshold};
	my $strands          = $parameters->{strands};
	my $matrix_parameter = $parameters->{matrix};
	my $matrix_mode      = $parameters->{matrix_mode};
	
        # Llama a MatScan en local
        my $_matscan_dir  = "/home/ug/gmaster/projects/Meta/";
        my $_matscan_bin  = "bin/matscan";
        my $_matscan_args = "-T $threshold";
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
	
	if (defined $matrix_parameter) {

	    if ($debug) {
	      print STDERR "matrix as a parameter...\n";
            }

	    if ($matrix_mode eq "raw format") {

	        if ($debug) {
  		  print STDERR "raw mode\n";
                }
                
	      SWITCH: {
		  if ($matrix_parameter eq "Transfac") { $matrix_file = "$_matscan_dir/matrices/Transfac_raw_format.matrices"; last SWITCH; }
		  if ($matrix_parameter eq "MEME")     { $matrix_file = "$_matscan_dir/matrices/Promo_raw_format.matrices"; last SWITCH; }
		  if ($matrix_parameter eq "Jaspar")   { $matrix_file = "$_matscan_dir/matrices/Jaspar_raw_format.matrices"; last SWITCH; }
		  # Default is Transfac
		  $matrix_file = "$_matscan_dir/matrices/Transfac_raw_format.matrices";
	      }
	    }
	    elsif ($matrix_mode eq "log-likelihood") {
	      if ($debug) {	
  	        print STDERR "log-likelihood mode\n";
              }
	
	      SWITCH: {
		  if ($matrix_parameter eq "Transfac") { $_matscan_args .= " -sm"; $matrix_file = "$_matscan_dir/matrices/Transfac_likelihood.matrices"; last SWITCH; }
		  if ($matrix_parameter eq "MEME")     { $_matscan_args .= " -sm"; $matrix_file = "$_matscan_dir/matrices/Promo_likelihood.matrices"; last SWITCH; }
		  if ($matrix_parameter eq "Jaspar")   { $_matscan_args .= " -sm"; $matrix_file = "$_matscan_dir/matrices/Jaspar_likelihood.matrices"; last SWITCH; }
		  # Default is Transfac
		  $_matscan_args .= " -sm";
		  $matrix_file = "$_matscan_dir/matrices/Transfac_likelihood.matrices";
	      }
	    }
	    else {
		print STDERR "don't know anything about matrix mode, $matrix_mode!\n";
		exit 0;
	    }
	}
	elsif (defined $matrix_input) {
	    if ($debug) {
	      print STDERR "matrix as an input...\n";
	    }
	    
	    # Make a temporary file with the matrix input
	    $_matscan_args .= " -sm";
	    ($matrix_fh, $matrix_file) = tempfile("/tmp/MATSCAN_MATRIX.XXXXXX", UNLINK => 0);
	    print $matrix_fh "$matrix_input";
	    close $matrix_fh;
	}
	else {
	    print STDERR "matrix_input neither matrix_parameter are defined!!\n";
	    exit 0;
	}
	
	if (not defined $matrix_file) {
	    print STDERR "Error, no defined matrix file!\n";
	    exit 0;
	}

	# Generate a temporary file locally with the sequence(s) in FASTA format
	# locally, ie not on a NFS mounted directory, for speed sake
	
	my ($seq_fh, $seqfile) = tempfile("/tmp/MATSCAN_SEQS.XXXXXX", UNLINK => 0);
	
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
	    print STDERR "Error, empty sequence file...\n";
	}
	
	if ($debug) {
	    print STDERR "Running Matscan, with this command:\n";
	    print STDERR "$_matscan_dir\/$_matscan_bin $_matscan_args $seqfile $matrix_file\n";
	}
	
        my $matscan_output = qx/$_matscan_dir\/$_matscan_bin $_matscan_args $seqfile $matrix_file | grep MatScan/;
        
	# Comment this line if you want to keep the file...
	unlink $seqfile;
	if (defined $matrix_input) {
	    unlink $matrix_file;
	}
	
        if (defined $matscan_output) {
		return $matscan_output;
	}
	else {
		# What else better to return ??
		return undef;
	}
}


sub MetaAlignment_call {
    my %args = @_;  

    # method output
    my ($meta_output, $note, $code) = ("", undef, undef);
    
    # relleno los parametros por defecto MetaAlignment_call

    my $map1 = $args{map1};
    my $map2 = $args{map2};
    my $parameters = $args{parameters} || undef;

    # Get the parameters

    my $alpha_penalty  = $parameters->{alpha_penalty};
    my $lambda_penalty = $parameters->{lambda_penalty};
    my $mu_penalty     = $parameters->{mu_penalty};
    
    my $output_format  = $parameters->{output_format};
    
    # Llama a Meta-alignment en local
    my $_meta_alignment_dir  = "/home/ug/gmaster/projects/Meta";
    my $_meta_alignment_bin  = "bin/meta";
    my $_meta_alignment_args = "-a $alpha_penalty -l $lambda_penalty -m $mu_penalty";

    # Check that the binary is in place
    if (! -f "$_meta_alignment_dir/$_meta_alignment_bin") {
	$note = "meta-alignment binary not found";
	print STDERR "$note\n";
	$code = 701;
	return ("", $note, $code);
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
	$note = "can't open meta-alignment input temporary file!\n";
	$code = 701;
	print STDERR "$note\n";
	return ("", $note, $code);
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
	$note = "can't open meta-alignment input temporary file!\n";
	$code = 701;
	print STDERR "$note\n";
	return ("", $note, $code);
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
	$note = "can't open a temporary file to store meta-alignment output!\n";
	$code = 701;
	print STDERR "$note\n";
	return ("", $note, $code);
    }
    
    my @args = ("$_meta_alignment_dir\/$_meta_alignment_bin $_meta_alignment_args $map1_file $map2_file > $stdout_file");
    
    my $failed = system (@args);
    if ($failed > 0) {
	$note = "meta-alignment system call died (with error code, $?).\n";
	$code = 701;
	
	if (($! != ENOTTY) || ($! ne "Inappropriate ioctl for device")) {
	    # This is not an error, just mean that stdout is a terminal !!
	    print STDERR "Error, '$!'\n";
	}
	print STDERR $note;
    }
    else {
	$meta_output = qx/cat $stdout_file/;
    }
    
    unlink $stdout_file;
    unlink $map1_file;
    unlink $map2_file;
    
    return ($meta_output, $note, $code);
    
}

sub generateScoreMatrix_call {
  my %args = @_;  
  
  # relleno los parametros por defecto generateScoreMatrix_call
  
  my $inputdata_arrayref = $args{similarity_results};
  my $parameters         = $args{parameters} || undef;
  
  # Get the parameters
  
  my $input_format   = $parameters->{input_format};
  my $output_format  = $parameters->{output_format};
  
  # Llama a generateScoreMatrix script en local
  my $_application_dir  = "/home/ug/gmaster/projects/generateScoreMatrices";
  my $_application_bin  = "generateScoreMatrix.pl";
  my $_application_args = "$input_format";
  
  # Generate a temporary file

  my ($meta_fh, $meta_file) = tempfile("/tmp/META_OUTPUT.XXXXXX", UNLINK => 0);
  print $meta_fh "@$inputdata_arrayref\n";
  close $meta_fh;
  
  # Run generateScoreMatrix
  
  # print STDERR "got " . @$inputdata_arrayref . " array references (meta alignments)\n";
  # print STDERR "Running generateScoreMatrix.pl, with this command:\n";
  # print STDERR "cat $meta_file | $_application_dir\/$_application_bin $_application_args\n";

  my $matrix = qx/cat $meta_file | $_application_dir\/$_application_bin $_application_args/;

  unlink $meta_file;
  
  return $matrix;

}

sub MEME_call {
    my %args = @_;
    
    my $sequences  = $args{sequences};
    my $format     = $args{format};
    my $parameters = $args{parameters};
    my $_debug     = $args{debug};

    my $motif_distribution    = $parameters->{motif_distribution};
    my $maximum_number_motifs = $parameters->{maximum_number_motifs};
    my $minimum_number_sites  = $parameters->{minimum_number_sites};
    my $maximum_number_sites  = $parameters->{maximum_number_sites};
    my $minimum_motif_width   = $parameters->{minimum_motif_width};
    my $maximum_motif_width   = $parameters->{maximum_motif_width};
    my $e_value_cutoff        = $parameters->{e_value_cutoff};
    
    # Llama a Meme en local
    my $_meme_dir   = "/usr/local/molbio/Install/meme-3.5.0";
    my $_meme_bin  = "bin/meme";
    my $_meme_args = "-nostatus -time 160 -maxiter 20 ";
    
    # Setting up the MEME parameters
    
    # Default is HTML
    
    if ($_debug) {
	print STDERR "format, $format\n";
    }
    
    if ($format eq "text-formatted") {
	$_meme_args .= "-text ";
    }
    
    if (not defined $motif_distribution) {
	$_meme_args .= "-mod zoops";
    }
    else {
      SWITCH: {
	  if ($motif_distribution eq "zero or one")               { $_meme_args .= "-mod zoops"; last SWITCH; }
	  if ($motif_distribution eq "one")                       { $_meme_args .= "-mod ops"; last SWITCH; }
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
    
    my ($seq_fh, $seqfile) = tempfile("/tmp/MEME.XXXXXX", UNLINK => 0);
    
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
	
	# bioperl object
	
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
	print STDERR "Error, empty sequence file...\n";
    }
    
    # MEME default is protein alphabet so if it is dna sequences, specify it, as the MEME model will be then different
    
    if (defined $alphabet && (($alphabet eq "dna") || ($alphabet eq "rna"))) {
	$_meme_args .= " -dna -revcomp";
    }

    if ($_debug) {
	print STDERR "Running Meme, with this command:\n";
	print STDERR "$_meme_dir\/$_meme_bin $seqfile $_meme_args\n";
    }
    
    my $meme_output = qx/source $_meme_dir\/etc\/meme.sh; $_meme_dir\/$_meme_bin $seqfile $_meme_args/;
    
    # Comment this line if you want to keep the file...
    unlink $seqfile;
    
    if (defined $meme_output) {
	return $meme_output;
    }	
    else {
	# What else better to return ??
	return undef;
    }
}


sub meme2matrix_call {
    my %args = @_;
    
    # relleno los parametros por defecto meme2matrix_call
    
    my $meme_predictions   = $args{meme_predictions} || undef;
    my $parameters         = $args{parameters}       || undef;
    my $_debug             = $args{debug};
    
    # Get the parameters
    
    my $matrix_mode = $parameters->{matrix_mode};
    
    # Llama a Meme2matrix en local
    my $_meme2matrix_dir  = "/home/ug/gmaster/projects/meme2matrix";
    my $_meme2matrix_bin  = "meme2matrix.pl";
    my $_meme2matrix_args = "";
    
    if ($matrix_mode eq "raw format") {
	$_meme2matrix_args = "score";
    }
    elsif ($matrix_mode eq "log-likelihood") {
	$_meme2matrix_args = "probability";
    }

    # Create the temp map files

    my ($meme2matrix_fh, $meme2matrix_file) = tempfile("/tmp/MEME2MATRIX.XXXXXX", UNLINK => 0);
    close ($meme2matrix_fh);
    
    open (FILE, ">$meme2matrix_file") or die "can't open MEME2MATRIX temp file, $meme2matrix_file!\n";
    print FILE $meme_predictions;
    close FILE;
    
    if ($_debug) {
	print STDERR "Running meme2matrix, with this command:\n";
	print STDERR "$_meme2matrix_dir\/$_meme2matrix_bin $meme2matrix_file $_meme2matrix_args\n";
    }
    
    my $matrices_output = qx/$_meme2matrix_dir\/$_meme2matrix_bin $meme2matrix_file $_meme2matrix_args/;
    
    # Comment this line if you want to keep the file...
    unlink $meme2matrix_file;
    
    my @matrices = ();
    
    # Return an array of matrices from the output (which is a string)
    # Because we want to return a collection of text-formatted objects (each object containing a matrix)
    # So better to return an array rather than just a string
    
    # To split the matrices, check which position they start, ie chech th position of 'MEME'
    my @positions = ();
    while ($matrices_output =~ m/MEME/g) {
	push (@positions, pos ($matrices_output) - 4);
    }
    
    if ($_debug) {
	print STDERR "positions: @positions\n";
    }
    
    if (@positions == 0) {
	return undef;
    }
    
    my $position_1 = $positions[0];
    my $index = 1;
    while (my $position_2 = $positions[$index]) {
	my $matrix = substr ($matrices_output, $position_1, $position_2 - $position_1);
	push (@matrices, $matrix);
	$position_1 = $position_2;
	$index++;
    }
    my $end = length ($matrices_output);
    my $matrix = substr ($matrices_output, $position_1, $end - $position_1);
    push (@matrices, $matrix);
    
    if ($_debug) {
	print STDERR "nb matrices: " . @matrices . ".\n";
        print STDERR "Dumping matrices array,\n" . Dumper (@matrices) . "\n";
    }
    
    if (@matrices > 0) {
	return \@matrices;
    }	
    else {
	# What else better to return ??
	return undef;
    }
    
}

1;

__END__

