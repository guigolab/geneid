# $Id: FactoryLog.pm,v 1.1 2007-07-15 12:34:59 arnau Exp $
#
# INBPerl module for INB::GRIB::Services::FactoryLog
#
# This file is an instance of a template written
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#

# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::FactoryLog

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

package INB::GRIB::Services::FactoryLog;

use strict;
use warnings;
use Data::Dumper;

use Carp;

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

# Aqui pondremos las funciones que exportaremos a otros módulos. Ejemplo
#
#  @EXPORT = qw( &func1 &func2);
#
our @EXPORT = qw(
  &LogReport_call
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
	
	use vars qw /$gmaster_home/;
	$gmaster_home = $ENV{HOME};
	# $gmaster_home = "/home/ug/gmaster";
	
}

###############################################################################

sub LogReport_call {
  my %args = @_;

  # relleno los parametros por defecto GeneID_call

  my $debug      = $args{debug}  || 0;
  my $queryID    = $args{queryID} || "";
  my $parameters = $args{parameters} || undef;
	
  # Get the parameters
	
  my $startTime    = $parameters->{startTime};
  my $endTime      = $parameters->{endTime};
  my $includeTests = $parameters->{includeTests};

  if ($debug) {	
    print STDERR "LogReport parameters (startTime, endTime, includeTests): $startTime, $endTime, $includeTests.\n";
  }
	
  # Output definition
	
  my $logEvents_href = {};
  my $moby_exceptions = [];

  # ...
	
  return ($logEvents_href, $moby_exceptions);
  
}

1;

__END__
