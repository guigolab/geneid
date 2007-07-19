# $Id: FactoryLog.pm,v 1.5 2007-07-19 20:18:57 arnau Exp $
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

# Date Format
use DateTime;
use DateTime::Format::W3CDTF;

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

sub _W3C_to_inb {
  my ( $w3c_time ) = @_;

  # W3C -> INB
  my $f = DateTime::Format::W3CDTF->new;
  my $dt = $f->parse_datetime( $w3c_time );
  my $time_inb = int(
    (sprintf "%04d",$dt->year) . (sprintf "%02d",$dt->month) . (sprintf "%02d",$dt->day) . (sprintf
    "%02d",$dt->hour) . (sprintf "%02d",$dt->minute) . (sprintf "%02d",$dt->second)
    );
  
  return $time_inb;
}
  
sub _inb_to_W3C {
  my ( $time_inb ) = @_;

  $time_inb =~ /(\d\d\d\d)(\d\d)(\d\d)(\d\d)(\d\d)(\d\d)/;
  my $dt = DateTime->new(
    year=> $1,
    month  => $2,
    day => $3,
    hour=> $4,
    minute => $5,
    second => $6,
    time_zone => 'Europe/Madrid',
  );
  my $f = DateTime::Format::W3CDTF->new;
  my $time_w3c = $f->format_datetime($dt);
 
 return $time_w3c;
}  

sub LogReport_call {
  my %args = @_;

  # relleno los parametros por defecto de LogReport_call

  my $debug      = $args{debug}  || 0;
  my $queryID    = $args{queryID} || "";
  my $parameters = $args{parameters} || undef;
	
  # Get the parameters
	
  my $startTime    = $parameters->{startTime};
  my $endTime      = $parameters->{endTime};
  my $includeTests = $parameters->{includeTests};

  my $inb_startTime = _W3C_to_inb ($startTime);
  my $inb_endTime   = _W3C_to_inb ($endTime);

  if ($debug) {	
    print STDERR "LogReport parameters (startTime, endTime, includeTests): $startTime, $endTime, $includeTests.\n";
    print STDERR "INB formatted start and end time, $inb_startTime, $inb_endTime\n";
  }

  # Output definition
	
  my $logEvents_href = {};
  my $moby_exceptions = [];

  my $log_file = $gmaster_home . "/projects/moby_logs/moby_services_statistics.log";
  # $log_file = "/home/ug/arnau/data/moby_services_statistics.log";
  
  if (! -f $log_file) {
  
    # raise an exception
    
    my $note = "Internal System Error. No log file was found!!\n";
    print STDERR "$note\n";
    my $code = 701;
    my $moby_exception = INB::Exceptions::MobyException->new (
                                                              code       => $code,
                                                              type       => 'error',
                                                              queryID    => $queryID,
                                                              message    => "$note",
                                                             );
    push (@$moby_exceptions, $moby_exception);
    
    return ($logEvents_href, 0, $moby_exceptions);
  }
  
  # always the case for us
  my $_number_of_cpus = 1;
  my $number_of_events = 0;
  
  open LOG, "$log_file";
  while (<LOG>) {
    my $line = $_;
    chomp $line;
    
    #2007/01/22 18:05:01> 116948550111510    START = 20070122180501
    #2007/01/22 18:05:01> 116948550111510    SERVICE = filterSequencesAndQualityDataByLength
    #2007/01/22 18:05:01> 116948550111510    URI = genome.imim.es
    #2007/01/22 18:05:01> 116948550111510    IP = 137.82.67.190
    #2007/01/22 18:05:01> 116948550111510    RESULT = filterSequencesAndQualityDataByLength terminated successfully
    #2007/01/22 18:05:01> 116948550111510    STATUS = 0
    #2007/01/22 18:05:01> 116948550111510    END = 20070122180501
    #2007/01/22 18:05:01> 116948550111510    TOTALEXECUTIONTIME =  0 wallclock secs ( 0.05 usr +  0.01 sys =  0.06 CPU)
    
    if (! ($line =~ /#/)) {
      if ($line =~/^.+> \d+\t\w+ = .+/) {
        $line =~ /^.+> (\d+)\t(\w+) = (.+)/;
        
        my $id    = $1;
        my $key   = $2;
        my $value = $3;
        
        if (!exists ($logEvents_href->{$id})) {
        
          # New entry
        
          # key must be 'START' !
          
          if ($key ne "START") {
          
            if ($debug) {
          
              # fine in the case that start was not within the input range
              # otherwise it is an error
          
              # print STDERR "Error, expected the key to be START as it is a new entry\n";
              # print STDERR "line, $line\n";
            }
          }
          else {
            
            if (($value >= $inb_startTime) && ($value <= $inb_endTime)) {
          
              # Log it
              
              my $formatted_startTime = _inb_to_W3C ($value);
          
              my $logEvent = {
                              ID => $id,
                              $key => $formatted_startTime,
                              NUMBER_CPUs => $_number_of_cpus
                             };
              $logEvents_href->{$id} = $logEvent;
              $number_of_events ++;
            }
          }
        }
        else {
        
          # Updating the hash data
        
          my $logEvent = $logEvents_href->{$id};
          
          if (defined $logEvent) {
          
            if ($key =~ /^SERVICE|IP|STATUS|END$/) {
            
              if ($debug) {
                # print STDERR "adding $key\n";
              }
            
              if ($key eq "END") {
                $value = _inb_to_W3C ($value);
              }
            
              $logEvent->{$key} = $value;
              if ($key eq "IP") {
                if ($value =~ /137.82.67.190/) {
                  $logEvent->{IS_TEST} = "true";
                }
                else {
                  $logEvent->{IS_TEST} = "false";
                }
              }
            }
          
            # no need to put it back, as it is by reference !
            # $logEvents_href->{$id} = $logEvent;
          
          }
          
        }
        
      }
      else {
        print STDERR "Error, parsing log file, line, $line\n";
        # for now, exit
        exit 1;
      }
    } # End line not matching '#'
  } # End processing current line
  close LOG;
	
  return ($logEvents_href, $number_of_events, $moby_exceptions);
  
}

1;

__END__
