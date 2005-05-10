#!/usr/local/bin/perl -w

# Issue warnings about suspicious programming.
use warnings 'all';

# Must declare and initialize all variables
use strict;

# Use the code module that contains MOBY Service.
use MOBY::Client::Central;
use MOBY::CommonSubs;

# be prepare for command-line options/arguments
use Getopt::Std;

sub help {
    return <<"END_HELP";
Description: Deregister a service in Moby Central
  Usage:
    
    deregisterService.pl [-h] -x {Moby Central} -s {Service Name}
    -h help
	-x MOBY Central: Chirimoyo, Xistral, Inab or BioMoby
	<1> or Chirimoyo
	<2> or Xistral
	<3> or Inab
	<4> or BioMoby
	-s Service Name
	
	Examples using some combinations:
	perl deregisterService.pl -x 1 -s runGeneIDGFF

END_HELP

}

BEGIN {
	
    # Determines the options with values from program
    use vars qw/$opt_h $opt_x $opt_s/;
    
    # these are switches taking an argument (a value)
    my $switches = 'hxs';
    
    # Get the switches
    getopt($switches);
    
    # If the user does not write nothing, skip to help
    if (defined($opt_h) || !defined($opt_x) || !defined($opt_s)){
	print STDERR help;
	exit 0;
    }
    
}

# URI
$::authURI = 'genome.imim.es';

# Service Name
my $serviceName = $opt_s;

# MOBY Central configuration

# Default registry server is Chirimoyo in Malaga

my $MOBY_URI    = $ENV{MOBY_URI}='http://chirimoyo.ac.uma.es/MOBY/Central';
my $MOBY_SERVER = $ENV{MOBY_SERVER}='http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';

if (defined($opt_x)) {
    # Delete spaces
    $opt_x =~ s/\s//g;
    
    # Assign the MOBY Server and MOBY URI
    if (($opt_x == 1) || ($opt_x eq 'Chirimoyo')) {
	
	$MOBY_URI    = $ENV{MOBY_URI}    = 'http://chirimoyo.ac.uma.es/MOBY/Central';
	$MOBY_SERVER = $ENV{MOBY_SERVER} = 'http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';
	
    }
    elsif (($opt_x == 2) || ($opt_x eq 'Xistral')) {
	
	$MOBY_URI    = $ENV{MOBY_URI}    = 'http://xistral/MOBY/Central';
	$MOBY_SERVER = $ENV{MOBY_SERVER} = 'http://xistral/cgi-bin/MOBY-Central.pl';
	
    }
    elsif (($opt_x == 3) || ($opt_x eq 'Inab')) {
	
	$MOBY_URI    = $ENV{MOBY_URI}    = 'http://www.inab.org/MOBY/Central';
	$MOBY_SERVER = $ENV{MOBY_SERVER} = 'http://www.inab.org/cgi-bin/MOBY-Central.pl';
	
    }
    elsif (($opt_x == 4) || ($opt_x eq 'BioMoby')) {
	
	$MOBY_URI    = $ENV{MOBY_URI}    = 'http://mobycentral.cbr.nrc.ca/cgi-bin/MOBY05/Central';
	$MOBY_SERVER = $ENV{MOBY_SERVER} = 'http://mobycentral.cbr.nrc.ca/cgi-bin/MOBY05/mobycentral.pl';
	
    }
    else {
	print STDERR help;
	exit 0;
    }
    
}
else {
    print STDERR help;
    exit 0;
}

# Connect to MOBY-Central registries for searching.
my $Central = MOBY::Client::Central->new (
					  Registries => {mobycentral => {URL => => $MOBY_SERVER, URI => $MOBY_URI}}
					  );

# Check if the service exist !!

my($sia,$ro);
($sia,$ro) = $Central->findService(serviceName=>$serviceName,authURI=>$::authURI);

if ((defined $ro) || (not defined $sia)) {
	print STDERR "No service, $serviceName, registred at $::authURI.\n";
	exit 0;
}

# Declare register variable.
my $REG = $Central->deregisterService (
				       serviceName  => "$serviceName",
				       authURI => $::authURI
				       );

# Check if the result has been successfully.
if ($REG->success) {
    
    # The result is valid.
    print "The 'runGeneIDGFF' service has been deregistered in Malaga successfully: ", $REG->success, "\n";
    
} else {
    
    # The result is valid.
    print "The 'runGeneIDGFF' service has failed: ", $REG->message,"\n"; 
    
}

