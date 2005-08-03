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
    
    deregisterService.pl [-h] -x {Moby Central} -s {Service Name} -a {Authoritative URI}
    -h help
	-x MOBY Central: Chirimoyo, Xistral, Inab or BioMoby
	<1> or Chirimoyo
	<2> or Xistral
	<3> or Inab
	<4> or BioMoby
	-s Service Name
	-a Authoritative URI

	Examples using some combinations:
	perl deregisterService.pl -x 1 -s runGeneIDGFF -a genome.imim.es

END_HELP

}

BEGIN {
	
    # Determines the options with values from program
    use vars qw/$opt_h $opt_x $opt_s $opt_a/;
    
    # these are switches taking an argument (a value)
    my $switches = 'hxsa';
    
    # Get the switches
    getopt($switches);
    
    # If the user does not write nothing, skip to help
    if (defined($opt_h) || !defined($opt_x) || !defined($opt_s) || !defined($opt_a)){
	print STDERR help;
	exit 0;
    }
    
}

# URI
$::authURI = $opt_a;

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
	
	print STDERR "It is not adviced to deregister a service in production!\n";
        print STDERR "Contact Oswaldo Trelles (ots\@ac.uma.es) or Sergio Ramirez (serr\@ac.uma.es) for updating a service\n";
        exit 0;
	
    }
    elsif (($opt_x == 4) || ($opt_x eq 'BioMoby')) {
	
	# Production

	$MOBY_URI    = $ENV{MOBY_URI}    = 'http://mobycentral.icapture.ubc.ca/MOBY/Central';
	$MOBY_SERVER = $ENV{MOBY_SERVER} = 'http://mobycentral.icapture.ubc.ca/cgi-bin/MOBY05/mobycentral.pl';

	# Development

	# $MOBY_URI    = $ENV{MOBY_URI}    = 'http://bioinfo.icapture.ubc.ca/MOBY/Central';
	# $MOBY_SERVER = $ENV{MOBY_SERVER} = 'http://bioinfo.icapture.ubc.ca/cgi-bin/mobycentral/MOBY-Central.pl';
	
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

print STDERR "Connecting to registry server, $MOBY_URI...\n";

# Connect to MOBY-Central registries for searching.
my $Central = MOBY::Client::Central->new (
					  Registries => {mobycentral => {URL => => $MOBY_SERVER, URI => $MOBY_URI}}
					  );

print STDERR "Connection done.\n";

# Check if the service exist !!

my($sia,$ro);
($sia,$ro) = $Central->findService(serviceName=>$serviceName,authURI=>$::authURI);

if ((defined $ro) || (not defined $sia)) {
	print STDERR "No service, $serviceName, registered at $::authURI.\n";
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
    print "The '$serviceName' service has been deregistered successfully: ", $REG->success, "\n";
    
} else {
    
    # The result is valid.
    print "The '$serviceName' service deregistration has failed: ", $REG->message,"\n"; 
    
}

