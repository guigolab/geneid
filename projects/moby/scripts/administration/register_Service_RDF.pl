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

# For debugging
use Data::Dumper;

sub help {
    return <<"END_HELP";
Description: Register a service in Moby Central
  Usage:
    
    registerService.pl [-h] -x {Moby Central} -s {Service Name}
    -h help
	-x MOBY Central: Chirimoyo, Mobydev, Inab or BioMoby
	<1> or Chirimoyo
	<2> or Mobydev
	<3> or Inab
	<4> or BioMoby
	-s Signature URL (where is the RDF document)
	
	Examples using some combinations:
	production version :
	   => perl registerService.pl -x 4 -s http://genome.imim.es/webservices/RDFs/runGeneIDGFF.rdf
	development version:
	   => perl registerService.pl -x 4 -s http://genome.imim.es/webservices/RDFs/devel/runGeneIDGFF.rdf

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

# MOBY Central configuration

# Default registry server is Chirimoyo in Malaga

my $MOBY_URI = $ENV{MOBY_URI}    ='http://chirimoyo.ac.uma.es/MOBY/Central';
my $MOBY_SERVER = $ENV{MOBY_SERVER} ='http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';

# URL

if (defined($opt_x)) {
    # Delete spaces
    $opt_x =~ s/\s//g;
    
    # Assign the MOBY Server and MOBY URI
    if (($opt_x == 1) || ($opt_x eq 'Chirimoyo')) {
	
	$MOBY_URI    = $ENV{MOBY_URI}    = 'http://chirimoyo.ac.uma.es/MOBY/Central';
	$MOBY_SERVER = $ENV{MOBY_SERVER} = 'http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';

    }
    elsif (($opt_x == 2) || ($opt_x eq 'Mobydev')) {
	
	$MOBY_URI    = $ENV{MOBY_URI}    = 'http://moby-dev.inab.org/MOBY/Central';
	$MOBY_SERVER = $ENV{MOBY_SERVER} = 'http://moby-dev.inab.org/cgi-bin/MOBY-Central.pl';

    }
    elsif (($opt_x == 3) || ($opt_x eq 'Inab')) {

	$MOBY_URI    = $ENV{MOBY_URI}    = 'http://www.inab.org/MOBY/Central';
	$MOBY_SERVER = $ENV{MOBY_SERVER} = 'http://www.inab.org/cgi-bin/MOBY-Central.pl';

    }
    elsif (($opt_x == 4) || ($opt_x eq 'BioMoby')) {
	
	$MOBY_URI    = $ENV{MOBY_URI}    = 'http://mobycentral.icapture.ubc.ca/MOBY/Central';
	$MOBY_SERVER = $ENV{MOBY_SERVER} = 'http://mobycentral.icapture.ubc.ca/cgi-bin/MOBY05/mobycentral.pl';
	
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

# URI
$::authURI = 'genome.imim.es';

# Contact e-mail
$::contactEmail = 'akerhornou@imim.es';

# Service Name

my $signatureURL = $opt_s;

# Connect to MOBY-Central registries for searching.
my $Central = MOBY::Client::Central->new (
					  Registries => {mobycentral => {URL => $MOBY_SERVER, URI => $MOBY_URI}}
					 );

print STDERR "Registrying service, defined by document at URL, $signatureURL, in registry, $MOBY_URI...\n";

my ($REG) = $Central->registerService(
				      signatureURL => $signatureURL,
				     );

# Check if the result has been registered successfully.
if ($REG->success) {
    # The result is valid.
    print STDERR "The '$signatureURL' service has been registered in successfully: ", $REG->success, "\n";
    
} else {
    # The result is valid.
    print STDERR "The '$signatureURL' service has failed: ", $REG->message,"\n"; 
}
