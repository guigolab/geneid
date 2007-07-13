#!/usr/bin/perl -w
# Header to make the script executable.

# Issue warnings about suspicious programming.
use warnings 'all';

# Must declare and initialize all variables
use strict;

# Use the code module that contains MOBY Service.
use MOBY::Client::Central;
use MOBY::CommonSubs;

# Client and server side SOAP implementation
use SOAP::Lite + 'trace';
 
# be prepare for command-line options/arguments
use Getopt::Std;

sub help {
return <<"END_HELP";
Description: Degister the INB Log Services.
Usage:

deregister-INB-Log-Services.pl
	[-h help]
	-x MOBY Central: inab, moby-dev or BioMoby
	-a authURI for Moby Service. For example, "cnio.es"

END_HELP

}

BEGIN {
	
	# Determines the options with values from program
	use vars qw/$opt_h $opt_x $opt_a/;
	   
	# these are switches taking an argument (a value)
	my $switches = 'hxa';
	   
	# Get the switches
	getopt($switches);
	
	# If the user does not write nothing, skip to help
	if ( defined($opt_h) || !defined($opt_x) || !defined($opt_a) ){
		print help;
		exit 0;
	}
	
}

#############
# Constants #
#############
my $mobycentral;
my $authURI = $opt_a ;

#######################################
# Assign the MOBY URI and MOBY SERVER #
#######################################
if ( defined($opt_x) ) {

	# Delete spaces
	$opt_x =~ s/\s//g;
	
	my(%MCENTRAL)=(

		'inab' => "inab",
		'Inab' => "inab",

		'moby-dev' => "moby-dev",
		'Moby-dev' => "moby-dev",
		'Moby-Dev' => "moby-dev",

		'biomoby' => "biomoby",
		'Biomoby' => "biomoby",
		'bioMoby' => "biomoby",
		'BioMoby' => "biomoby",

	);
		
	# Store the new database name.
	$mobycentral = $MCENTRAL{$opt_x};

	# Store the new database name.
	$mobycentral = $opt_x;


	# Assign the MOBY Server and MOBY URI
	if ($mobycentral eq 'inab') {
	
		# export MOBY_URI
		# export MOBY_SERVER
		$ENV{MOBY_URI}='http://www.inab.org/MOBY/Central';
		$ENV{MOBY_SERVER}='http://www.inab.org/cgi-bin/MOBY-Central.pl';
	
	}elsif ($mobycentral eq 'moby-dev') {

		# export MOBY_URI
		# export MOBY_SERVER
		$ENV{MOBY_URI}='http://moby-dev.inab.org/MOBY/Central';
		$ENV{MOBY_SERVER}='http://moby-dev.inab.org/cgi-bin/MOBY-Central.pl';

	}elsif ($mobycentral eq 'biomoby') {

		# export MOBY_URI
		# export MOBY_SERVER
		$ENV{MOBY_URI}='http://mobycentral.icapture.ubc.ca/MOBY/Central';
		$ENV{MOBY_SERVER}='http://mobycentral.icapture.ubc.ca/cgi-bin/MOBY05/mobycentral.pl';

	}else {
		print help;
		exit 0;
	}
	
}else {
	print help;
	exit 0;
}



# Connect to MOBY-Central registries for searching.
my $Central = MOBY::Client::Central->new();

# The new MOBY Service Instance is registered.
my ($REG) = $Central->deregisterService(
			serviceName  => "getStatisticalLog",
			authURI => $authURI
		);
	
# Check if the result has been registered successfully.
if ($REG->success) {
	# The result is valid.
	print "The 'getStatisticalLog' service has been deregistered in $mobycentral successfully: ", $REG->success, "\n";
} else {
	# The result is valid.
	print "The 'getStatisticalLog' service has failed: ", $REG->message,"\n"; 
}
