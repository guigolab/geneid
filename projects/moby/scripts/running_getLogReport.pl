#!/usr/local/bin/perl -w
# Header to make the script executable.

# Issue warnings about suspicious programming.
use warnings 'all';

# Must declare and initialize all variables
use strict;

# MOBY libraries
use MOBY::Client::Central;
use MOBY::Client::Service;

# SOAP Lite (if you want not to see the SOAP trace delete 'trace' comment)
use SOAP::Lite + 'trace';

# be prepare for command-line options/arguments
use Getopt::Std;

sub help {
return <<"END_HELP";
Description: Register the INB Log Services.
Usage:

Register_INBServices.pl
	[-h help]
	-x MOBY Central: inab, moby-dev or BioMoby
	
	-a AuthURI for Moby Service	(eg, cnio.es)
	-s Start time			(date format: YYYYMMDDhhmmss)
	-e End time			(date format: YYYYMMDDhhmmss)
	-i Include test			(true or false)
	-o Output file result

Example:
	perl running_getLogReport.pl -x moby-dev -a cnio.es -s 2007-03-22T17:00:00Z -e 2007-03-27T18:00:00Z -i false -o /tmp/t.xml

END_HELP

}

BEGIN {
	# Determines the options with values from program
	use vars qw/$opt_h $opt_x $opt_a $opt_s $opt_e $opt_i $opt_o/;
	   
	# these are switches taking an argument (a value)
	my $switches = 'hxaseio';
	   
	# Get the switches
	getopt($switches);
	
	# If the user does not write nothing, skip to help
	if ( defined($opt_h) || !defined($opt_x) || !defined($opt_a)  || !defined($opt_s)  || !defined($opt_e) || !defined($opt_i) || !defined($opt_o) ){
		print help;
		exit 0;
	}
	
}

#############
# Constants #
#############
my $authURI = $opt_a ;
my $start_time = $opt_s ;
my $end_time = $opt_e ;
my $include_tests = $opt_i ;
my $output_file = $opt_o ;

#######################################
# Assign the MOBY URI and MOBY SERVER #
#######################################
my $mobycentral;
my $URL;
my $URI;
if ( defined( $opt_x ) ) {

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
	
		$URL = 'http://www.inab.org/cgi-bin/MOBY-Central.pl';
		$URI = 'http://www.inab.org/MOBY/Central';

	
	}elsif ($mobycentral eq 'moby-dev') {

		$URL = 'http://moby-dev.inab.org/cgi-bin/MOBY-Central.pl';
		$URI = 'http://moby-dev.inab.org/MOBY/Central';


	}elsif ($mobycentral eq 'biomoby') {

		$URL = 'http://mobycentral.icapture.ubc.ca/cgi-bin/MOBY05/mobycentral.pl';
		$URI = 'http://mobycentral.icapture.ubc.ca/MOBY/Central';

	}else {
		print help;
		exit 0;
	}
	
}else {
	print help;
	exit 0;
}


# Get instance against one MOBY Central
my $Central = MOBY::Client::Central->new(
        Registries => {mobycentral => {URL => $URL,URI => $URI}});

# Get WSDL of MOBY Service
my ($ServiceInstances, $RegObject) = $Central->findService(
                serviceName=> 'getStatisticalLog',
                authURI => $authURI
               );
unless ( defined ( $ServiceInstances->[0] ) ) {
	print STDERR "This MOBY service is not register into given MOBYCentral\n";
	exit 1;
}
my $wsdl = $Central->retrieveService( $ServiceInstances->[0] );


# Execute MOBY Service
my $Service = MOBY::Client::Service->new( service => $wsdl );

my $mobySimple_start_time = qq{ <moby:DateTime moby:namespace="" moby:id="" moby:articleName="">$start_time</moby:DateTime> };

my $mobySimple_end_time = qq{ <moby:DateTime moby:namespace="" moby:id="" moby:articleName="">$end_time</moby:DateTime> };


my $response = $Service->execute(XMLinputlist => [
		[
			'end_time', $mobySimple_end_time,
			'start_time', $mobySimple_start_time,
			'includeTests', "<Value>$include_tests</Value>",
		]
]);

# Save result into file
if ( defined( $response ) ) {
	open (OUTPUT_FILE, ">$output_file") or die "Can not create output file: $!\n";
	print OUTPUT_FILE $response;
	close OUTPUT_FILE;
} else {
	print STDERR "\nThere was not answer\n\n";
	exit 1;
}

# RUN PYTHON SCRIPT
