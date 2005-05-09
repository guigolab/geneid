#!/usr/local/bin/perl -w

# Issue warnings about suspicious programming.
use warnings 'all';

# Must declare and initialize all variables
use strict;

# Use the code module that contains MOBY Service.
use MOBY::Client::Central;
use MOBY::CommonSubs;

# URI
$::authURI = 'genome.imim.es';
# Service Name
my $serviceName = shift || die "Specify a service name\nderegisterService.pl \"ServiceName\"\n";

# MOBY Central configuration

my $MOBY_URI = $ENV{MOBY_URI}='http://chirimoyo.ac.uma.es/MOBY/Central';
my $MOBY_URL = $ENV{MOBY_SERVER}='http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';

# Connect to MOBY-Central registries for searching.
my $Central = MOBY::Client::Central->new (
					  Registries => {mobycentral => {URL => => $MOBY_URL, URI => $MOBY_URI}}
					 );

# Check if the service exist !!

# ...

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

