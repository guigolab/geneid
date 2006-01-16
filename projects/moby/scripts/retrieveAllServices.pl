#!/usr/local/bin/perl -w

#edu mola mil
# Issue warnings about suspicious programming.
use warnings 'all';

# Must declare and initialize all variables
use strict;

# Use the code module that contains MOBY Service.
use MOBY::Client::Central;

# use SOAP::Lite + 'trace';

# be prepare for command-line options/arguments
use Getopt::Std;

use Data::Dumper;

sub help {
return <<"END_HELP";
Description: Get the WSDL definition of the service from Moby Central
Usage:

retrieveService.pl [-h] -x {Moby Central}
	-h help
	-x MOBY Central: Chirimoyo, Xistral, Inab or BioMoby
		<1> or Chirimoyo
		<2> or Xistral
		<3> or Inab
		<4> or BioMoby
	-t Sorted by service type, instead by service provider (by default)
	
Examples using some combinations:
	perl retrieveAllServices.pl -x Chirimoyo -t 1

END_HELP

}


BEGIN {
	
	# Determines the options with values from program
	use vars qw/$opt_h $opt_x $opt_t/;
	
	# these are switches taking an argument (a value)
	my $switches = 'hxt';
	   
	# Get the switches
	getopt($switches);
	
	# If the user does not write nothing, skip to help
	if (defined($opt_h) || !defined($opt_x)){
		print STDERR help;
		exit 0;
	}
	
}

##############################################
#	ASSIGN THE MOBY URI AND MOBY SERVER
##############################################
if (defined($opt_x)) {

	# Delete spaces
	$opt_x =~ s/\s//g;

	# Assign the MOBY Server and MOBY URI
	if (($opt_x =~ /chirimoyo/i) || (($opt_x =~ /\d/) && ($opt_x == 1))) {
	
		# export MOBY_URI
		# export MOBY_SERVER
		$ENV{MOBY_URI}='http://chirimoyo.ac.uma.es/MOBY/Central';
		$ENV{MOBY_SERVER}='http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';
	
	}elsif (($opt_x =~ /xistral/i) || (($opt_x =~ /\d/) && ($opt_x == 2))) {
	
		# export MOBY_URI
		# export MOBY_SERVER
		$ENV{MOBY_URI}='http://xistral/MOBY/Central';
		$ENV{MOBY_SERVER}='http://xistral/cgi-bin/MOBY-Central.pl';

	}elsif (($opt_x =~ /inab/i) || (($opt_x =~ /\d/) && ($opt_x == 3))) {

		# export MOBY_URI
		# export MOBY_SERVER
		$ENV{MOBY_URI}='http://www.inab.org/MOBY/Central';
		$ENV{MOBY_SERVER}='http://www.inab.org/cgi-bin/MOBY-Central.pl';
	
	}elsif (($opt_x =~ /biomoby/i) || (($opt_x =~ /\d/) && ($opt_x == 4))) {
	    
	    # export MOBY_URI
	    # export MOBY_SERVER
	    
	    # Production Canada server
	    
	    $ENV{MOBY_URI}='http://mobycentral.icapture.ubc.ca/MOBY/Central';
	    $ENV{MOBY_SERVER}='http://mobycentral.icapture.ubc.ca/cgi-bin/MOBY05/mobycentral.pl';
	    
	    # Testing Canada server
	    
	    # $ENV{MOBY_URI}='http://bioinfo.icapture.ubc.ca/MOBY/Central';
	    # $ENV{MOBY_SERVER}='http://bioinfo.icapture.ubc.ca/cgi-bin/mobycentral/MOBY-Central.pl';
	    
	}else {
		print STDERR help;
		exit 0;
	}

}else {
	print STDERR help;
	exit 0;
}

my $sorted_by_providers = 1;
(defined ($opt_t) && $opt_t == 1) && ($sorted_by_providers = 0);

# Connect to MOBY-Central registries for searching.
my $MOBYCentral = MOBY::Client::Central->new();

if ($sorted_by_providers) {
    
    # Sorted by service provider

    # Retrieve service providers first
    
    my @service_provider_URIs = $MOBYCentral->retrieveServiceProviders ();
    
    foreach my $uri (@service_provider_URIs) {
	
	print "*******************************************\n";
	print "*\n";
	print "* Service provider: $uri\n";
	print "*\n";
	print "*******************************************\n\n";
	
	my($si_aref,$ro);
	($si_aref,$ro) = $MOBYCentral->findService(authURI=> $uri);
	
	foreach my $si (@$si_aref) {
	    print "\t* Name: ",$si->name,"\n";
	    print "\t* Contact e-mail: ",$si->contactEmail,"\n";
	    print "\t* service type: ",$si->type,"\n";
	    my $desc = $si->description;
	    print "\t* Description: $desc\n";
	}
    }
}
else {
    # Sorted by service type
    
    my %services_by_type;
    
    my($si_aref,$ro);
    ($si_aref,$ro) = $MOBYCentral->findService();
    
    foreach my $si (@$si_aref) {
	my $type = $si->type;
	
	if (exists $services_by_type{$type}) {
	    my $si_tmp = $services_by_type{$type};
	    push (@$si_tmp, $si);
	    $services_by_type{$type} = $si_tmp;
	}
	else {
	    $services_by_type{$type} = [$si];
	}
    }
    
    my @types = keys (%services_by_type);
    
    foreach my $type (@types) {
	
	print "*******************************************\n";
	print "*\n";
	print "* MOBY service type: $type\n";
	print "*\n";
	print "*******************************************\n\n";
	
	foreach my $si (@{$services_by_type{$type}}) {
	    print "\t* Name: ",$si->name,"\n";
	    print "\t* Authority: ",$si->authority,"\n";
	    print "\t* Contact e-mail: ",$si->contactEmail,"\n";
	    # print "\t* Service type: ",$si->type,"\n";
	    my $desc = $si->description;
	    print "\t* Description: $desc\n";
	}
    }
    
}
