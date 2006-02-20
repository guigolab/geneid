#!/usr/local/bin/perl -w

# Issue warnings about suspicious programming.
use warnings 'all';

# Must declare and initialize all variables
use strict;

# Use the code module that contains MOBY Service.
use MOBY::Client::Central;

# be prepare for command-line options/arguments
use Getopt::Std;

sub help {
return <<"END_HELP";
Description: Retrieve an object with its relationships from Moby Central
Usage:

retrieveObject.pl [-h] -x {Moby Central} -o {Object Name}
	-h help
	-x MOBY Central: Chirimoyo, Mobydev, Inab or BioMoby
		<1> or Chirimoyo
		<2> or mobydev
		<3> or Inab
		<4> or BioMoby
	-o Object Name
	
Examples using some combinations:
	perl retrieveObject.pl -x Inab -o Object

END_HELP

}


BEGIN {
	
	# Determines the options with values from program
	use vars qw/$opt_h $opt_x $opt_o/;
	   
	# these are switches taking an argument (a value)
	my $switches = 'hxo';
	   
	# Get the switches
	getopt($switches);
	
	# If the user does not write nothing, skip to help
	if (defined($opt_h) || !defined($opt_x) || !defined($opt_o)){
		print help;
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
	if (($opt_x == 1) || ($opt_x eq 'Chirimoyo')) {
	
		# export MOBY_URI
		# export MOBY_SERVER
		$ENV{MOBY_URI}='http://chirimoyo.ac.uma.es/MOBY/Central';
		$ENV{MOBY_SERVER}='http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';
	
	}elsif (($opt_x == 2) || ($opt_x eq 'Mobydev')) {
	
		# export MOBY_URI
		# export MOBY_SERVER
		$ENV{MOBY_URI}='http://moby-dev.inab.org/MOBY/Central';
		$ENV{MOBY_SERVER}='http://moby-dev.inab.org/cgi-bin/MOBY-Central.pl';

	}elsif (($opt_x == 3) || ($opt_x eq 'Inab')) {

		# export MOBY_URI
		# export MOBY_SERVER
		$ENV{MOBY_URI}='http://www.inab.org/MOBY/Central';
		$ENV{MOBY_SERVER}='http://www.inab.org/cgi-bin/MOBY-Central.pl';
	
	}elsif (($opt_x == 4) || ($opt_x eq 'BioMoby')) {

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
my $MOBYCentral = MOBY::Client::Central->new();

if (! defined $MOBYCentral) {
	print STDERR "Error, could not instanciate MOBY Central!\n";
	exit 1;
}

# Retrieve the definition (type, description, relations, ...) from the input object.
my $objDefinition = $MOBYCentral->retrieveObjectDefinition($opt_o);

if (defined($objDefinition)) {
print<<EOF;
\t\tObjectType:\t$objDefinition->{objectType}
\t\tDescription:\t$objDefinition->{description}
\t\tcontact E-mail:\t$objDefinition->{contactEmail}
\t\tAuthURI:\t$objDefinition->{authURI}
EOF
	# Object relationships
	my($reltype,$prelarr,$prel);
	print "\t\tRelationships:\n";
	
	while(($reltype,$prelarr)=each(%{$objDefinition->{Relationships}})) {
		print "\t\t* $reltype to:\n";
		foreach $prel (@$prelarr) {
			print "\t\t    ",$prel->[0],' (',$prel->[1],")\n";
		}
	}
	
	#Object XML
	print "\t\tXML info:\n$objDefinition->{XML}\n";
}
