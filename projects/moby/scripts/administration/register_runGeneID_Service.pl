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
	-x MOBY Central: Chirimoyo, Xistral, Inab or BioMoby
	<1> or Chirimoyo
	<2> or Xistral
	<3> or Inab
	<4> or BioMoby
	-s Service Name
	
	Examples using some combinations:
	perl registerService.pl -x 1 -s runGeneID

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

# Development by default, but if the registry server is the production one,
# then the URL will be the production server at genome.imim.es
$::URL = 'http://genome.imim.es/cgi-bin/moby/devel/MobyServices.cgi';

if (defined($opt_x)) {
    # Delete spaces
    $opt_x =~ s/\s//g;
    
    # Assign the MOBY Server and MOBY URI
    if (($opt_x == 1) || ($opt_x eq 'Chirimoyo')) {
	
	$MOBY_URI    = $ENV{MOBY_URI}    = 'http://chirimoyo.ac.uma.es/MOBY/Central';
	$MOBY_SERVER = $ENV{MOBY_SERVER} = 'http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';

	$::URL = 'http://genome.imim.es/cgi-bin/moby/devel/MobyServices.cgi';
	
    }
    elsif (($opt_x == 2) || ($opt_x eq 'Xistral')) {
	
	$MOBY_URI    = $ENV{MOBY_URI}    = 'http://xistral/MOBY/Central';
	$MOBY_SERVER = $ENV{MOBY_SERVER} = 'http://xistral/cgi-bin/MOBY-Central.pl';
	
    }
    elsif (($opt_x == 3) || ($opt_x eq 'Inab')) {

	$MOBY_URI    = $ENV{MOBY_URI}    = 'http://www.inab.org/MOBY/Central';
	$MOBY_SERVER = $ENV{MOBY_SERVER} = 'http://www.inab.org/cgi-bin/MOBY-Central.pl';
	
	# Production
	$::URL = 'http://genome.imim.es/cgi-bin/moby/MobyServices.cgi';

    }
    elsif (($opt_x == 4) || ($opt_x eq 'BioMoby')) {

	$MOBY_URI    = $ENV{MOBY_URI}    = 'http://mobycentral.icapture.ubc.ca/MOBY/Central';
	$MOBY_SERVER = $ENV{MOBY_SERVER} = 'http://mobycentral.icapture.ubc.ca/cgi-bin/MOBY05/mobycentral.pl';

	# Production
	$::URL = 'http://genome.imim.es/cgi-bin/moby/MobyServices.cgi';

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

my $serviceName = $opt_s;

# Connect to MOBY-Central registries for searching.
my $Central = MOBY::Client::Central->new (
					  Registries => {mobycentral => {URL => $MOBY_SERVER, URI => $MOBY_URI}}
					 );

# Check that the service doesn't already exist !

my($sia,$ro);
($sia,$ro) = $Central->findService(serviceName=>$serviceName,authURI=>$::authURI);

if ((defined $sia) && (@$sia > 0)) {
	print STDERR "The registration of the service, $serviceName, has failed. A service already exists at $::authURI.\n";
	exit 0;
}

print STDERR "Registrying service, $serviceName, $::URL from this server, $::URL ...\n";

# Declare register variable.
my ($REG) = $Central->registerService(
				      serviceName  => $serviceName,
				      serviceType  => "Analysis",
				      authURI      => $::authURI,
				      contactEmail => $::contactEmail,
				      description  => "Ab initio gene prediction tool.",
				      category     => "moby",
				      URL          => $::URL,
				      input		=> [
							    ['sequences', ["NucleotideSequence" => []]],
							    ],
				      output		=> [
							    ['geneid_predictions', ["text-html" => []]],
							    ],
				      secondary	=> {
					  'profile' => {
					      datatype => 'String',
					      enum => ['Human','Tetraodon','Drosophila','Celegans','Wheat','Arabidopsis','Rice','Plasmodium','Dictyostelium','Aspergillus','Neurospora','Cryptococcus','Coprinus'],
					      default => 'Human',
					      max => 'MAX',
					      min => 'MIN',
					  },
					  'strands' => {
					      datatype => 'String',
					      enum     => ['Forward', 'Reverse','Both'],
					      default  => 'Both',
					      max      => 'MAX',
					      min      => 'MIN',
					  },
					  'engine' => {
					      datatype => 'String',
					      default  => 'Normal',
					      enum     => ['Normal','Exon Mode','Gene Assembly Mode'],
					      max      => 'MAX',
					      min      => 'MIN',
					  },
				}
				      );

# Check if the result has been registered successfully.
if ($REG->success) {
    
    # The result is valid.
    print STDERR "The '$serviceName' service has been registered in successfully: ", $REG->success, "\n";
    
} else {
    
    # The result is valid.
    print STDERR "The '$serviceName' service has failed: ", $REG->message,"\n"; 
    
}

