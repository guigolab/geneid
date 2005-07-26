#!/usr/local/bin/perl -w

# Issue warnings about suspicious programming.
use warnings 'all';

# Must declare and initialize all variables
use strict;

# Use the code module that contains MOBY Service.
use MOBY::Client::Central;
use MOBY::CommonSubs;

# MOBY Central configuration

my $MOBY_URI = $ENV{MOBY_URI}='http://chirimoyo.ac.uma.es/MOBY/Central';
my $MOBY_URL = $ENV{MOBY_SERVER}='http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';

# URI
$::authURI = 'genome.imim.es';

# Contac e-mail
$::contactEmail = 'akerhornou@imim.es';

# URL
$::URL = 'http://genome.imim.es/cgi-bin/moby/devel/MobyServices.cgi';

# Connect to MOBY-Central registries for searching.
my $Central = MOBY::Client::Central->new (
					  Registries => {mobycentral => {URL => $MOBY_URL, URI => $MOBY_URI}}
					 );

# Check that the service doesn't already exist !

# ...

# Declare register variable.
my ($REG) = $Central->registerService(
				serviceName  => "runSGP2GFF",
				serviceType  => "Analysis",
				authURI      => $::authURI,
				contactEmail => $::contactEmail ,
				description  => "Gene prediction software that works by comparing anonymous genomic sequences from two different species. It combines tblastx, a sequence similarity search program, with geneid, an \"ab initio\" gene prediction program",
				category     => "moby",
				URL          => $::URL,
		input		=> [
				['sequences', ["NucleotideSequence" => []]],
				['tblastx', ["BLAST-Text" => []]],
				],
		output		=> [
				    ['geneid_predictions', ["GFF" => []]],
				    ],
		);
		
		# Check if the result has been registered successfully.
		if ($REG->success) {
			
			# The result is valid.
			print "The 'runSGP2GFF' service has been registered in successfully: ", $REG->success, "\n";
		
		} else {
			
			# The result is valid.
			print "The 'runSGP2GFF' service has failed: ", $REG->message,"\n"; 
			
		}

