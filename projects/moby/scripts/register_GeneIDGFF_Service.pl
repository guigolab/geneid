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
$::URL = 'http://cel.imim.es/cgi-bin/moby/MobyServices.cgi';

# Connect to MOBY-Central registries for searching.
my $Central = MOBY::Client::Central->new (
					  Registries => {mobycentral => {URL => $MOBY_URL, URI => $MOBY_URI}}
					 );

# Check that the service doesn't already exist !

# ...

# Declare register variable.
my ($REG) = $Central->registerService(
				serviceName  => "runGeneIDGFF",
				serviceType  => "Analysis",
				authURI      => $::authURI,
				contactEmail => $::contactEmail ,
				description  => "Ab initio gene prediction tool - Return the output predictions in GFF format.",
				category     => "moby",
				URL          => $::URL,
		input		=> [
				['', ["NucleotideSequence" => []]],
				],
		output		=> [
				['', ["GFF" => []]],
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
						enum => ['Forward', 'Reverse','Both'],
						default => 'Both',
						max => 'MAX',
						min => 'MIN',
					},
					'engine' => {
						datatype => 'String',
						default => 'Normal',
						enum => ['Normal','Nogenes','Genamic'],
						max => 'MAX',
						min => 'MIN',
					},
				}
		);
		
		# Check if the result has been registered successfully.
		if ($REG->success) {
			
			# The result is valid.
			print "The 'runGeneIDGFF' service has been registered in successfully: ", $REG->success, "\n";
		
		} else {
			
			# The result is valid.
			print "The 'runGeneIDGFF' service has failed: ", $REG->message,"\n"; 
			
		}

