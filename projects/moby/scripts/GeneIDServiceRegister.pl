#!/usr/bin/perl -w
use warnings 'all';		# Issue warnings about suspicious programming
use strict;				# Must declare and initialize all variables
use MOBY::Client::Central;
use MOBY::CommonSubs;

#$ENV{MOBY_URI}='http://mobycentral.cbr.nrc.ca:8080/MOBY/Central';
#$ENV{MOBY_SERVER}='http://mobycentral.cbr.nrc.ca:8080/cgi-bin/MOBY05/mobycentral.pl';

$ENV{MOBY_SERVER}='http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';
$ENV{MOBY_URI}='http://chirimoyo.ac.uma.es/MOBY/Central'; 



#$ENV{MOBY_URI}='http://xistral.cnb.uam.es/MOBY/Central';
#$ENV{MOBY_SERVER}='http://xistral.cnb.uam.es/cgi-bin/MOBY-Central.pl';

my $Central = MOBY::Client::Central->new();
 
############ REGISTER SERVICE #######################################

my ($REG)=$Central->registerService(
		    serviceName  => "GeneID",  
		    serviceType  => "Analysis",  
		    authURI      => "genome.imim.es",      
		    contactEmail => 'rguigo@imim.es',      
		    description =>  "Ab initio gene prediction tool", 
		    category  	=>  "moby",
		    URL    	=>  "http://cel.imim.es/cgi-bin/Services.cgi",
		    input 	=> [
				['', ["NucleotideSequence" => []]],  
                     ],
		    output	=> [
				['', ["text-html" => []]], # Simple
			   ]
			);
print $REG->message,"\n";
