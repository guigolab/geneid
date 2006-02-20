#!/usr/bin/perl -w
# Header to make the script executable.

# Issue warnings about suspicious programming.
use warnings 'all';

# Must declare and initialize all variables
use strict;

# Use the code module that contains MOBY Service.
use MOBY::Client::Central;
use MOBY::CommonSubs;

# URI
$::authURI = 'genome.imim.es';

# Contac e-mail
$::contactEmail = 'akerhornou@imim.es';


# Connect to MOBY-Central registries for searching.
my $Central = MOBY::Client::Central->new();


# The new MOBY Service Instance is registered.
		my ($REG) = $Central->registerObject(
				objectType  => "GFF6",
				description => "the last version",
				authURI      => $::authURI,
				contactEmail => $::contactEmail ,
		Relationships 	=> {
			ISA	=> [
				['GFF', ""]
				],
			}
		);
		
		# Check if the result has been registered successfully.
		if ($REG->success) {
			
			# The result is valid.
			print "The 'GFF6' object has been registered in $opt_x successfully: ", $REG->success, "\n";
		
		} else {
			
			# The result is valid.
			print "The 'GFF6' object has failed: ", $REG->message,"\n"; 
			
		}
		
