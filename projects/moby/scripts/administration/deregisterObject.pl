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
$::authURI = 'pdg.cnb.uam.es';

# Connect to MOBY-Central registries for searching.
my $Central = MOBY::Client::Central->new();


# The next registered MOBY Object is deregistered.
		my ($REG) = $Central->deregisterObject(
					objectType  => "GFF6", 
					authURI => $::authURI
		);
		
		# Check if the result has been successfully.
		if ($REG->success) {
			
			# The result is valid.
			print "The 'SwissProt_Text' object has been deregistered in $opt_x successfully: ", $REG->success, "\n";
		
		} else {
			
			# The result is valid.
			print "The 'SwissProt_Text' object has failed: ", $REG->message,"\n"; 
			
		}
