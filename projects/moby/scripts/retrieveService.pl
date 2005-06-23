#!/usr/local/bin/perl -w

#edu mola mil
# Issue warnings about suspicious programming.
use warnings 'all';

# Must declare and initialize all variables
use strict;

# Use the code module that contains MOBY Service.
use MOBY::Client::Central;

# be prepare for command-line options/arguments
use Getopt::Std;

use Data::Dumper;

sub help {
return <<"END_HELP";
Description: Get the WSDL definition of the service from Moby Central
Usage:

retrieveService.pl [-h] -x {Moby Central} -s {Service Name} -w {Host Server}
	-h help
	-x MOBY Central: Chirimoyo, Xistral, Inab or BioMoby
		<1> or Chirimoyo
		<2> or Xistral
		<3> or Inab
		<4> or BioMoby
	-s Service Name
	-w Host Server
	
Examples using some combinations:
	perl retrieveService.pl -x Xistral -s runGeneIDGFF -w genome.imim.es

END_HELP

}


BEGIN {
	
	# Determines the options with values from program
	use vars qw/$opt_h $opt_x $opt_s $opt_w/;
	   
	# these are switches taking an argument (a value)
	my $switches = 'hxsw';
	   
	# Get the switches
	getopt($switches);
	
	# If the user does not write nothing, skip to help
	if (defined($opt_h) || !defined($opt_x) || !defined($opt_s) || !defined($opt_w)){
		print help;
		exit 0;
	}
	
}

#########################################################################
#	CONSTANTS
#########################################################################
# URI
$::authURI = 'pdg.cnb.uam.es';
$::authURI = 'genome.imim.es';

if (defined($opt_w)) {
	$::authURI = $opt_w;
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
		$ENV{MOBY_URI}='http://mobycentral.cbr.nrc.ca/cgi-bin/MOBY05/Central';
		$ENV{MOBY_SERVER}='http://mobycentral.cbr.nrc.ca/cgi-bin/MOBY05/mobycentral.pl';

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

my($sia,$si,$ro);
($sia,$ro) = $MOBYCentral->findService(serviceName=>$opt_s,authURI=>$::authURI);

if ((defined $ro) || (not defined $sia) || (not defined $sia->[0])) {
	print STDERR "Can't get information about the given service !!\n";
	if (defined $ro) {
	    print STDERR $ro->message . "\n";
	}
	exit 1;
}

# print STDERR "Dumping sia, " . Dumper ($sia) . "\n";

$si=$sia->[0];
print "\t* Name: ",$si->name,"\n";
print "\t* Contact e-mail: ",$si->contactEmail,"\n";
print "\t* Provider: ",$si->authority,"\n";
print "\t* URL: ", $si->URL,"\n";
print "\t* Type: ",$si->type,"\n";
print "\t* Category: ",$si->category,"\n";
print "\t* Description: ",$si->description,"\n";
print "\t* Input:\n";
foreach my $article (@{$si->input}) {
	if(ref($article) eq 'MOBY::Client::SimpleArticle') {
		print "\t\tarticleName: ",$article->articleName," objectType: ",$article->objectType,"\n";
		my @namespaces = @{$article->namespaces};
		print "\t\tnamespaces: ",join (', ', @namespaces),"\n";
	} else {
		print "\t\tarticleName: ",$article->articleName," collection of=>\n";
		foreach my $sarticle (@{$article->Simples()}) {
			print "\t\t\tarticleName: ",$sarticle->articleName," objectType: ",$sarticle->objectType,"\n";
			my @namespaces = @{$article->namespaces};
			print "\t\tnamespaces: ",join (', ', @namespaces),"\n";
		}
	}
}
print "\t* Output:\n";
foreach my $article (@{$si->output}) {
	if(ref($article) eq 'MOBY::Client::SimpleArticle') {
		print "\t\tarticleName: ",$article->articleName," objectType: ",$article->objectType,"\n";
	} else {
		print "\t\tarticleName: ",$article->articleName," collection of=>\n";
		foreach my $sarticle (@{$article->Simples()}) {
		print "\t\t\tarticleName: ",$sarticle->articleName," objectType: ",$sarticle->objectType,"\n";
		}
	}
}

my @secondaryArticles = @{$si->secondary};
if (@secondaryArticles > 0) {
    foreach my $sa (@secondaryArticles) {
	print "Secondary Article, \"" . $sa->articleName . "\":\n\t\t\t" . $sa->XML . "\n";
    }
}
