#!/usr/local/bin/perl -w
# Header to make the script executable.

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

retrieveObject.pl [-h] -x {Moby Central} [-b {Object Name}] -o {Output Image file}
	-h help
	-x MOBY Central: Chirimoyo, Moby-dev, Inab or BioMoby
		<1> or Chirimoyo
		<2> or Moby-dev
		<3> or Inab
		<4> or BioMoby
	-b Object Name. If you use this command option the ontology diagram is given from object name. Otherwise the whole ontology is represented.
	-o Output Image File (PNG, JPG, etc) that will be represents ontology diagram. This file is obtained by Dot program (filter for drawing directed graphs) -by default PNG file-

NOTE: We have to be installed dot program.

Examples using some combinations:
	perl retrieveObjectOntology.pl -x Inab -b VirtualSequence -o VirtualSequence_Ontology.png
	perl retrieveObjectOntology.pl -x Inab -o INB_Ontology.png

END_HELP

}


BEGIN {
	
	# Determines the options with values from program
	use vars qw/$opt_h $opt_x $opt_b $opt_o/;
	   
	# these are switches taking an argument (a value)
	my $switches = 'hxbo';
	   
	# Get the switches
	getopt($switches);
	
	# If the user does not write nothing, skip to help
	if (defined($opt_h) || !defined($opt_x)){
		print help;
		exit 0;
	}
	
}

##############################################
#	GLOBAL VARIABLES
##############################################
$::TMP_DOT_FILE = '/tmp/tmp_Ontology_Diagram.dot';

##############################################
#	PROTOTYPE FUNCTIONS
##############################################
sub assignMobyURIandMobyServer();
sub createObjectRelationShipsDirectly($$);
sub recursivePrintOntology($\%$);
sub printOntology($\%);
sub recursiveCreationDotFile_UpTree($$$);
sub recursiveCreationDotFile_DownTree($$\%);
sub createImageOntologyFileByDotProgram($$$\%);


###############################################################################
# NAME: assignMobyURIandMobyServer
#
# DESCRIPTION: Assigns the moby URI and Moby Server from input parameter
#		giving by user
#
# INPUTS: 	- Input parameter giving by user ($opt_x)
#
# OUTPUT: 	- Moby Central name (Chirimoyo, Inab, Xistral, ...)
#
# AUTHOR: Jose Manuel Rodriguez Carrasco -jmrc@cnb.uam.es- (INB-CNB)
###############################################################################
sub assignMobyURIandMobyServer() {

	my ($mobycentral);
	
	if (defined($opt_x)) {
	
		# Delete spaces
		$opt_x =~ s/\s//g;
		
		my(%MCENTRAL)=(
			'1' => "chirimoyo",
			'Chirimoyo' => "chirimoyo",
			'chirimoyo' => "chirimoyo",
			'2' => "mobydev",
			'Mobydev' => "mobydev",
			'mobydev' => "mobydev",
			'3' => "inab",
			'Inab' => "inab",
			'inab' => "inab",
			'4' => "biomoby",
			'BioMoby' => "biomoby",
			'biomoby' => "biomoby",
		);
			
		# Store the new database name.
		($mobycentral) = $MCENTRAL{$opt_x};
	
		# Assign the MOBY Server and MOBY URI
		if ($mobycentral eq 'chirimoyo') {
		
			# export MOBY_URI
			# export MOBY_SERVER
			$ENV{MOBY_URI}='http://chirimoyo.ac.uma.es/MOBY/Central';
			$ENV{MOBY_SERVER}='http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';
		
		}elsif ($mobycentral eq 'mobydev') {
		
			# export MOBY_URI
			# export MOBY_SERVER
			$ENV{MOBY_URI}='http://moby-dev.inab.org/MOBY/Central';
			$ENV{MOBY_SERVER}='http://moby-dev.inab.org/cgi-bin/MOBY-Central.pl';
	
		}elsif ($mobycentral eq 'inab') {
	
			# export MOBY_URI
			# export MOBY_SERVER
			$ENV{MOBY_URI}='http://www.inab.org/MOBY/Central';
			$ENV{MOBY_SERVER}='http://www.inab.org/cgi-bin/MOBY-Central.pl';
			
		}elsif ($mobycentral eq 'biomoby') {
		    
		    $ENV{MOBY_SERVER}='http://mobycentral.icapture.ubc.ca/cgi-bin/MOBY05/mobycentral.pl';
		    $ENV{MOBY_URI}='http://mobycentral.icapture.ubc.ca/MOBY/Central';
		    
		}else { print help; exit 0;}

	}else { print help; exit 0; }
	
	# Return name of moby central 
	return $mobycentral;

} # End assignMobyURIandMobyServer

###############################################################################
# NAME: createObjectRelationShipsDirectly
#
# DESCRIPTION: Creates a new hash variable which is going to contains:
#		For each object
#		{"Name of object"} -> [Objects that have relationships among them]
#
# INPUTS: 	- Moby Central
#		- List of all registered object types
#
# OUTPUT: 	- Hash variable that contains the whole object ontology of
#		current moby Central
#
# AUTHOR: Jose Manuel Rodriguez Carrasco -jmrc@cnb.uam.es- (INB-CNB)
###############################################################################
sub createObjectRelationShipsDirectly($$) {

	my ($MOBYCentral, $refObjectNameList) = @_;
	my (%ontology);

	# For each object inside of Object List
	while( my($objName, $objDefinition) = each(%{$refObjectNameList})) {

		# Retrieve the definition (type, description, relations, ...) from the input object.
		my ($actualObjDefinition) = $MOBYCentral->retrieveObjectDefinition($objName);

		if (defined($actualObjDefinition)) {
		
			# Create the new hash that will be contains the whole object ontology
			while( my ($reltype, $refRelarr) = each(%{$actualObjDefinition->{Relationships}})) {
				if ($reltype eq 'urn:lsid:biomoby.org:objectrelation:isa') {
				
					my ($objRelation) = ${@{$refRelarr}}[0][0];
					$objRelation =~ s/urn:lsid:biomoby.org:objectclass://g;
					
					# the ISA object has been added before
					if (exists $ontology{$objRelation}) {
						push(@{$ontology{$objRelation}}, $objName);
					} else { # New ISA object
						$ontology{$objRelation} = [$objName];
					}
				}
			}
			
		} else { print "\nERROR: Returning the ontology\n"; exit 1; }
	}
	print "OK\n";

	return %ontology;

} # End createObjectRelationShipsDirectly

###############################################################################
# NAME: recursivePrintOntology
#
# DESCRIPTION: Goes throw recursively the ontology printing every object and 
#		its children
#
# INPUTS: 	- Actual object name
#		- Ontology tree
#		- Number of spaces. It is used to represents more clear the tree
#
# AUTHOR: Jose Manuel Rodriguez Carrasco -jmrc@cnb.uam.es- (INB-CNB)
###############################################################################
sub recursivePrintOntology($\%$){
	my ($objName, $refOntology, $numSpaces) = @_;
	
	if ($objName eq 'Object') {
		print "\t$objName\n";
	} else {
		print "\t"."    " x $numSpaces;
		print "|-->$objName\n";
	}
		
	($numSpaces) ++;
	foreach my $v (@{$refOntology->{$objName}}) {
		if (exists $refOntology->{$v}) {
			recursivePrintOntology($v, %{$refOntology}, $numSpaces);
		} else {
			print "\t"."    " x $numSpaces;
			print "|-->$v\n";
		}
	}

} # End recursivePrintOntology

###############################################################################
# NAME: printOntology
#
# DESCRIPTION: Prints into standard output the ontology tree. The funcion is 
#		called recursively to print each node
#
# INPUTS: 	- Actual object name
#		- Ontology tree
#
# AUTHOR: Jose Manuel Rodriguez Carrasco -jmrc@cnb.uam.es- (INB-CNB)
###############################################################################
sub printOntology($\%){

	my ($objName, $refOntology) = @_;

	print "Printing ontology diagram...";
	print "\n--------------------------------------------------------------------\n\n";

	# Print the ontology tree in standard output
	recursivePrintOntology($objName, %{$refOntology}, 0);

	print "\n\n--------------------------------------------------------------------\n";

} # End printOntology

###############################################################################
# NAME: recursiveCreationDotFile_UpTree
#
# DESCRIPTION: Write recursively into dot file the down tree from giving object name
#
# INPUTS: 	- DOT File
#		- Input object Name
#		- Ontology tree (by reference)
#
# AUTHOR: Jose Manuel Rodriguez Carrasco -jmrc@cnb.uam.es- (INB-CNB)
###############################################################################
sub recursiveCreationDotFile_UpTree($$$) {

	my ($OUTPUTFILE, $MOBYCentral, $objName) = @_;

	unless ($objName eq 'Object') { # Leaf nodes
		# Retrieve the definition (type, description, relations, ...) from the input object.
		my ($actualObjDefinition) = $MOBYCentral->retrieveObjectDefinition($objName);
	
		if (defined($actualObjDefinition)) {
			
			while( my ($reltype, $refRelarr) = each(%{$actualObjDefinition->{Relationships}})) {
				if ($reltype eq 'urn:lsid:biomoby.org:objectrelation:isa') {
					
					my ($objFather) = ${@{$refRelarr}}[0][0];
					$objFather =~ s/urn:lsid:biomoby.org:objectclass://g;
						
					# The ISA object is added into dot file
					print $OUTPUTFILE "\t\"$objFather\" -> \"$objName\"\n";
					print $OUTPUTFILE "\t\"$objName\" [color=black, fillcolor=burlywood, style=filled]\n";
					recursiveCreationDotFile_UpTree($OUTPUTFILE, $MOBYCentral, $objFather);

				}
			}
		}
	} else { # Root node
		print $OUTPUTFILE "\t\"$objName\" [color=black, fillcolor=burlywood, style=filled]\n";
	}

} # End recursiveCreationDotFile_UpTree

###############################################################################
# NAME: recursiveCreationDotFile_DownTree
#
# DESCRIPTION: Write recursively into dot file the up tree from giving object name
#
# INPUTS: 	- DOT File
#		- Input object Name
#		- Ontology tree (by reference)
#
# AUTHOR: Jose Manuel Rodriguez Carrasco -jmrc@cnb.uam.es- (INB-CNB)
###############################################################################
sub recursiveCreationDotFile_DownTree($$\%) {

	my ($OUTPUTFILE, $objName, $refOntology) = @_;

	foreach my $objChild (@{$refOntology->{$objName}}) {
		if (exists $refOntology->{$objChild}) {
			print $OUTPUTFILE "\t\"$objName\" -> \"$objChild\"\n";
			print $OUTPUTFILE "\t\"$objName\" [color=black, fillcolor=burlywood, style=filled]\n";
			recursiveCreationDotFile_DownTree($OUTPUTFILE, $objChild, %{$refOntology});
		} else { # leaf node
			print $OUTPUTFILE "\t\"$objName\" -> \"$objChild\"\n";
			print $OUTPUTFILE "\t\"$objName\" [color=black, fillcolor=burlywood, style=filled]\n";
			print $OUTPUTFILE "\t\"$objChild\" [color=black, fillcolor=burlywood, style=filled]\n";
		}
	}
} # End recursiveCreationDotFile_DownTree

###############################################################################
# NAME: createImageOntologyFileByDotProgram
#
# DESCRIPTION: Creates file that has "dot" format. 
#		"dot" - filter for drawing directed graphs (see manual)
#		It is scanned every node from ontology and creates 
#		each relation among objects. Therefore, the dot file is created 
#		in two parts. First, is created the relationships between input
#		object and its fathers (up-tree) and then, among its children
#		(down-tree)
#
# INPUTS: 	- Name of Moby Central
#		- Moby Central Instance
#		- Input object Name
#		- Ontology tree (by reference)
#
# AUTHOR: Jose Manuel Rodriguez Carrasco -jmrc@cnb.uam.es- (INB-CNB)
###############################################################################
sub createImageOntologyFileByDotProgram($$$\%){
	my ($mobycentralName, $MOBYCentral, $objName, $refOntology) = @_;
	
	print "Creating dot file...";

	local(*DOTFILE);
	open(DOTFILE, ">$::TMP_DOT_FILE") || die "\nCannot create dot file $::TMP_DOT_FILE: $!\n";
	print DOTFILE "digraph $mobycentralName"."_MobyCentralOntology{\n";
	print DOTFILE "\tbgcolor=grey\n";

	recursiveCreationDotFile_UpTree(\*DOTFILE, $MOBYCentral, $objName);
	recursiveCreationDotFile_DownTree(\*DOTFILE, $objName, %{$refOntology});

	print DOTFILE "}\n";
	close(DOTFILE);
	print "OK\n";

	# Execute dot program to draw ontology graph
	if (defined($opt_o)) {
		if (system("dot -o $opt_o -Tjpg $::TMP_DOT_FILE")) {
			print "\nERROR: There was an error during creation of image file\n";
			exit 1;
		}
	}

	# Delete temporal dot file
	unlink($::TMP_DOT_FILE);

} # End createImageOntologyFileByDotProgram

###############################################################################
# 				MAIN PROGRAM
###############################################################################
sub main() { 

	# Check if dot program exists
	print "Dot Program:\n";
	unless (system("which dot") == 0) {print "\nERROR: Not exist Dot Program\n"; exit 1;}

	# Check if input image file has image extension correctly
	unless (defined($opt_o) && (($opt_o =~ /.png/) || ($opt_o =~ /.jpg/) || ($opt_o =~ /.gif/))) {print "\nERROR: Invalid format of input image file\n"; exit 1;}

	# Assign the moby URI and Moby Server from input values giving by user
	my ($mobyCentralName) = assignMobyURIandMobyServer();

	# Connect to MOBY-Central registries for searching.
	print "Connecting to Moby Central...";
	my ($MOBYCentral) = MOBY::Client::Central->new();

	# Check if connection to Moby Central was wrong
	unless(defined($MOBYCentral)) {print "\nERROR: Not available the connection to Moby Central\n"; exit 1;}
	else { print "OK\n"; }

	# Check if input object name is defined and it is valid
	my ($actualObjectName);
	if (defined($opt_b)) { 
		my ($objectDefinition) = $MOBYCentral->retrieveObjectDefinition($opt_b);
		unless (defined(scalar(each(%{$objectDefinition->{Relationships}})))) {print "\nERROR: The input object name $opt_b is not valid\n"; exit 1;}
		($actualObjectName) = $opt_b;
	} else {
		($actualObjectName) = 'Object';
	}

	# Get the list of all registered object types
	print "Creating ontology diagram...\n";
	my ($objectNameList) = $MOBYCentral->retrieveObjectNames();

	# Create a new hash variable which is going to contains all relationships
	# between object and all his children
	my (%ontology);
	if (defined($objectNameList)) { 
		(%ontology) = createObjectRelationShipsDirectly($MOBYCentral, $objectNameList);
	} else { print "\nERROR: Not defined $objectNameList\n"; exit 1; }
		
	# Print the ontology tree in standard output
	printOntology($actualObjectName, %ontology);

	# Create output image file by means of Dot program
	createImageOntologyFileByDotProgram($mobyCentralName, $MOBYCentral, $actualObjectName, %ontology);

}

main();
