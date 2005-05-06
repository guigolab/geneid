# Name of the package
package INB::Utils::TestAsync;

# Perl pragma to restrict unsafe constructs
use strict;

# Issue warnings about suspicious programming.
use warnings;

# Client and server side SOAP implementation
use SOAP::Lite;

# Return name and handle of a temporary file safely
#use File::Temp qw/ tempfile tempdir /;
require File::Temp;
use File::Temp ();

# A set of exportable subroutines that are useful in clients and 
# services to deal with the input/output from MOBY Services
use MOBY::CommonSubs qw(:all);

# Use XML library
use XML::DOM;
use XML::LibXML;

sub runTestAsync($$) {

	# Get current input, dividing itself into user from calling and data.
	my ($caller, $MOBY_message) = @_;
	
	# Declare variable that contains MOBY response.
	my $MOBY_RESPONSE = 'Something';
	
	# Declare ID from MOBY input
	my $inputID;
	
	# Declare namespace variables.
	my $inputNamespace;
	
	# Obtain the array of MOBY Data inputs.
	my (@MOBYDatainputs) = getInputs($MOBY_message);

	# If there are not input, return an empty response.
	return SOAP::Data->type('string' => responseHeader('pdg.cnb.uam.es') . responseFooter ())
	unless (scalar(@MOBYDatainputs));

	# Every MOBY Data is scanned.
	foreach my $MOBYDatainput (@MOBYDatainputs) { 
	
		# Get query ID from current MOBY Data.
		my $queryID = getInputID($MOBYDatainput);
		
		# Get the Simple|Collection|Secondary articles lists making up this MOBY Data.
		my @articlesList = getArticles($MOBYDatainput);
		
		# Article List is scanned getting corresponding input.
		foreach my $article (@articlesList) {
		
			# The article name and input article are obtained.
			my ($articleName, $inputArticle) = @{$article};
		
			# Check if type if input is SIMPLE
			if (isSimpleArticle($inputArticle)) {
		
				# Create MOBY Response joining sequence response.
				$MOBY_RESPONSE .= simpleResponse($SEQUENCE_RESPONSE, "", $queryID);
			
			# Check if type if input is SECONDARY PARAMETER
			} elsif (isSecondaryArticle($inputArticle)) {
			
			
			# Check if type if input is COLLECTION
			} elsif (isCollectionArticle($inputArticle)) {
			
				# Create MOBY Response joining sequence response.
				$MOBY_RESPONSE .= collectionResponse($Result, "", $queryID);
			}
			
			
			
		}
		
	} # End of query loop
	
	# Return the ERROR RESPONSE firstly and then the MOBY RESPONSE
	return SOAP::Data->type('string' => (responseHeader('pdg.cnb.uam.es') . $MOBY_RESPONSE . responseFooter()));
	
} # End runFunCUT

1;