# $Id: EnvServices.pm,v 1.1 2005-05-06 10:53:53 arnau Exp $
#
# INBPerl module for INB::UPC::NCBI_BLAST::MobyParser
#
# Cared for by Roman Roset Mayals, rroset@lsi.upc.es
# For copyright and disclaimer see below.
#

# POD documentation - main docs before the code

package EnvServices;

use strict;
use warnings;
use Env; 

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [ qw() ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
        &load_environment 
);

our $VERSION = '1.00';

###############################################################################
sub load_environment {

	# Load the Global Environment from a file

	# This file tells where are the INB libraries (Development or Production location)

	if (-f "config.pl") {
		require "config.pl";
	}
	else {
		print STDERR "Error - can't find configuration file, config.pl!\n";
		exit 1;
	}
}
###############################################################################

1;
