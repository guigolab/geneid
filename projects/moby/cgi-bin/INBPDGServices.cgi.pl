#!/usr/local/bin/perl -w

# Must declare and initialize all variables
use strict;

# Issue warnings about suspicious programming.
use warnings 'all';

# Use the SOAP::Transport::HTTP module - this is a CGI SOAP server
use SOAP::Transport::HTTP;

# Add the path of your Services subroutines into your perl lib path
use FindBin;
use lib "$FindBin::Bin"; # the location of your Service modules

# Use the code module that contains our service subroutines.
use Services::SwissProtServices;
use Services::ResidueSumService;
use Services::GenBankServices;
use Services::EMBLServices;
use Services::XNUServices;
#use Services::PDGServices;
use Services::FunCUTServices;


#$ENV{MOBY_SERVER}='http://www.inab.org/cgi-bin/MOBY-Central.pl';
#$ENV{MOBY_URI}='http://www.inab.org/MOBY/Central';
#$ENV{MOBY_PROXY}='No Proxy Server';


# Create an instance of a SOAP CGI server
my $Server= new SOAP::Transport::HTTP::CGI || die "Can't get SOAP: $!\n";

# Make a mapping between SOAP Action headers and code modules calling to next methods.
$Server->dispatch_with({
	'http://biomoby.org/#getAASeqfromSwissProt' => 'Services::SwissProtServices',
	'http://biomoby.org/#getDescriptionfromSwissProt' => 'Services::SwissProtServices',
	'http://biomoby.org/#getEntryfromSwissProt' => 'Services::SwissProtServices',
	'http://biomoby.org/#getKeywordfromSwissProt' => 'Services::SwissProtServices',
	'http://biomoby.org/#getNbResfromSwissProt' => 'Services::ResidueSumService',
	'http://biomoby.org/#getGenericSeqfromGenBank' => 'Services::GenBankServices',
	'http://biomoby.org/#getNucleotideSeqfromEMBL' => 'Services::EMBLServices',
	'http://biomoby.org/#runXNU' => 'Services::XNUServices',
#	'http://biomoby.org/#runFunCUT' => 'Services::PDGServices',
	'http://biomoby.org/#runFunCUT' => 'Services::FunCUTServices',
#	'http://biomoby.org/#runISS' => 'Services::PDGServices',
	'http://biomoby.org/#runISS' => 'Services::FunCUTServices',
	'http://biomoby.org/#runISSComplete' => 'Services::FunCUTServices',
	'http://biomoby.org/#parserISS_Output_Into_Text_Xml' => 'Services::FunCUTServices',
	'http://biomoby.org/#parserISS_Output_Into_NCBI_Blast_Text' => 'Services::FunCUTServices',
	'http://biomoby.org/#parserISS_Output_Into_NCut_Input' => 'Services::FunCUTServices',
	'http://biomoby.org/#runNCut' => 'Services::FunCUTServices',
	'http://biomoby.org/#runOFunCUT' => 'Services::FunCUTServices',
	'http://biomoby.org/#fromFunCUTtoGFF' => 'Services::FunCUTServices',
	'http://biomoby.org/#getAccNumberfromID' => 'Services::FunCUTServices'
});

# The handle method starts the script listening.
$Server->handle;
