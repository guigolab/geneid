#!/usr/local/bin/perl -w
#
# INBPerl module 
#
# Initially written by Roman Roset Mayals, rroset@lsi.upc.es
# Cared by Arnaud Kerhornou, akerhornou@imim.es
# For copyright and disclaimer see below.
# 

use strict;
use FindBin qw($Bin);
use lib "$Bin";
use SOAP::Transport::HTTP;
use EnvServices;
use POSIX qw(setsid);

###############################################################################
BEGIN {
	load_environment();
}

###############################################################################

# To know where are installed INB libraries, See EnvServices package
use INB::GRIB::Services::GeneIDServices;
use INB::GRIB::Services::SGP2Services;
use INB::GRIB::Services::GOstatServices;
use INB::GRIB::Services::UtilsServices;
use INB::GRIB::Services::PromoterExtractionServices;
use INB::GRIB::Services::MatScanServices;
use INB::GRIB::Services::MetaAlignmentServices;
###############################################################################

sub daemonize {

	my $port = shift;

	open STDIN, '/dev/null' or die "Can't read /dev/null: $!";
	defined(my $pid = fork) or die "Can't fork: $!";
	exit if $pid;
	POSIX::setsid or die "Can't start a new session: $!";
	print "[$$]Contact to SOAP server at port $port\n";
	open STDERR, '>&STDOUT' or die "Can't dup stdout: $!";
	umask(0);
}
###############################################################################

my $is_daemon = 0;
my $port      = 8081;
my $x;

if ($ARGV[0] and $ARGV[0] =~ /^--daemon$/) {
	$port = $ARGV[1] || 8081;
	daemonize($port);
	$x = new SOAP::Transport::HTTP::Daemon(
		  LocalPort => $port
		, host      => 'localhost') or die "Can't get SOAP: $!\n";
} else {
	$x = new SOAP::Transport::HTTP::CGI || die "Can't get SOAP: $!\n";
}

$x->dispatch_with({
    'http://biomoby.org/#runGeneID'    => 'INB::GRIB::Services::GeneIDServices',
    'http://biomoby.org/#runGeneIDGFF' => 'INB::GRIB::Services::GeneIDServices',
    'http://biomoby.org/#runSGP2GFF'   => 'INB::GRIB::Services::SGP2Services',
    'http://biomoby.org/#runGOstat'    => 'INB::GRIB::Services::GOstatServices',
    'http://biomoby.org/#translateGeneIDGFFPredictions' => 'INB::GRIB::Services::UtilsServices',
    'http://biomoby.org/#getUpstreamSeqfromEnsembl'     => 'INB::GRIB::Services::PromoterExtractionServices',
    'http://biomoby.org/#runMatScanGFF'                 => 'INB::GRIB::Services::MatScanServices',
    'http://biomoby.org/#runMatScanGFFCollection'       => 'INB::GRIB::Services::MatScanServices',
    'http://biomoby.org/#runMetaAlignment'              => 'INB::GRIB::Services::MetaAlignmentServices',
    'http://biomoby.org/#runMetaAlignmentGFF'           => 'INB::GRIB::Services::MetaAlignmentServices',
    'http://biomoby.org/#runMultiMetaAlignment'         => 'INB::GRIB::Services::MetaAlignmentServices',
    'http://biomoby.org/#runMultiMetaAlignmentGFF'      => 'INB::GRIB::Services::MetaAlignmentServices',
    'http://biomoby.org/#fromGenericSequencetoFASTA'    => 'INB::GRIB::Services::UtilsServices',
    'http://biomoby.org/#fromGenericSequenceCollectiontoFASTA' => 'INB::GRIB::Services::UtilsServices',
    'http://biomoby.org/#generateScoreMatrix'           => 'INB::GRIB::Services::UtilsServices',
});
$x->handle;
