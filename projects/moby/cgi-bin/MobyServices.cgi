#!/usr/local/bin/perl -w
#
# MobyServices dispatcher CGI script
#
# Initially written by Roman Roset Mayals, rroset@lsi.upc.es
# Cared by Arnaud Kerhornou, akerhornou@imim.es
# For copyright and disclaimer see below.
# 

use strict;
use EnvServices;
use FindBin qw($Bin);
use lib "$Bin";
use SOAP::Transport::HTTP;
use POSIX qw(setsid);

###############################################################################
BEGIN {
    load_environment();
}

# Statistics and logs managment
use Benchmark;
use Log::Log4perl qw(get_logger :levels);
use Log::Dispatch::File;

###############################################################################

# To know where are installed INB libraries, See EnvServices package
use INB::GRIB::Services::GeneIDServices;
use INB::GRIB::Services::SGP2Services;
use INB::GRIB::Services::GOstatServices;
use INB::GRIB::Services::UtilsServices;
use INB::GRIB::Services::PromoterExtractionServices;
use INB::GRIB::Services::MatScanServices;
use INB::GRIB::Services::MetaAlignmentServices;
use INB::GRIB::Services::MemeServices;
use INB::GRIB::Services::ParsingServices;
use INB::GRIB::Services::ConversionServices;
use INB::GRIB::Services::MaskingServices;
use INB::GRIB::Services::VectorScreeningServices;
use INB::GRIB::Services::AssemblyServices;
use INB::GRIB::Services::BaseCallingServices;

use INB::GRIB::Services::OldServices;
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
					   LocalPort => $port,
					   host      => 'localhost') or die "Can't get SOAP: $!\n";
} else {
    $x = new SOAP::Transport::HTTP::CGI || die "Can't get SOAP: $!\n";
}
# Compression handling
# How does that work ????????
# $x->options({compress_threshold => 1});

##############################################################################
# Stats reporting into a file

###############################################################################
#
# Moby Logger initialisation
#
###############################################################################

my $logfile = $ENV{STATS_FILE};

# No longer using easy_init, because it was duplicating logs in error_log file !!!!
# Log::Log4perl->easy_init($INFO);
my $appender = Log::Log4perl::Appender->new(
					    "Log::Dispatch::File",
					    filename => "$logfile",
					    mode     => "append",
					    );

my $conf = qq(
    log4perl.logger                    = INFO, FileApp
    log4perl.appender.FileApp          = Log::Log4perl::Appender::File
    log4perl.appender.FileApp.filename = $logfile
    log4perl.appender.FileApp.layout   = PatternLayout
    log4perl.appender.FileApp.layout.ConversionPattern = %d> %m%n
    );

# Initialize logging behaviour
Log::Log4perl->init( \$conf );

my $moby_logger = get_logger ("MobyServices");
# No need to append any longer as it is already done !!!
# $moby_logger->add_appender ($appender);

###############################################################################

my $starttime_benchmark = Benchmark->new ();
my $starttime;
{
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime;
  # month = [0..11] so increment it
  $mon++;
  $year += 1900;	
  $starttime = sprintf "%s%2.2d%2.2d%2.2d%2.2d%2.2d", $year, $mon, $mday, $hour, $min, $sec;
}
my $serviceName = "";
my $URI         = "genome.imim.es";
my $IP_address  = $ENV{REMOTE_ADDR};
my $remote_host = $ENV{REMOTE_HOST};
my $unique_id   = int (rand (9)) . int (rand (9));
my $request_id  = $starttime . "_" . $unique_id;

# Disable for the time being until i find out how to propagate the request id...

# $moby_logger->info ("$request_id: User request from remote host, $remote_host ($IP_address)");
# $moby_logger->info ("$request_id: Started at, $starttime");

$moby_logger->info ("User request from remote host, $remote_host ($IP_address)");
$moby_logger->info ("Started at, $starttime");

# Get the service name

$x->on_action(sub {
    my $action = shift;
    $action =~ /^([^#]+)#(\w+)/;
    # die "SOAPAction shall match 'uri#method'\n" if $action ne join '#', @_;
    $serviceName = $2;

    # $moby_logger->info ("$request_id: Executing $serviceName service hosted by service provider authority, $URI");
    $moby_logger->info ("Executing $serviceName service hosted by service provider authority, $URI");
});

$x->dispatch_with({
    'http://biomoby.org/#runGeneID'    => 'INB::GRIB::Services::GeneIDServices',
    'http://biomoby.org/#runGeneIDGFF' => 'INB::GRIB::Services::GeneIDServices',
    
    'http://biomoby.org/#runSGP2GFF'   => 'INB::GRIB::Services::SGP2Services',
    
    'http://biomoby.org/#runGOstat'    => 'INB::GRIB::Services::GOstatServices',

    'http://biomoby.org/#translateGeneIDGFFPredictions' => 'INB::GRIB::Services::UtilsServices',
    
    'http://biomoby.org/#getUpstreamSeqfromEnsembl'     => 'INB::GRIB::Services::PromoterExtractionServices',
    
    'http://biomoby.org/#runMatScanGFF'                 => 'INB::GRIB::Services::MatScanServices',
    'http://biomoby.org/#runMatScanGFFCollection'       => 'INB::GRIB::Services::MatScanServices',
    'http://biomoby.org/#runMatScanGFFCollectionVsInputMatrices' => 'INB::GRIB::Services::MatScanServices',

    'http://biomoby.org/#runMetaAlignment'              => 'INB::GRIB::Services::MetaAlignmentServices',
    'http://biomoby.org/#runMetaAlignmentGFF'           => 'INB::GRIB::Services::MetaAlignmentServices',
    'http://biomoby.org/#runMultiMetaAlignment'         => 'INB::GRIB::Services::MetaAlignmentServices',
    'http://biomoby.org/#runMultiMetaAlignmentGFF'      => 'INB::GRIB::Services::MetaAlignmentServices',
    'http://biomoby.org/#runMultiPairwiseMetaAlignment'          => 'INB::GRIB::Services::MetaAlignmentServices',
    'http://biomoby.org/#runMultiPairwiseMetaAlignmentGFF'       => 'INB::GRIB::Services::MetaAlignmentServices',

    # 'http://biomoby.org/#fromGenericSequencetoFASTA'             => 'INB::GRIB::Services::UtilsServices',
    # 'http://biomoby.org/#fromGenericSequenceCollectiontoFASTA'   => 'INB::GRIB::Services::UtilsServices',
    # 'http://biomoby.org/#fromFASTAtoDNASequenceCollection'       => 'INB::GRIB::Services::UtilsServices',
    # 'http://biomoby.org/#fromFASTAtoAminoAcidSequenceCollection' => 'INB::GRIB::Services::UtilsServices',
    # 'http://biomoby.org/#fromFASTAtoGenericSequenceCollection'   => 'INB::GRIB::Services::UtilsServices',
    
    'http://biomoby.org/#fromGenericSequenceToFASTA'             => 'INB::GRIB::Services::ConversionServices',
    'http://biomoby.org/#fromGenericSequenceCollectionToFASTA'   => 'INB::GRIB::Services::ConversionServices',
    'http://biomoby.org/#fromFASTAToDNASequence'                 => 'INB::GRIB::Services::ConversionServices',
    'http://biomoby.org/#fromFASTAToDNASequenceCollection'       => 'INB::GRIB::Services::ConversionServices',
    'http://biomoby.org/#fromFASTAToAminoAcidSequence'           => 'INB::GRIB::Services::ConversionServices',
    'http://biomoby.org/#fromFASTAToAminoAcidSequenceCollection' => 'INB::GRIB::Services::ConversionServices',
    'http://biomoby.org/#fromFASTAToGenericSequenceCollection'   => 'INB::GRIB::Services::ConversionServices',
    
    'http://biomoby.org/#generateScoreMatrix'                    => 'INB::GRIB::Services::OldServices',
    'http://biomoby.org/#fromMetaAlignmentsToScoreMatrix'        => 'INB::GRIB::Services::ParsingServices',
    'http://biomoby.org/#fromMetaAlignmentsToTextScoreMatrix'    => 'INB::GRIB::Services::ParsingServices',
    
    'http://biomoby.org/#runMemeHTML'                            => 'INB::GRIB::Services::MemeServices',
    'http://biomoby.org/#runMemeText'                            => 'INB::GRIB::Services::MemeServices',
    
    'http://biomoby.org/#parseMotifMatricesFromMEME'             => 'INB::GRIB::Services::ParsingServices',
    
    'http://biomoby.org/#runDust'                                => 'INB::GRIB::Services::MaskingServices',
    'http://biomoby.org/#runDustCollection'                      => 'INB::GRIB::Services::MaskingServices',
    'http://biomoby.org/#runRepeatMasker'                        => 'INB::GRIB::Services::MaskingServices',
    'http://biomoby.org/#runRepeatMaskerCollection'              => 'INB::GRIB::Services::MaskingServices',
    
    'http://biomoby.org/#runCrossMatchToScreenVector'            => 'INB::GRIB::Services::VectorScreeningServices',
    'http://biomoby.org/#runCrossMatchToScreenVectorCollection'  => 'INB::GRIB::Services::VectorScreeningServices',
    
    'http://biomoby.org/#runPhrap'                               => 'INB::GRIB::Services::AssemblyServices',
    'http://biomoby.org/#runPhrapWithQualityData'                => 'INB::GRIB::Services::Assemblyservices',
    
    'http://biomoby.org/#runPhred'                               => 'INB::GRIB::Services::BaseCallingServices',
    'http://biomoby.org/#runPhredCollection'                     => 'INB::GRIB::Services::BaseCallingServices',
    
    'http://biomoby.org/#generateScoreMatrix'                    => 'INB::GRIB::Services::OldServices',
    'http://biomoby.org/#parseMotifMatricesfromMEME'             => 'INB::GRIB::Services::OldServices',
    'http://biomoby.org/#runMatScanGFFCollectionVsInputMatrix'   => 'INB::GRIB::Services::OldServices',
});
$x->handle;

# what else, execution status => el codigo del peor error o OK (700) si OK

my $endtime_benchmark = Benchmark->new ();
my $endtime;
{
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime;
  # month = [0..11] so increment it
  $mon++;
  $year += 1900;	
  $endtime = sprintf "%s%2.2d%2.2d%2.2d%2.2d%2.2d", $year, $mon, $mday, $hour, $min, $sec;
}

# $moby_logger->info ("$request_id: Ending at, $endtime");
# $moby_logger->info ("$request_id: Total execution time: ", timestr (timediff ($endtime_benchmark, $starttime_benchmark)));
$moby_logger->info ("Ending at, $endtime");
$moby_logger->info ("Total execution time: ", timestr (timediff ($endtime_benchmark, $starttime_benchmark)));
$moby_logger->info ("#");
