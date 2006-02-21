# Production Moby GRIB libraries location
use lib "/home/ug/gmaster/projects/moby/prod/lib";

# BioMoby core libraries
use lib "/home/ug/gmaster/projects/moby/biomoby.0.8.2a/Perl";

# Also need SOAP and Bioperl libraries that are already installed in the default perl library path
# /usr/local/Install/perl-5.8.5/lib/site_perl/5.8.5/Bio

# SOAP v0.67
# use lib "/home/ug/gmaster/projects/lib/5.8.5/i686-linux-thread-multi";
# use lib "/home/ug/gmaster/projects/lib/site_perl/5.8.5";

# Log::Log4perl and depedencies - temporary place !
# Installed Globally now
# use lib "/home/ug/gmaster/projects/lib/site_perl/5.8.5";
# use lib "/home/ug/gmaster/projects/lib/site_perl/5.8.5/i686-linux-thread-multi";

# GOstat setup

# GOstat libraries
use lib "/home/ug/gmaster/projects/gostat/lib";
$ENV{GOSTAT} = "/home/ug/gmaster/projects/gostat";

# SGP2 setup
$ENV{SGP2} = "/home/ug/gmaster/sgp2/sgp2_2003/";

# Statistic reporting file name
$ENV{STATS_FILE}="/home/ug/gmaster/projects/moby_logs/moby_services_statistics.log";

1;
