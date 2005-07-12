# Production Moby GRIB libraries location
use lib "/home/ug/gmaster/projects/moby/prod/lib";

# BioMoby core libraries
use lib "/home/ug/gmaster/projects/moby/biomoby.0.8.1/Perl";

# Also need SOAP and Bioperl libraries that are already installed in the default perl library path
# /usr/local/Install/perl-5.8.5/lib/site_perl/5.8.5/Bio/

# GOstat setup

# GOstat libraries
use lib "/home/ug/gmaster/projects/gostat/lib";
$ENV{GOSTAT} = "/home/ug/gmaster/projects/gostat";

# SGP2 setup
$ENV{SGP2} = "/home/ug/gmaster/sgp2/sgp2_2003/";

1;
