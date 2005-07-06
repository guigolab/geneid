# MOBY GRIB libraries
# They are in GRIB CVS, under projects
use lib "/home/ug/arnau/cvs/GRIB/projects/moby/lib";

# BioMOBY libraries
use lib "/home/ug/arnau/cvs/moby-live/Perl";

# Also need SOAP and Bioperl libraries that are already installed in the default perl library path
# /usr/local/Install/perl-5.8.5/lib/site_perl/5.8.5/Bio/

# GOstat setup

# GOstat libraries
use lib "/home/ug/gmaster/projects/gostat/lib";
$ENV{GOSTAT} = "/home/ug/gmaster/projects/gostat";

# SGP2 setup
$ENV{SGP2} = "/home/ug/gmaster/sgp2/sgp2_2003/";

1;