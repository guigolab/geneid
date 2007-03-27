# Production Moby GRIB libraries location
use lib "$HOME/projects/moby/prod/lib";

# print STDERR "HOME: $HOME\n";

# BioMoby core libraries
use lib "$HOME/projects/moby/biomoby.0.8.2a/Perl";
# use lib "$HOME/projects/moby/biomoby.perl.22.02.2006/Perl";
# use lib "$HOME/projects/moby/Moby_Perl.07.02.2007/Perl";
# use lib "$HOME/projects/moby/Moby_Perl.27.03.2007/Perl";

# Also need SOAP and Bioperl libraries that are already installed in the default perl library path
# /usr/local/Install/perl-5.8.5/lib/site_perl/5.8.5/Bio

# SOAP v0.60
use lib "$HOME/projects/lib/5.8.8/x86_64-linux-thread-multi";
use lib "$HOME/projects/lib/site_perl/5.8.8";

# Log::Log4perl and depedencies - temporary place !
# Also CGI::Ajax HTML::Template Class::Accessor
# Installed Globally now
# use lib "$HOME/lib/site_perl/5.8.5";
# use lib "$HOME/lib/site_perl/5.8.5/i686-linux-thread-multi";

# GOstat setup

# GOstat libraries
use lib "$HOME/projects/gostat/lib";
$ENV{GOSTAT} = "$HOME/projects/gostat";

# SGP2 setup
$ENV{SGP2} = "$HOME/projects/sgp2";

# Statistic reporting file name
$ENV{STATS_FILE}="$HOME/projects/moby_logs/moby_services_statistics.log";

1;
