##################################################################################
#
# Configuration file for promoter sequence extraction from Ensembl Database
#
##################################################################################

# bioperl libraries - at least bioperl v1.2 is required
use lib "/usr/local/Install/perl-5.8.6/lib/site_perl/5.8.6/";

# Data::UUID
# required for Ensembl compara
use lib "/home/ug/gmaster/lib/site_perl/5.8.5/i686-linux-thread-multi";

#########################################
#
# Libraries Path - Modify to fit your system !
#
#########################################

$_libraries_path = "/home/ug/gmaster/projects/promoter_extraction/lib";

# Ensembl perl API is release 32 by default
# Initialised here because otherwise it breaks !!
# but actually the script will be using the perl API related to the database release given by the user !

# Libraries Path - Modify to fit your system here as well !!

use lib "/home/ug/gmaster/projects/promoter_extraction/lib/ensembl-32/modules";
use lib "/home/ug/gmaster/projects/promoter_extraction/lib/ensembl-compara-32/modules";

# Ensembl database parameters

$dbhost = "ensembldb.ensembl.org";
$dbuser = "anonymous";
$dbpassword = "";

1;
