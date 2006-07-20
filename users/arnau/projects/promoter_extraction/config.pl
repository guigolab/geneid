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
# Libraries Path - Modify to fit your system and the Ensembl release you want to use!
#
#########################################

$_libraries_path = "/home/ug/gmaster/projects/promoter_extraction/lib";
# $_libraries_path = "/home/ug/arnau/no_backup/cvs";

use lib "/home/ug/gmaster/projects/promoter_extraction/lib/ensembl-39/modules";
use lib "/home/ug/gmaster/projects/promoter_extraction/lib/ensembl-compara-39/modules";

# Latest Ensembl release

#########################################################
#
# Modify every time there is a new Ensembl Release !!!!
#
#########################################################

$latest_release = "39";

# Ensembl database connection parameters
# Two sets of parameters to allow connecting to a different server depending on which release the user is asking !!!
# The default is ensembldb server but you can specify another server below

$dbhost_default = "ensembldb.ensembl.org";
$dbuser_default = "anonymous";
$dbpassword_default = "";

$dbhost_latest = "ensembldb.ensembl.org";
$dbuser_latest = "anonymous";
$dbpassword_latest = "";

1;
