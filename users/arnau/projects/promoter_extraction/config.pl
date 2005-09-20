# Ensembl perl API is release 31 by default 
# but actually the script will be using the perl API related to the database release given by the user !

use lib "/home/ug/gmaster/projects/promoter_extraction/lib/ensembl-32/modules";
use lib "/home/ug/gmaster/projects/promoter_extraction/lib/ensembl-compara-32/modules";

# Data::UUID
use lib "/home/ug/gmaster/lib/site_perl/5.8.5/i686-linux-thread-multi";

# Ensembl database parameters

$dbhost = "ensembldb.ensembl.org";
$dbuser = "anonymous";
