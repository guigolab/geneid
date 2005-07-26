# Ensembl perl API is release 31 by default 
# but actually the script will be using the perl API related to the database release given by the user !

use lib "/home/ug/gmaster/projects/promoter_extraction/lib/ensembl-31/modules";

# Ensembl database parameters

$dbhost = "ensembldb.ensembl.org";
$dbuser = "anonymous";
