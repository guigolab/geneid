#!/bin/csh

echo ""
echo "ENSEMBL CVS PASSWORD IS 'CVSUSER'"
echo ""

set release = $1

if ("$release" == "") then
	echo "set the Ensembl release number you wish to retrieve"
	echo "e.g. ensembl_cvs_update.csh 43"
	exit 1;
endif

cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl login
# password is CVSUSER
cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl checkout -r branch-ensembl-${release} ensembl/modules

cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl checkout -r branch-ensembl-${release} ensembl-compara/modules

mv ensembl ensembl-${release}
mv ensembl-compara ensembl-compara-${release}
