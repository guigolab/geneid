#
# $Id: runsgp2-a.sh,v 1.4 2000-10-19 11:23:25 jabril Exp $
#

SGP2="$SGP/bin/sgp2.pl";

ISEQ=$1;
ODIR=$2;
IFILE=$3;

while read gene locus1 locus2;
do
    echo $gene $locus1 $locus2;
    $SGP2 -v -1 "$ISEQ/$locus1" -2 "$ISEQ/$locus2"\
          -k "$IFILE/$gene." -p "$IFILE/$gene." > "$ODIR/$gene.sgp";
done;

exit 0;
