#
# $Id: runsgp2-a.sh,v 1.3 2000-10-16 18:02:03 jabril Exp $
#

ISEQ=$1;
ODIR=$2;

IFILE=$3;

SGP2=$SGP/bin/sgp2.pl;
while read gene locus1 locus2;
do
    echo $gene $locus1 $locus2;
    $SGP2 -v -1 "$ISEQ/$locus1" -2 "$ISEQ/$locus2"\
          -k "$IFILE/$gene." -p "$IFILE/$gene." > "$ODIR/$gene.sgp";
done;

exit 0;
