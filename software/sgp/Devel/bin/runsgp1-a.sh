#
# $Id: runsgp1-a.sh,v 1.1 2000-10-17 18:44:42 jabril Exp $
#

SGP1="$SGP_DIR/sgp-1";

ISEQ=$1;
ODIR=$2;
IFILE=$3;

while read gene locus1 locus2;
do
    echo $gene $locus1 $locus2;
    $SGP1 --sgp -d -- \
          "$IFILE/${gene}.${locus1}_${locus2}.tbx" \
          "$ISEQ/$locus1" "$ISEQ/$locus2" > "$ODIR/$gene.sgp1"; # | \
          # perl -ne 'BEGIN{$p=1} /^Predicted/ && ($p=0); print $_ if $p' - > "$ODIR/$gene.sgp1";
done;

exit 0;
