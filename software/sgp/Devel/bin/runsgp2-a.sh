ISEQ=$1
ODIR=$2

IFILE=$3

SGP2=sgp2
while read gene locus1 locus2
do
      echo $gene $locus1 $locus2
      $SGP2 -v -1 $ISEQ/$locus1 -2 $ISEQ/$locus2 -k $IFILE/ > $ODIR/$gene
done

