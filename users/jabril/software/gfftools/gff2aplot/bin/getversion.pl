($PROGRAM,$VERSION,$REVISION,$REVISION_DATE,$LAST_UPDATE) = 
   ( 'gff2aplot','v2.0',
     '$Revision: 1.1 $', #'
     '$Date: 2003-03-04 18:05:37 $', #'
     '$Id: getversion.pl,v 1.1 2003-03-04 18:05:37 jabril Exp $', #'
    );
$REVISION =~ s/\$//og;
$REVISION_DATE =~ s/\$//og;
$tex = (scalar(@ARGV) > 0) ? 1 : 0; # pretty simple...take care
if ($tex) {
    print '\def\PROGname{'.$PROGRAM.'}'.
          '\def\PROGversion{'.$VERSION.'}%'."\n";
} else {
    print "$PROGRAM $VERSION\n";
};
exit(0);
