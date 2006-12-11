#!/bin/bash

# Moby services testing script

export PERL5LIB=$HOME/lib/biomoby.0.8.2a/Perl

server=$1
errfile='/home/ug/arnau/cron/logs/log.err'
outfile='/home/ug/arnau/cron/logs/log.out'

/home/ug/arnau/cvs/GRIB/projects/moby/scripts/testServices/testServices.pl -x $server 2> $errfile > $outfile

result=`find . -name "log.err" -empty`

if [ "$result" != "" ]; then
    echo "fine" > /dev/null
else
    echo "send an email to report the problems"
    mail -s "Moby services report" akerhornou@imim.es < $errfile
fi
