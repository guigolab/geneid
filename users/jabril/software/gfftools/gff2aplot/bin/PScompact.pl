#!/usr/local/bin/perl -w
# This is perl, v5.6.1 built for i686-linux
# /usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
# 
# PScompact.pl
#
# $Id: PScompact.pl,v 1.1 2002-12-24 13:07:36 jabril Exp $
#
use strict;
my ($MaxLen,$clflg,$fnflg,$ln,$cln) = (25,0,0,0,0);

while (<STDIN>) {
    my $tt;
    chomp;
    $_ =~ /^\%/o && do { # raw commentss are thrown to output directly...
        $tt = $clflg ? "\n" : '';
        $fnflg = $clflg = $ln = 0;
        print STDOUT "$tt$_\n";
        next;
    }; 
    $_ =~ s/\s+/ /og;
    $clflg && ($_ =~ s/^\s//o);
  TYPE: {
      ($_ =~ s/\s\%\-\>.*$//og) && do {
          $fnflg = $clflg = 1;
          last TYPE;
      };
      $fnflg = 0;
      $clflg = ($_ =~ s/\s\%\-\%.*$//og) ? 1 : 0;
  }; # TYPE
    $cln = length($_);
    $ln > $MaxLen && do {
        $fnflg && ($_ = "  $_"); # indenting opened functions
        $_ = "\n$_";
        $ln = 0;
    };
    $ln += $cln;
    print STDERR "$cln ($ln) : $_\n";
    $tt = $clflg ? " " : "\n";
    print STDOUT "$_$tt";
}; # while

exit(0);
