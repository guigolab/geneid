#!/bin/bash

toerr () {
  perl -e '
    $cdon=0;
    while (<STDIN>) {
      $_ =~ /^@(begin|end) code/o && do {
        $cdon = 1 - $cdon;
        print STDERR $_ unless $cdon;
      };
      $cdon && do { print STDERR $_ };
      print STDOUT $_;
    };
  ' $1 2> $2;
}

toerr - $1 |
  $BIN/PScompactnw.pl "POSTSCRIPT|PSFunction|PSVariables" - | \
  toerr - $2 ; 
#
