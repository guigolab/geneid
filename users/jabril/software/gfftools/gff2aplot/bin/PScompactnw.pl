#!/usr/local/bin/perl -w
# This is perl, v5.6.1 built for i686-linux
# /usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
# 
# USAGE:  PScompactnw.pl "regexp" <STDIN >STDOUT
#
# $Id: PScompactnw.pl,v 1.1 2002-12-24 13:07:36 jabril Exp $
#
use strict;
# my $MaxLen = 80; # Not Used now
my ($clflg,$usflg,$ln,$cln,$go) = (0,0,0,0,0);

my $regexp = shift @ARGV;
my $thori = '@defn';
my $colin = '@text';
my $nw    = '@nl';
my $nwlin = "$nw\n";
my $uslin = '@use';
my $idxln = '@index';
my $thend = '@end code';

(defined($regexp) && $regexp ne '') || do {
    while (<STDIN>) {
        print STDOUT $_;
    };
    exit(0);
};

while (<STDIN>) {
    $_ =~ /^$thori\s+($regexp)/o && do {
        $go = 1;
        print STDOUT "$_$nwlin";
        next;
    };
    $_ =~ /^$thend/o && do {
        $go = 0;
        print STDOUT ($clflg ? "\n$nwlin$_" : "$_");
        next;
    };
    $go || do {
        print STDOUT $_;
        next;
    }; # !$go
    $_ =~ /^$nw/o && next;
    chomp;
    $usflg = 0;
    $_ =~ /^$uslin/o && do {
        print STDOUT ($clflg ? "\n$nwlin" : "")."$_\n$nwlin";
        $usflg = 1;
        $clflg = 0;
        next;        
    }; # $_ =~ /^$uslin/
    $_ =~ /^$idxln/o && do {
        print STDOUT ($clflg ? "\n$nwlin" : '')."$_\n";
        $clflg = 0;
        next;        
    }; # $_ =~ /^$idxln/
    $_ =~ /^$colin/o && do {
        $_ =~ s/^$colin //o;
        $usflg && do {
            $usflg = 0;
            next;
        };
    #    $_ =~ /^\s*$/o && next;
        $_ =~ s/^\s*//o;
        $_ =~ /^\%/o && do { # raw comments are thrown to output directly...
            $_ =~ s/\s*\%\:\%.*$//o && do {
                $clflg = 0;
            };
            $_ =~ s/\s*\%(\-(\%|\>)).*$//o;
            $_ eq '' && do {
                next;
            };
            print STDOUT ($clflg ? "\n$nwlin" : "")."$colin $_\n$nwlin" ;
            $clflg = 0;
            next;
        }; # $_ =~ /^\%/o
        $_ =~ s/\s+/ /og;
        $_ eq '' && do {
            next;
        };
        TYPE: {
            $clflg || ($_ = "$colin $_");
            ($_ =~ s/\s*\%\-(\>|\%).*$//og) && do {
                $_ .= " ";
                $clflg = 1;
                last TYPE;
            };
            $_ =~ s/\s*\%\:\%.*$//og && do {
                $_ .= "\n$nwlin";
                $clflg = 0;
                last TYPE;
            };
            $clflg && ($_ = "\n$nwlin$colin $_");
            $_ .= "\n$nwlin";
            $clflg = 0;
        }; # TYPE
        # $cln = length($_);
        # $cln > 0 || next;
        # $_ =~ /^\s*$/o && next;
        # SIZES: {
        #    $ln == 0 && ( $_ = "$colin $_" );
        #    $clflg || ( $ln = 0, last SIZES );
        ##    $clflg && ( $_ = "\n$nwlin$colin $_" );
        #    $ln += $cln + $clflg;
        #}; # SIZES
        # $_ .= $clflg ? " " : "\n$nwlin";
        print STDOUT "$_";
        # print STDERR "$cln ($ln) : \n$_\n";
    }; # $_ =~ /^$colin/
}; # while

exit(0);
