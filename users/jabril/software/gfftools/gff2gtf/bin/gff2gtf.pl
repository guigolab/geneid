#!/usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
#line 658 "/home/ug/jabril/development/softjabril/gfftools/gff2gtf/gff2gtf.nw"
# #------------------------------------------------------------#
# #                        gff2gtf                             #
# #------------------------------------------------------------#
# 
#         Converting GFF formated records into GTF2.
# 
#     Copyright (C) 2001 - Josep Francesc ABRIL FERRANDO  
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
# 
# #----------------------------------------------------------------#
#line 652 "/home/ug/jabril/development/softjabril/gfftools/gff2gtf/gff2gtf.nw"
# $Id: gff2gtf.pl,v 1.1 2001-09-06 18:27:41 jabril Exp $
#line 500 "/home/ug/jabril/development/softjabril/gfftools/gff2gtf/gff2gtf.nw"
#
use strict;
#line 250 "/home/ug/jabril/development/softjabril/gfftools/gff2gtf/gff2gtf.nw"
#
# MODULES
#

#
# VARIABLES
#

#
# MAIN LOOP
#
#line 291 "/home/ug/jabril/development/softjabril/gfftools/gff2gtf/gff2gtf.nw"
#
# gff2gtf.pl < infile > outfile
#
#     Converting GFF formated records into GTF
#
# use lib qw( /usr/lib/perl5/site_perl/5.005/ /usr/lib/perl5/5.00503/ ) ;
use Benchmark;
# use Data::Dumper;
# local $Data::Dumper::Purity = 0;
# local $Data::Dumper::Deepcopy = 1;
my ($T,$F) = (1,0);
my @Timer = (new Benchmark);
my $PROGRAM = 'gff2gtf.pl';
(my $VERSION = ' $Name: not supported by cvs2svn $ ') =~ s/^.*?:\s*|(\s*\$?\s*)?$//og;
my $DATE = localtime;
my $USER = defined($ENV{USER}) ? $ENV{USER} : '??????';
my $host = `hostname`;
chomp($host);
my $line = ('#' x 80)."\n";
my $s = '### ';
#
print STDERR << "+++EOR+++";
$line$s\n$s Running $PROGRAM\n$s
$s HOST: $host
$s USER: $USER
$s DATE: $DATE\n$s\n$line$s
+++EOR+++

my %genes = ();
my %gene_list = ();
my @f = ();
my $gcnt = 0;
while (<>) {
#     /^\# Gene/o && do {
#         $gcnt++;
#         $genes{$gcnt}{COMMENT} = $_;
#         next;
#     };
    my ($_cnt,$key);
    next if /^\#/o;
    next if /^\s*$/o;
    chomp;
    @f = split /\s+/og, $_;
    $key = $f[8].$f[6];
    print STDERR "### $key\n";
    defined($gene_list{$key}) || do { # group + strand
        $gcnt++;
        $gene_list{$key} = $gcnt;
        $genes{$gcnt}{NAME} = $f[8];
        $genes{$gcnt}{START} = $f[3];
        $genes{$gcnt}{STOP} = $f[4];
        $genes{$gcnt}{STRAND} = $f[6];
        $genes{$gcnt}{LEN} = ($f[4] - $f[3] + 1);
        push @{ $genes{$gcnt}{FEATURES} }, [ @f[0..8] ];
        next;
    };
    $_cnt = $gene_list{$key};
    $genes{$_cnt}{LEN} += ($f[4] - $f[3] + 1);
    $genes{$_cnt}{START} = &min($f[3],$genes{$_cnt}{START});
    $genes{$_cnt}{STOP}  = &max($f[4],$genes{$_cnt}{STOP});
    push @{ $genes{$_cnt}{FEATURES} }, [ @f[0..8] ];
}; # while

# print STDERR Dumper(\%genes);

my $GTFrec = ("\%s\t" x 8)."\%s\n";
foreach my $i (1..$gcnt) {
    my ($_gnam,$_gnum,$_glen,$_gstr,$_gori,$_gend);
    $_gnam = $genes{$i}{NAME};
    $_gnum = scalar(@{ $genes{$i}{FEATURES} });
    $_gstr = $genes{$i}{STRAND};
    $_glen = $genes{$i}{LEN};
    $_gori = ($_gstr ne "-") ? $genes{$i}{START} : $genes{$i}{STOP};
    $_gend = ($_gstr ne "-") ? $genes{$i}{STOP} : $genes{$i}{START};
    print STDOUT "\# Gene:  $_gnam  $_gstr  $_gori  $_gend  (${_glen}bp)".
                 ",  $_gnum element".(($_gnum > 1) ? "s" : "")."\n";
    my ($main,$tail_O,$tail_E,$groupflg,$gtfgroup) = ('','','',$T,'');
    foreach my $j (0..$#{ $genes{$i}{FEATURES} }) {
        my ($seq,$source,$feat,$start,$end,
            $score,$strand,$frame,$group,
            $cdsfeat,$cdsori,$cdsend,$t);
        ($seq,$source,$feat,$start,$end,
            $score,$strand,$frame,$group) =
                @{ $genes{$i}{FEATURES}[$j] };
        $score = '.';
        $groupflg && do {
            $group =~ s/gene_//o;
            $t = &fill_left($group,3,"0");
            $gtfgroup = "gene_id $t; transcript_id $t.1";
            $groupflg = $F;
        }; # $groupflg
      FEATS: {
          $feat eq 'Single' && do {
              ($cdsfeat,$cdsori,$cdsend) = &get_start($start,$end,$strand);
              $tail_O = sprintf($GTFrec,$seq,$source,$cdsfeat,$cdsori,$cdsend,
                                $score,$strand,'0',$gtfgroup);
              ($cdsfeat,$cdsori,$cdsend) = &get_final($start,$end,$strand);
              $tail_E = sprintf($GTFrec,$seq,$source,$cdsfeat,$cdsori,$cdsend,
                                $score,$strand,'0',$gtfgroup);
              last FEATS;
          }; # Single
          $feat eq 'First' && do {
              ($cdsfeat,$cdsori,$cdsend) = &get_start($start,$end,$strand);
              $tail_O = sprintf($GTFrec,$seq,$source,$cdsfeat,$cdsori,$cdsend,
                                $score,$strand,'0',$gtfgroup);              
              last FEATS;
          }; # First
          $feat eq 'Terminal' && do {
              ($cdsfeat,$cdsori,$cdsend) = &get_final($start,$end,$strand);
              $tail_E = sprintf($GTFrec,$seq,$source,$cdsfeat,$cdsori,$cdsend,
                                $score,$strand,'0',$gtfgroup);
                  # only if nucleotides of stop codon not included
                  #     $strand eq '+' && ($end = $end + 3);
                  #     $strand eq '-' && ($start = $start - 3);
              # last FEATS;
          }; # Terminal
        }; # FEATS
        $feat = 'CDS';
        $main .= sprintf($GTFrec,$seq,$source,$feat,$start,$end,
                         $score,$strand,$frame,$gtfgroup); 
    }; # foreach $j
    print STDOUT "$main$tail_O$tail_E";
}; # foreach $i

my $total_time = &timing($T);
print STDERR << "+++EOR+++";
$s\n$line$s\n$s $PROGRAM FINISHED\n$s
$s TOTAL TIME: $total_time\n$line
+++EOR+++

#line 278 "/home/ug/jabril/development/softjabril/gfftools/gff2gtf/gff2gtf.nw"
exit(0);
#line 262 "/home/ug/jabril/development/softjabril/gfftools/gff2gtf/gff2gtf.nw"
#
# FUNCTIONS
#
#line 423 "/home/ug/jabril/development/softjabril/gfftools/gff2gtf/gff2gtf.nw"
#
sub get_start() {
    my ($o,$e,$s) = @_;
    my $str = "start_codon";
    $s eq '+' && do {
       return $str, $o, ($o + 2);
    }; # forward
    $s eq '-' && do {
       return $str, ($e - 2), $e;
    }; # reverse
    die("### ERROR: Strand not defined... ($o $e : $s)... $!");
} # get_start
#
sub get_final() {
    my ($o,$e,$s) = @_;
    my $str = "stop_codon";
    $s eq '+' && do {
        return $str, ($e - 2), $e;
           # only if nucleotides of stop codon not included
           #     return $str, ($e + 1), ($e + 3);
    }; # forward
    $s eq '-' && do {
        return $str, $o, ($o + 2);
           # only if nucleotides of stop codon not included
           #     return $str, ($o - 3), ($o -1);
    }; # reverse
    die("### ERROR: Strand not defined... ($o $e : $s)... $!");
} # get_final
#
sub fill_left() { ($_[2] x ($_[1] - length($_[0]))).$_[0] } 
#
sub timing() {
    push @Timer, (new Benchmark);
    # partial time
    $_[0] ||
        (return timestr(timediff($Timer[$#Timer],$Timer[($#Timer - 1)])));
    # total time
    return timestr(timediff($Timer[$#Timer],$Timer[0]));
} # timing
#line 565 "/home/ug/jabril/development/softjabril/gfftools/gff2gtf/gff2gtf.nw"
#
sub max() {
    my $z = shift @_;
    foreach my $l (@_) { $z = $l if $l > $z };
    return $z;
} # max
sub min() {
    my $z = shift @_;
    foreach my $l (@_) { $z = $l if $l < $z };
    return $z;
} # min
