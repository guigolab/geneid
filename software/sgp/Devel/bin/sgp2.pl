#!/usr/bin/perl -w
#line 93 "./sgp2.nw"
#line 542 "./sgp2.nw"
######################################################################
#                               sgp2                                 #
######################################################################
#
#     Sinteny based Gene Prediction tool.
#
#     Copyright (C) 2000 - Josep Francesc ABRIL FERRANDO
#                                 Roderic GUIGO SERRA
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
######################################################################
#
#line 94 "./sgp2.nw"
my $Start = time;
#line 103 "./sgp2.nw"
my $PROGRAM = "sgp2";
my @tmp_ver = split / +/, ' $Id: sgp2.pl,v 1.5 2000-10-11 19:10:27 jabril Exp $ ';
my $VERSION = "v$tmp_ver[3] [$tmp_ver[4] $tmp_ver[5] $tmp_ver[7]]";
#line 111 "./sgp2.nw"
use strict;
#line 186 "./sgp2.nw"
use Getopt::Long;
Getopt::Long::Configure qw/ bundling pass_through /;
#line 267 "./sgp2.nw"
use Pod::Text;

#line 65 "./sgp2.nw"
### VARS ###

#line 203 "./sgp2.nw"
# GetOptions Variables
my ( $seq1, $seq2, $geneid_opt, $geneid_param, 
     $blast_opt, $score_cutoff, $shrink, $tbx,
     $hsp, $ofn, $ps_output, $verbose_flg, $help_flg );
#line 212 "./sgp2.nw"
my ($blast_flg, $hsps_flg, $geneid_flg, $plots_flg) = (1, 1, 1, 1);
#line 362 "./sgp2.nw"
# Processes
my $status
# Blast 
my $BlastProgram = "tblastx";
my $PressdbProgram = "pressdb";
my $BlosumMatrix = "$SGP2param/blosum62mod"; 
my $BlastOptions = "-matrix $BlosumMatrix -hspmax=10000 -nogap";
my @pressdb = ( "$PressdbProgram" );
my @blast = ( "$BlastProgram" "$BlastOptions" );
#line 411 "./sgp2.nw"
# geneid
my $GeneidProgram = "$SGP2bin/geneid.v1.0-sgp/bin/geneid";
my $GeneidParam   = "$SGP2bin/geneid.v1.0-sgp/param/human3iso.param";
my $GeneidOptions = "-GP $GeneidParam";
my @geneid = ( "$GeneidProgram" "$GeneidOptions" );
#line 439 "./sgp2.nw"
# gff2ps and aplot
#line 136 "./sgp2.nw"
# some defaults
my $file_number = 0;
my $S_CUTOFF = 50;
my $SCF = 12; # substract to tblastx scores S_CUTOFF - SCF;
#line 145 "./sgp2.nw"
# setting paths
my $SGP2   = '';
my $SGP2bin = "$SGP2/bin";
my $SGP2param = "$SGP2bin/param";
my $TMP    = '/tmp';
my $SGPTMP = "$TMP/sgp2_$$";
#line 164 "./sgp2.nw"
$SIG{INT}  = \&trap_signals;
$SIG{QUIT} = \&trap_signals;
$SIG{TERM} = \&trap_signals;
$SIG{CHLD} = 'IGNORE';
sub trap_signals() {
    # opendir(TDIR, $TMP) or die "Can't open directory $TMP: $!";
    # while (defined(my $file = readdir(TDIR)))
    &clean_tmp();
    die "WARNING !!! $PROGRAM has been stopped by user.";
}

#line 70 "./sgp2.nw"
### SUBS ###

#line 273 "./sgp2.nw"
# Prints help 
sub prt_Help() {
    my $tmp_pod_file = "$TMP/sgp.pod";
    open(KI, "> $tmp_pod_file");
    while (<DATA>) {
        s/\$PROGRAM/$PROGRAM/g ;
        s/\$VERSION/$tmp_ver[3]/g ;
        print KI $_ ;
    };
    close(KI);
    pod2text($tmp_pod_file);
    unlink($tmp_pod_file) or die "Can't delete $tmp_pod_file: $!\n";
    exit(1);
}
#line 297 "./sgp2.nw"
# Deleting temporary files on TMP
sub clean_tmp() {
    my @files = glob($SGPTMP."*");
    unlink(@files) or warn "Can't unlink @files : $? $!";
}
#line 307 "./sgp2.nw"
# Reporting IN/OUT progress.
sub prt_progress() {
    $verbose_flg && do {
        print STDERR ".";
        (($_[0] % 50) == 0) && print STDERR "[".&fill_left($_[0],6,"0")."]\n";
    };
}
#
sub prt_foeprg() {
    $verbose_flg && ((($_[0] % 50) != 0) && print STDERR "[".&fill_left($_[0],6,"0")."]\n" );
}
# Get a fixed length string from a given string and filling char/s.
sub fill_right() { $_[0].($_[2] x ($_[1] - length($_[0]))) }
sub fill_left()  { ($_[2] x ($_[1] - length($_[0]))).$_[0] }
# returns the max value from input array
sub max() { my ($z) = shift @_; my $l; foreach $l (@_) { $z = $l if $l > $z ; }; $z; } 
#line 334 "./sgp2.nw"
# Timing.
sub get_exec_time() {
    $verbose_flg && do {
        my $End = $_[0];
        my ($c,$s,$m,$h,$r);
        $r = $End - $Start;
        $s = $r % 60;
        $r = ($r - $s) / 60;
        $m = $r % 60;
        $r = ($r - $m) / 60;
        $h = $r % 24;
        ($s,$m,$h) = (&fill_left($s,2,"0"),&fill_left($m,2,"0"),&fill_left($h,2,"0"));
print STDERR <<EOF;
##
##########################################################
## \"$PROGRAM\"  Execution Time:  $h:$m:$s
##########################################################
EOF
    };
}
#line 218 "./sgp2.nw"
# Parsing command-line options and processing its parameters.
#line 242 "./sgp2.nw"
# Checking input sequence files
sub exits_file() {
    my @files = @_;
    my ($n, $r) = ('', 0);
    foreach $n @files {
        $r++ if (-f $n && -T _);
    };
    return $r;
}
#line 220 "./sgp2.nw"
sub Which_Options() {
    GetOptions( 
                "1=s"      => \$seq1         , # seqfile_1
                "2=s"      => \$seq2         , # seqfile_2
                "g"        => \$geneid_opt   , # geneid options      
                "P"        => \$geneid_param , # geneid parameter file 
                "o"        => \$blast_opt    , # tblastx options 
                "c"        => \$score_cutoff , # tblastx score cutoff 
                "s"        => \$shrink       , # shrink hsp's by
                "t"        => \$tbx          , # read tblastx from file
                "f"        => \$hsp          , # read HSP files in directory
                "k"        => \$ofn          , # intermediate filename
                "p"        => \$ps_output    , # postscript output 
                "v"        => \$verbose_flg  , # verbose    
                "h|help|?" => \$help_flg     , # print help
                );
    &prt_Help if $help_flg;
    
#line 254 "./sgp2.nw"
$file_number = &exits_file($seq1, $seq2);
die "FATAL ERROR!!! Two sequences are needed (-1 and -2 options are mandatory): $!"
    unless $file_number == 2; 
$blast_flg  = 0;
$hsps_flg   = 0;
$geneid_flg = 0;
$plots_flg  = 0;
#line 238 "./sgp2.nw"
}; # sub Which_Options

#line 75 "./sgp2.nw"
### MAIN SUBS ###

#line 390 "./sgp2.nw"
# Runnig blast
sub Run_Blast() {
    $status = system { $pressdb[0] } @pressdb;
    die "FATAL ERROR !!! $? $!\n  Unsuccessful PRESSDB command : @pressdb \n" unless $status == 0;
    $status = system { $blast[0] } @blast;
    die "FATAL ERROR !!! $? $!\n  Unsuccessful BLAST command : @blast \n" unless $status == 0;
} # END_SUB: Run_Blast
#line 403 "./sgp2.nw"
# Extracting HSPs
sub Extract_HSP() {
} # END_SUB: Extract_HSP
#line 429 "./sgp2.nw"
# Running geneid
sub Run_geneid() {
    $status = system { $geneid[0] } @geneid;
    die "FATAL ERROR !!! $? $!\nUnsuccessfull GENEID command : @geneid \n" unless $status == 0;
} # END_SUB: Run_geneid
#line 442 "./sgp2.nw"
# Running gff2ps and aplot
sub Make_Plots() {
    $status = system { $gff2ps[0] } @gff2ps;
    die "FATAL ERROR !!! $? $!\nUnsuccessfull GFF2PS command : @gff2ps \n" unless $status == 0;
    $status = system { $aplot[0] } @aplot;
    die "FATAL ERROR !!! $? $!\nUnsuccessfull APLOT command : @aplot \n" unless $status == 0;
} # END_SUB: Make_Plots

#line 82 "./sgp2.nw"
### MAIN LOOP ###

#line 119 "./sgp2.nw"
# get options
&Which_Options();
# main external program calls
&Run_Blast()   if $blast_flg;
&Extract_HSP() if $hsps_flg;
&Run_geneid()  if $geneid_flg;
&Make_Plots()  if $plots_flg; 
# timing
&get_exec_time(time);
# removing temporary files and exit
&clean_tmp();
exit(1);

#line 86 "./sgp2.nw"
### EOF ###

#line 579 "./sgp2.nw"
__DATA__
=head1 NAME

    
$PROGRAM ($VERSION) - Improving Gene Prediction with Sinteny.

=head1 SYNOPSIS

    
    $PROGRAM [-hv] [-o 'options'] [-g 'options'] \
      [-P filename] [-p filename] [-k filename] \
      [-c value] [-s value] -1 seqfile_1 -2 seqfile_2

=head1 DESCRIPTION

=head1 OPTIONS

    

=over 4

=item B<-1> I<seqfile_1>

input file for first species.

=item B<-2> I<seqfile_2>

input file for second species.

=item B<-g>

geneid options

=item B<-o>

tblastx options

=item B<-c> I<value>

tblastx score cuttof

=item B<-s> I<value>

shrink hsp\'s by value

=item B<-t> I<filename>

read tblastx file

=item B<-f> I<prefix>

read hsp gff files with in directory prefix and extension .hsp-rs

=item B<-k> I<prefix>

keep intermediate files with prefix

=item B<-p> I<filename>

ps output in filename file 

=item B<-P> I<filename>

geneid parameter file

=item B<-v>

verbose mode

=item B<-h>

produces this message

=back

=head1 FILES

=head1 DIAGNOSTICS

=head1 REQUIRES

=head1 BUGS

    
Report any problem to: B<jabril@imim.es>

=head1 AUTHOR

    
Roderic Guigo   : B<rguigo@imim.es>

Josep F. Abril  : B<jabril@imim.es>

B<$PROGRAM> is under GNU-GPL (C) 2000
