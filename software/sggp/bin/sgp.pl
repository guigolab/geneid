#!/usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
#
# $Id: sgp.pl,v 1.2 2001-06-06 19:00:58 jabril Exp $
#
use strict;
#
######################################################################
#                               SGP                                  #
######################################################################
#
#     Sinteny based Gene Prediction tool.
#
#     Copyright (C) 2001 - Josep Francesc ABRIL FERRANDO
#                                   Genis PARRA FARRE
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

#### MODULES ####

use Benchmark;
  my @Timer = (new Benchmark);
use Data::Dumper;
local $Data::Dumper::Purity = 0;
local $Data::Dumper::Deepcopy = 1; 

#### CONSTANTS DEFINITION ####

my $PROGRAM = "SGP";
my @tmp_ver = split / +/, 
    ': $Id: sgp.pl,v 1.2 2001-06-06 19:00:58 jabril Exp $ :';
my $VERSION = "v$tmp_ver[3] [@tmp_ver[4,5,7]]";
my ($T,$F) = (1,0); # for 'T'rue and 'F'alse
my $line = ("#" x 60)."\n";
my ($s,$r) = ("###","###\n");
my $Error = "\<\<\<  ERROR  \>\>\> ";
my $Warn  = "\<\<\< WARNING \>\>\> ";
my $spw   = "\<\<\<         \>\>\> ";
my $spl   = "\<\<\<\-\-\-\-\-\-\-\-\-\>\>\>\n";

#### GLOBAL VARIABLES ####

my $DATE = localtime;
my $USER = $ENV{USER};
#
# Error/Warning Message List
my %ErrorList = (
                 
                 UNKNOWN_CL_OPTION => "!!! Error trapped while processing command-line:\n",
                 CMD_LINE_ERROR => "!!! Please, check your command-line options!!!\n",
                 );
#
# Verbose Message List
my %MessageList = (
                   
                   HEADER => $line."$s    RUNNING: $PROGRAM\n".
                                   "$s       USER: $USER\n".
                                   "$s       DATE: $DATE\n".$line,
                   JOB_DONE => $line."$s  TOTAL TIME SPENT --> \%s\n".$line,
                   );

#### MAIN LOOP ####

    &report('HEADER');
    # get options
    &Which_Options();
    #
    &report('JOB_DONE',&timing($T));
    exit(0);

##### SUBS #####

# Parsing command-line options and processing its parameters.
# Checking input sequence files
sub exists_file() {
    my @files = @_;
    my ($n, $r) = (' ', 0);
    foreach $n (@files) {
        $r++ if (-e "$n");
    };
    return $r;
}
sub check_fasta_format() {
    my $file = $_[0];
    my ($n, $c) = (undef, 0);
    open(TMP,"< $file");
    while (<TMP>) {
        next unless /^>/;
        />(\S+)\b/ && do {
            $n = $1;
            $c++;
            next unless $c>1;
        };
        &go_to_die("FATAL ERROR !!! Multiple locus names found.\n  File \'$file\' must contain only one sequence definition.\n");
    }
    &go_to_die("FATAL ERROR !!! There is no '>' line, locus name not found.\n  Please, verify your fasta file \'$file\'\n") unless defined($n);
    return $n;
}
sub Which_Options() {
    $SIG{__WARN__} = sub { &warn('UNKNOWN_CL_OPTION',$_[0]) };
    GetOptions( 
                "1=s"      => \$Seq1         , # seqfile_1
                "2=s"      => \$Seq2         , # seqfile_2
                "g=s"      => \$geneid_opt   , # geneid options      
                "P=s"      => \$geneid_param , # geneid parameter file 
                "o=s"      => \$blast_opt    , # tblastx options 
                "c=f"      => \$score_cutoff , # tblastx score cutoff
                "s=f"      => \$shrink       , # shrink hsp's by
                "t=s"      => \$tbx          , # read tblastx from file
                "f=s"      => \$hsp          , # read HSP files in directory
                "k=s"      => \$ofn          , # intermediate filename
                "p=s"      => \$ps_output    , # postscript output 
                "v"        => \$verbose_flg  , # verbose    
                "h|help|?" => \$help_flg     , # print help
                ) || (&warn('CMD_LINE_ERROR'), exit(1));
    $SIG{__WARN__} = 'DEFAULT';
    &prt_Help if $help_flg;
    &prt_Header("Processsing Command-Line Options") if $verbose_flg;
    # do seq files exists
    my $file_number = &exists_file($Seq1, $Seq2);
    &go_to_die("FATAL ERROR!!! Two sequences are needed (-1 and -2 options are mandatory).\n")
        unless $file_number == 2; 
    # extract basenames
    $Seq1_Name = basename($Seq1);
    $Seq2_Name = basename($Seq2);
    # check if files are provided in fasta format and get locus names
    $Loc1 = &check_fasta_format($Seq1);
    $Loc2 = &check_fasta_format($Seq2);
    &go_to_die("FATAL ERROR!!! Locus1($Loc1) have the same name as Locus2($Loc2).\n  Sequences \'$Seq1_Name\' and \'$Seq2_Name\' must have different locus names.\n")
        if $Loc1 eq $Loc2;
    # tmpfiles prefix
    $SGPtmp1 = "$SGPTMP.$Seq1_Name";
    $SGPtmp2 = "$SGPTMP.$Seq2_Name";
    $SGPtmpG = "$SGPTMP.${Loc1}_${Loc2}"; 
    do {
        &copy_files("$hsp${Loc1}_${Loc2}.srQ", "$SGPtmpG.srQ");
        &copy_files("$hsp${Loc1}_${Loc2}.srS", "$SGPtmpG.srS");
        &copy_files("$hsp${Loc1}_${Loc2}.aln", "$SGPtmpG.aln");
        $hsps_flg = 0;
    } if defined($hsp);
    do {
        &copy_files("$tbx", "$SGPTMP.tbx");
        $blast_flg = 0;
    } if defined($tbx);
    # $geneid_flg = 0;
    $savefiles_flg = 1 if defined($ofn);
    $plots_flg = 0 unless defined($ps_output);
    # setting other variables
    $GeneidOptions = "$geneid_opt"   if defined($geneid_opt);
    $GeneidParam   = "$geneid_param" if defined($geneid_param);
    $BlastOptions  = "$blast_opt"    if defined($blast_opt);
    $verbose_str = "" if $verbose_flg;
} # sub Which_Options

$start = new Benchmark;
$sys_rpt = system("$cmdline");
$end = new Benchmark;
&check_sys_result($sys_rpt);
$ttime = timestr(timediff($end,$start));
print STDERR "         Time Spent: $ttime \n";



#
sub check_sys_result() {
    my $runOK = $F;
    my $prog_exit = 0xffff & $_[0];
    printf STDERR "### Command returned %#04x : ", $prog_exit;
    if ($prog_exit == 0) {
        print STDERR "ran with normal exit ...\n";
        $runOK = $T;
    } elsif ($prog_exit == 0xff00) {
        print STDERR "command failed: $! ...\n";
    } elsif (($prog_exit & 0xff) == 00) {
        $prog_exit >>= 8;
        print STDERR "ran with non-zero exit status $prog_exit ...\n";
    } else {
        print STDERR "ran with ";
        if ($prog_exit &   0x80) {
            $prog_exit &= ~0x80;
            print STDERR "coredump from ";
        };
        print STDERR "signal $prog_exit ...\n";
    };
    return $runOK;
} # check_sys_result()
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
#
sub fill_right() { $_[0].($_[2] x ($_[1] - length($_[0]))) }
sub fill_left()  { ($_[2] x ($_[1] - length($_[0]))).$_[0] }
sub fill_mid()   { 
    my $l = length($_[0]);
    my $k = int(($_[1] - $l)/2);
    ($_[2] x $k).$_[0].($_[2] x ($_[1] - ($l+$k)));
} # fill_mid
sub report() { print STDERR sprintf($MessageList{ shift @_ },@_) }
sub warn() { print STDERR sprintf($ErrorList{ shift @_ }, @_) }

#### EOF ####

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

Genis Parra     : B<gparra@imim.es>

B<$PROGRAM> is under GNU-GPL (C) 2000
