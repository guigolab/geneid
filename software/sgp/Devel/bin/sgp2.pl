#!/usr/bin/perl -w
#line 96 "./sgp2.nw"
#line 785 "./sgp2.nw"
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
#line 97 "./sgp2.nw"
my $Start = time;
#line 106 "./sgp2.nw"
my $PROGRAM = "sgp2";
my @tmp_ver = split / +/, ' $Id: sgp2.pl,v 1.7 2000-10-13 19:56:06 jabril Exp $ ';
my $VERSION = "v$tmp_ver[3] [$tmp_ver[4] $tmp_ver[5] $tmp_ver[7]]";
#line 114 "./sgp2.nw"
use strict;
#line 191 "./sgp2.nw"
use Getopt::Long;
Getopt::Long::Configure qw/ bundling pass_through /;
#line 279 "./sgp2.nw"
use File::Basename;
#line 347 "./sgp2.nw"
use Pod::Text;
#line 375 "./sgp2.nw"
use File::Copy;

#line 65 "./sgp2.nw"
### VARS ###

#line 141 "./sgp2.nw"
# some defaults
my $S_CUTOFF = 50;
my $SCF = 12; # substract to tblastx scores S_CUTOFF - SCF;
#line 149 "./sgp2.nw"
# setting paths
my $SGP2   = '/home3/jabril/development/sgp/Devel';
my $SGP2bin = "$SGP2/bin";
my $SGP2param = "$SGP2bin/param";
my $TMP    = '/tmp';
my $TMPROOT = "sgp2_$$";
my $SGPTMP = "$TMP/$TMPROOT";
#line 169 "./sgp2.nw"
$SIG{INT}  = \&trap_signals;
$SIG{QUIT} = \&trap_signals;
$SIG{TERM} = \&trap_signals;
$SIG{CHLD} = 'IGNORE';
sub trap_signals() {
    # opendir(TDIR, $TMP) or &go_to_die("Can't open directory $TMP: $!");
    # while (defined(my $file = readdir(TDIR)))
    &clean_tmp();
    &go_to_die("WARNING !!! $PROGRAM has been stopped by user.");
}
#line 208 "./sgp2.nw"
# GetOptions Variables
my ( $Seq1, $Seq2, $geneid_opt, $geneid_param, 
     $blast_opt, $score_cutoff, $shrink, $tbx,
     $hsp, $ofn, $ps_output, $savefiles_flg, $verbose_flg, $help_flg
     ) = ( undef, undef, undef, undef, undef, undef,
           undef, undef, undef, undef, undef, 0, 0, 0 );
#line 285 "./sgp2.nw"
# sequence files
my ($Seq1_Name, $Seq2_Name, $Loc1, $Loc2, $SGPtmp1, $SGPtmp2, $SGPtmpG);
#line 311 "./sgp2.nw"
# file flags
my ($blast_flg, $hsps_flg, $geneid_flg, $plots_flg) = (1, 1, 1, 1);
#line 460 "./sgp2.nw"
# Processes
my $status;
# Blast 
my $BlastProgram = "tblastx";
my $PressdbProgram = "pressdb";
my $BlastMatrix = "$SGP2param/blosum62mod"; 
my $BlastOptions = "-warnings -hspmax 10000 -nogap";
#line 524 "./sgp2.nw"
# Extract HSPs
my $GetSRsProgram = "$SGP2bin/GetSRsAln.pl";
#line 619 "./sgp2.nw"
# geneid
my $GeneidProgram = "$SGP2bin/geneid.v1.0-sgp/bin/geneid";
my $GeneidParam   = "$SGP2bin/geneid.v1.0-sgp/param/human3iso.param";
my $GeneidOptions = "-G";
#line 662 "./sgp2.nw"
# gff2ps and aplot
my $Gff2psProgram = "$SGP2bin/gff2ps";
my $AplotProgram  = "$SGP2bin/gff2aplot";

#line 70 "./sgp2.nw"
### SUBS ###

#line 318 "./sgp2.nw"
# copying files
sub copy_files() {
  	&go_to_die("FATAL ERROR!!! Couln't find '$_[0]' file. $!\n")
         unless &exists_file($_[0]);
    copy(@_) or &copy_error(@_);
}
sub copy_error() { &go_to_die("FATAL ERROR !!! Could not copy file '$_[0]' to '$_[1]': $!\n") }
#line 353 "./sgp2.nw"
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
    unlink($tmp_pod_file) or &go_to_die("Can't delete $tmp_pod_file: $!\n");
    exit(1);
}
#line 382 "./sgp2.nw"
# Deleting temporary files on TMP
sub clean_tmp() {
    opendir(DIR,"$TMP");
    my @files = map { "$TMP/$_" } grep { /^$TMPROOT/ } readdir(DIR);
    closedir(DIR);
    $#files and 
        ( unlink(@files)
          and die "********* Temporary files were deleted *********\n"
          or die "********* Can't unlink @files : $! *********"
        );
    die "********* There are no temporary files in $TMP *********";
}
# writing die messages to STDERR and clean_tmp before exit.
sub go_to_die() { (print STDERR "@_") && &clean_tmp() }
#line 402 "./sgp2.nw"
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
sub fill_mid()   { my $l = length($_[0]); my $k = int(($_[1] - $l)/2); ($_[2] x $k).$_[0].($_[2] x ($_[1] - ($l+$k))) }
# returns the max value from input array
sub max() { my ($z) = shift @_; my $l; foreach $l (@_) { $z = $l if $l > $z ; }; $z } 
# section headers to STDERR
sub prt_Header() { print STDERR ("*" x 80)."\n** ".&fill_mid("@_",74," ")." **\n".("*" x 80)."\n" }
#line 432 "./sgp2.nw"
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
#line 219 "./sgp2.nw"
# Parsing command-line options and processing its parameters.
#line 247 "./sgp2.nw"
# Checking input sequence files
sub exists_file() {
    my @files = @_;
    my ($n, $r) = ('', 0);
    foreach $n (@files) {
        $r++ if (-e $n);
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
		&go_to_die("FATAL ERROR !!! Multiple locus names found.\n  File '$file' must contain only one sequence definition.\n");
    }
    &go_to_die("FATAL ERROR !!! There is no '>' line, locus name not found.\n  Please, verify your fasta file '$file'\n") unless defined($n);
    return $n;
}
#line 221 "./sgp2.nw"
sub Which_Options() {
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
                );
    &prt_Help if $help_flg;
    &prt_Header("Processsing Command-Line Options") if $verbose_flg;
    
#line 290 "./sgp2.nw"
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
&go_to_die("FATAL ERROR!!! Locus1($Loc1) have the same name as Locus2($Loc2).\n  Sequences '$Seq1_Name' and '$Seq2_Name' must have different locus names.\n")
    if $Loc1 eq $Loc2;
# tmpfiles prefix
$SGPtmp1 = $SGPTMP.$Seq1_Name;
$SGPtmp2 = $SGPTMP.$Seq2_Name;
$SGPtmpG = "${SGPTMP}${Loc1}_${Loc2}"; 
#line 328 "./sgp2.nw"
defined($hsp) and do {
    &copy_files("$hsp${Loc1}_${Loc2}.srQ", "$SGPtmpG.srQ");
    &copy_files("$hsp${Loc1}_${Loc2}.srS", "$SGPtmpG.srS");
    &copy_files("$hsp${Loc1}_${Loc2}.aln", "$SGPtmpG.aln");
    $hsps_flg = 0;
};
defined($tbx) and do {
    &copy_files("$tbx", "$SGPTMP.tbx");
    $blast_flg = 0;
};
# $geneid_flg = 0;
defined($ofn) and ($savefiles_flg = 1);
defined($ps_output) or ($plots_flg = 0);
#line 240 "./sgp2.nw"
} # sub Which_Options

#line 75 "./sgp2.nw"
### MAIN SUBS ###

#line 484 "./sgp2.nw"
# Runnig blast
sub Run_Blast() {
    # building DB on SEQ1
    # convert the first fasta file to database
    &prt_Header("Compiling BLAST DB on ${Loc1}") if $verbose_flg;
    my $Seq1DB = "$SGPTMP.$Seq1_Name.DB";   #temporal name for database
    copy("$Seq1", "$Seq1DB")
        or &go_to_die("FATAL ERROR !!! Could not copy file '$Seq1' to '$Seq1DB' : $!\n");
    my $pressdb = "$PressdbProgram $Seq1DB";
    $status = system $pressdb;
    &go_to_die("FATAL ERROR !!! $? $!\n  Unsuccessful PRESSDB command :\n  $pressdb \n")
        unless $status == 0;
    # running blast program: Seq1 as DB, Seq2 as query.
    &prt_Header("Running BLAST on ${Loc1}(DB) : ${Loc2}(Query)") if $verbose_flg;
    my ($dir, $ext);
    ($Seq1DB, $ENV{'BLASTDB'}, $ext) = fileparse($Seq1DB,'');
    ($BlastMatrix, $dir, $ext) = fileparse($BlastMatrix,'');
    ($ENV{'BLASTMAT'} = $dir) =~ s%/$%% ;
    my $blast = "$BlastProgram $Seq1DB $Seq2 -matrix $BlastMatrix $BlastOptions > $SGPTMP.tbx";
    $status = system $blast;
    &go_to_die("FATAL ERROR !!! $? $!\n  Unsuccessful BLAST command :\n  $blast \n")
        unless $status == 0;
    unlink("$ENV{'BLASTDB'}$Seq1DB");
    # saving files
    copy("$SGPTMP.tbx", "$ofn${Loc1}_${Loc2}.tbx") if $savefiles_flg;
} # END_SUB: Run_Blast
#line 537 "./sgp2.nw"
# Extracting and processing HSPs
#line 566 "./sgp2.nw"
sub write_HSP_files() {
    my @names = @_;
    my ($DSC, $SHSP, $scnd, $thrd) = (38, 0, 0, 0);
    open(SR,"< ".(shift @names));
    open(S1,"> $names[0]");
    @names >= 2 and (open(S2,"> $names[1]"), $scnd = 1);
    @names == 3 and (open(S3,"> $names[2]"), $thrd = 1);
    while (<SR>) {
		next if /^\#/;
        next if /^[ \t]*$/;
        chomp;
        my @ln = split;
        $ln[10] =~ s/\"//go;
        print S1 "$ln[1] $ln[2] $ln[3] ".($ln[4]+$SHSP)." ".($ln[5]-$SHSP)." ".($ln[6]-$DSC)." $ln[7] $ln[8]\n";
        print S2 "$ln[10] $ln[2] $ln[3] ".($ln[11]+$SHSP)." ".($ln[12]-$SHSP)." ".($ln[6]-$DSC)." $ln[16] $ln[18]\n"
            if $scnd;
        print S3 "$ln[1]:$ln[10] $ln[2] hsp:hsp ".($ln[4]+$SHSP).":".($ln[11]+$SHSP)." ".($ln[5]-$SHSP).":".($ln[12]-$SHSP)." ".($ln[6]-$DSC)." $ln[7]:$ln[16] $ln[8]:$ln[18]\n"
            if $thrd;
    };
    close(S3) if $thrd;
    close(S2) if $scnd;
    close(S1);
    close(SR);
}
#line 539 "./sgp2.nw"
sub Extract_HSP() {
    &prt_Header("Extracting SRs from BLAST output") if $verbose_flg;
    #
    my $getSR = "$GetSRsProgram -bC $S_CUTOFF -H $SGPtmpG -W $SGPtmpG -- $SGPTMP.tbx";
    $status = system $getSR;
    &go_to_die("FATAL ERROR !!! $? $!\n  Unsuccessful GETSRs command :\n  $getSR \n")
        unless $status == 0;
    #
    &write_HSP_files("$SGPtmpG.srS", "$SGPtmp1.hsp-sr"); # I suposse that $SGPtmp1 is Subject
    &write_HSP_files("$SGPtmpG.srS", "$SGPtmp2.hsp-sr"); # and $SGPtmp2 is query, else modify this lines.
    &write_HSP_files("$SGPtmpG.hsp", "$SGPtmp2.hsp", "$SGPtmp1.hsp", "$SGPTMP.aln"); 
    # 
  SAVEFILES: do {
      copy("$SGPtmpG.srQ",    "$ofn${Loc1}_$Loc2.srQ"); # Those files have the relationships between
      copy("$SGPtmpG.srS",    "$ofn${Loc1}_$Loc2.srS");  # Query and Subject coordinates so you can
      copy("$SGPtmpG.hsp",    "$ofn${Loc1}_$Loc2.hsp");  # recover the sequence substring for both.
      copy("$SGPtmp1.hsp-sr", "$ofn$Loc1.hsp-sr");
      copy("$SGPtmp2.hsp-sr", "$ofn$Loc2.hsp-sr");
      copy("$SGPtmp1.hsp",    "$ofn$Loc1.hsp");
      copy("$SGPtmp2.hsp",    "$ofn$Loc2.hsp");
      copy("$SGPTMP.aln",     "$ofn${Loc1}_$Loc2.aln");
    } if $savefiles_flg;
} # END_SUB: Extract_HSP
#line 636 "./sgp2.nw"
# Running geneid
sub Run_geneid() {
	my $geneid;
    # running geneid in $Seq1 with tblastx output
    &prt_Header("Running GENEID program on $Loc1") if $verbose_flg;
    $geneid = "$GeneidProgram $GeneidOptions -P $GeneidParam -S  $SGPtmp1.hsp-sr $Seq1 >  $SGPtmp1.sgp";
    $status = system $geneid;
    &go_to_die("FATAL ERROR !!! $? $!\nUnsuccessfull GENEID command on '$Loc1':\n  $geneid \n") unless $status == 0;
    # running geneid in $Seq2 with tblastx output
    &prt_Header("Running GENEID program on $Loc2") if $verbose_flg;
    $geneid = "$GeneidProgram $GeneidOptions -P $GeneidParam -S  $SGPtmp2.hsp-sr $Seq2 >  $SGPtmp2.sgp";
    $status = system $geneid;
    &go_to_die("FATAL ERROR !!! $? $!\nUnsuccessfull GENEID command on '$Loc2':\n  $geneid \n") unless $status == 0;
  SAVEFILES: do {
      copy("$SGPtmp1.sgp", "$ofn$Loc1.sgp");
      copy("$SGPtmp2.sgp", "$ofn$Loc2.sgp");
    } if $savefiles_flg;
  ECHOFILES: do {
      system "cat $SGPtmp1.sgp $SGPtmp2.sgp 1>&2"; 
    } if $verbose_flg;
} # END_SUB: Run_geneid
#line 667 "./sgp2.nw"
# Running gff2ps and aplot
sub Make_Plots() {
    my ($gff2ps, $aplot);
    # predicted genes separately -- gff2ps on Loc1
    &prt_Header("Running GFF2PS program on $Loc1") if $verbose_flg;
    $gff2ps = "$Gff2psProgram $SGPtmp1.sgp $SGPtmp1.hsp-sr > $ps_output$Loc1.ps";
    $status = system $gff2ps;
    &go_to_die("FATAL ERROR !!! $? $!\nUnsuccessfull GFF2PS command : $gff2ps \n")
        unless $status == 0;
    # predicted genes separately -- gff2ps on Loc2
    &prt_Header("Running GFF2PS program on $Loc2") if $verbose_flg;
    $gff2ps = "$Gff2psProgram $SGPtmp2.sgp $SGPtmp2.hsp-sr > $ps_output$Loc2.ps";
    $status = system $gff2ps;
    &go_to_die("FATAL ERROR !!! $? $!\nUnsuccessfull GFF2PS command : $gff2ps \n")
        unless $status == 0;
    
#line 686 "./sgp2.nw"
# run aplot
&prt_Header("Running APLOT program on $Loc1:$Loc2") if $verbose_flg;
$aplot = "$AplotProgram ";
$status = system $aplot;
&go_to_die("FATAL ERROR !!! $? $!\nUnsuccessfull APLOT command : $aplot \n")
    unless $status == 0;
#line 683 "./sgp2.nw"
} # END_SUB: Make_Plots

#line 82 "./sgp2.nw"
### MAIN LOOP ###

#line 122 "./sgp2.nw"
# get options
&Which_Options();
# main external program calls
DOBLAST: do {
    &Run_Blast() if $blast_flg;
    &Extract_HSP();
} if $hsps_flg;
&Run_geneid()  if $geneid_flg;
&Make_Plots()  if $plots_flg; 
# removing temporary files
&clean_tmp();
# timing and exit
&get_exec_time(time);
exit(1);

#line 86 "./sgp2.nw"
### EOF ###

#line 822 "./sgp2.nw"
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
