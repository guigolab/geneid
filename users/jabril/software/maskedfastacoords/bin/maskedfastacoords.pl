# This is perl, version 5.005_03 built for i386-linux
#
# 
#  #$PROGRAM [options] < fasta_file > coords_file
#
#  #Retrieving masked regions coords from fasta files.
#
#
# #----------------------------------------------------------------#
# #                          maskedfastacoords                         #
# #----------------------------------------------------------------#
# 
#    Remember to put a short description of your script here...
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
#
# $Id: maskedfastacoords.pl,v 1.1 2001-09-25 17:34:08 jabril Exp $
#
use strict;
#
# BEGIN {
    my $PROGRAM = 'getfastamasked.pl';
    my $VERSION = '0.1_alpha';
    my $DATE = localtime;
    my $USER = defined($ENV{USER}) ? $ENV{USER} : 'Child Process';
    my $host = `hostname`;
    chomp($host);
# } # BEGIN
#
# MODULES
#
use Benchmark;
  my @Timer = (new Benchmark);
use Getopt::Long;
Getopt::Long::Configure qw/ bundling /;
use Bio::Seq;
use Bio::SeqIO;
#
# VARIABLES
#
my ($T,$F) = (1,0); # for 'T'rue and 'F'alse
my $line = ('#' x 80)."\n";
my $s = '### ';
my $Error = "\<\<\<  ERROR  \>\>\> ";
my $Warn  = "\<\<\< WARNING \>\>\> ";
my $spl   = "\<\<\<\-\-\-\-\-\-\-\-\-\>\>\>\n";
my $spw   = "\<\<\<         \>\>\> ";
my %Messages = (
    # ERROR MESSAGES
    'UNKNOWN_CL_OPTION' =>
      $Warn."Error trapped while processing command-line:\n".(" "x16)."\%s\n",
    'CMD_LINE_ERROR' =>
      $spl.$spw." Please, check your command-line options!!!\n".$Error."\n".
      $spw." ".("."x12)." Type \"maskedfastacoords.pl -h\" for help.\n".$spl,
    # WORKING MESSAGES
    'PROG-START'  => "$line$s\n$s Running $PROGRAM\n$s".
                     "$s HOST: $host".
                     "$s USER: $USER".
                     "$s DATE: $DATE\n$s\n$line$s" ,
    'PROG-FINISH' => "$s\n$line$s\n$s $PROGRAM FINISHED\n$s".
                     "$s TOTAL TIME: \%s\n$line" ,
   ); # %Messages 
my ($id,$seq) = ('','');
#
# MAIN LOOP
#
&parse_cmdline(); # PROG-START

&main();

&report('PROG-FINISH',&timing($T));

exit(0);
#
# FUNCTIONS
#
sub parse_cmdline() {

    $SIG{__WARN__} = sub { &warn('UNKNOWN_CL_OPTION',$T,$_[0]) };
    GetOptions(
               "g|gff"        => ,
               "seq-name=s"   => ,
               "source=s"     => ,
               "f|feature=s"  => ,
               "strand=s"     => ,
               "group=s"      => ,
               "version"   => \&prt_version, 
               "h|help|?"  => \&prt_help,
               ) || (&warn('CMD_LINE_ERROR',$T), exit(1));
    $SIG{__WARN__} = 'DEFAULT';

    &report("PROG-START");
} # parse_cmdline
sub prt_version() {
    &report('SHOW_VERSION',$PROGRAM,$VERSION);
    exit(1);
} # prt_version
sub prt_help() {
    print STDERR <<"+++EndOfHelp+++";
PROGRAM:
                        $PROGRAM - $VERSION

    #Retrieving masked regions coords from fasta files.

USAGE:    #$PROGRAM [options] < fasta_file > coords_file


DESCRIPTION:

    Retrieving positions for masked regions of masked sequences 
    in fasta format. Program can also output those coords in GFF
    format, where you can define several fields from command-line.


REQUIRES:

    $PROGRAM needs the following Perl modules 
    installed in your system, we used those available 
    from the standard Perl distribution. Those that 
    are not in the standard distribution are marked 
    with an '(*)', in such cases make sure that you 
    already have downloaded them from CPAN 
    (http://www.perl.com/CPAN) and installed.

      "Bio::Seq"
      "Bio::SeqIO" - BioPerl modules to handle sequence objects. (*)
                     You can download directly from CPAN or 
                     from BioPerl web site at "http://bioperl.org/".
      "Getopt::Long" - processing command-line options.
      "Benchmark" - checking and comparing running times of code.


COMMAND-LINE OPTIONS:

    A double dash on itself "--" signals end of the options
    and start of file names (if present). After double dash,
    you can use a single dash "-" as STDIN placeholder. 
    Available options and a short description are listed here:

    + General options:

      -h, --help            Shows this help.
      --version             Shows current version and exits.  

    + GFF output:

      -g, --gff              Output in GFF format (default coords).
      --seq-name <string>    Sets sequence GFF field (default "noname").
      --source <string>      Sets source GFF field (default "masked").
      -f, --feature <string> Sets feature GFF field (default "masked").
      --strand <string>      Sets strand GFF field (default ".").
      --group  <string>      Sets group GFF field (default none).


BUGS:    Report any problem to 'jabril\@imim.es'.

AUTHOR:  $PROGRAM is under GNU-GPL (C) 2000 - Josep F. Abril

+++EndOfHelp+++
    exit(1);
} # prt_help
sub main() {
    my $seqin = Bio::SeqIO->new(-format => 'FASTA', -fh => \*STDIN);
    while (my $sequence = $seqin->next_seq()) {
        my ($sid,$len,$seq,@nuc,@coords,$masked_flg,$match,$msk_num);
        @coords = ();
        print STDERR "### READING FASTA............\n";
        $sid  = $sequence->display_id();
        $len  = $sequence->length();
        $seq  = $sequence->seq();
        print STDERR "###         PARSING SEQUENCE: $sid ($len bp)\n";
        @nuc = split //og, $seq;
        ($masked_flg,$match) = ($F,$F) ;
        for (my $n = 0; $n <= $#nuc; $n++) {
            $match = ( $nuc[$n] =~ /[NnXx]/o ) ? $T : $F;
            ( !$masked_flg && $match ) && do {
                $masked_flg = $T ;
                # $n contains the last non-masked nucleotide
                push @coords, ($n + 1);
                next;
            };
            $masked_flg && do {
                $match && (next);
                $masked_flg = $F ;
                # $n contains the last masked nucleotide now
                push @coords, $n;
            };
        }; # for nuc in $seq
        # if last nucleotide is masked, previous loop not includes its coord. 
        $masked_flg &&( push @coords, $len);
        $msk_num = scalar(@coords) / 2;
        print STDERR "###         WRITING GFF COORDS: $msk_num masked regions found.\n";
        for (my $n = 0; $n <= $#coords; $n+=2) {
            my $GFFstring = ("\%s\t" x 5).(".\t" x 3).".\n";
            printf STDOUT $GFFstring, $sid, "masked", "masked", @coords[$n..($n + 1)];
        }; # for coords in @coords
    }; # while 
} # main
sub timing() {
    push @Timer, (new Benchmark);
    # partial time 
    $_[0] || 
        (return timestr(timediff($Timer[$#Timer],$Timer[($#Timer - 1)])));
    # total time
    return timestr(timediff($Timer[$#Timer],$Timer[0]));
} # timing
sub report() { print STDERR sprintf($Messages{ shift @_ },@_) }
sub warn() { print STDERR sprintf($Messages{ shift @_ }, @_) }
