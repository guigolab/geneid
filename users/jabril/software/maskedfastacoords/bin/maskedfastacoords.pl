#!/usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
#
#line 371 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
# 
#  
#line 364 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
#$PROGRAM [options] < fasta_file > coords_file
#line 373 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
#
#  
#line 367 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
#Retrieving masked regions coords from fasta files.
#line 375 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
#
#line 741 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
#
#line 929 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
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
#line 743 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
#
#line 923 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
# $Id: maskedfastacoords.pl,v 1.2 2001-10-03 11:23:38 jabril Exp $
#line 745 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
#
use strict;
#
#line 360 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
my $PROGRAM = 'maskedfastacoords.pl';
my $VERSION = '0.1_alpha';
#line 749 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
my $DATE = localtime;
my $USER = defined($ENV{USER}) ? $ENV{USER} : 'Child Process';
my $host = `hostname`;
chomp($host);
#
#line 341 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
#
# MODULES
#
#line 765 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
use Benchmark;
  
#line 773 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
my @Timer = (new Benchmark);
#line 450 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
use Getopt::Long;
Getopt::Long::Configure qw/ bundling /;
#line 434 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
use Bio::Seq;
use Bio::SeqIO;
#line 345 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
#
# VARIABLES
#
#line 757 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
my ($T,$F) = (1,0); # for 'T'rue and 'F'alse
#line 873 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
my $line = ('#' x 80)."\n";
my $s = '### ';
my $sp = "###\n";
my $Error = "\<\<\<  ERROR  \>\>\> ";
my $Warn  = "\<\<\< WARNING \>\>\> ";
my $spl   = "\<\<\<\-\-\-\-\-\-\-\-\-\>\>\>\n";
my $spw   = "\<\<\<         \>\>\> ";
#line 885 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
my %Messages = (
    # ERROR MESSAGES
    
#line 478 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
'UNKNOWN_CL_OPTION' =>
  $Warn."Error trapped while processing command-line:\n".(" "x16)."\%s\n",
'CMD_LINE_ERROR' =>
  $spl.$spw." Please, check your command-line options!!!\n".$Error."\n".
  $spw." ".("."x12)." Type \"maskedfastacoords.pl -h\" for help.\n".$spl,
#line 888 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
    # WORKING MESSAGES
    
#line 406 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
'PROG-START'  => "$line$s\n$s Running $PROGRAM\n$s".
                 "$s HOST: $host".
                 "$s USER: $USER".
                 "$s DATE: $DATE\n$s\n$line$s" ,
'PROG-FINISH' => "$s\n$line$s\n$s $PROGRAM FINISHED\n$s".
                 "$s TOTAL TIME: \%s\n$line" ,
#line 890 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
    
#line 506 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
'SHOW_VERSION' => $sp."### \%s -- Version: \%s\n".$sp,
#line 891 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
   ); # %Messages 
#line 387 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
my ($id,$seq) = ('','');
my %CmdLineVar = (
                 
#line 589 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
MASKING  => 0,
#line 390 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
                 
#line 608 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
LARGE    => 0,
#line 391 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
                 
#line 626 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
GFF      => 0,
sequence => 'noname',
source   => 'masked',
feature  => 'masked',
strand   => '.',
group    => '',
#line 392 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
                 );
#line 659 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
my $seqin;
#line 349 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
#
# MAIN LOOP
#
#line 396 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
&parse_cmdline(); # PROG-START

&main();

&report('PROG-FINISH',&timing($T));

exit(0);
#line 353 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
#
# FUNCTIONS
#
#line 461 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
sub parse_cmdline() {

    $SIG{__WARN__} = sub { &warn('UNKNOWN_CL_OPTION',$T,$_[0]) };
    GetOptions(
               
#line 584 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
"a|all-masked"   => sub { $CmdLineVar{MASKING} = 2 },
"s|soft-masked"  => sub { $CmdLineVar{MASKING} = 4 },
"m|merge-masked" => sub { $CmdLineVar{MASKING} = 8 },
#line 466 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
               
#line 605 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
"l|large-fasta"  => \$CmdLineVar{LARGE},
#line 467 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
               
#line 618 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
"g|gff"        => \$CmdLineVar{GFF},
"seq-name=s"   => \$CmdLineVar{sequence},
"source=s"     => \$CmdLineVar{source},
"f|feature=s"  => \$CmdLineVar{feature},
"strand=s"     => \$CmdLineVar{strand},
"group=s"      => \$CmdLineVar{group},
#line 468 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
               
#line 490 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
"version"   => \&prt_version, 
"h|help|?"  => \&prt_help,
#line 469 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
               ) || (&warn('CMD_LINE_ERROR',$T), exit(1));
    $SIG{__WARN__} = 'DEFAULT';

    &report("PROG-START");
} # parse_cmdline
#line 499 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
sub prt_version() {
    &report('SHOW_VERSION',$PROGRAM,$VERSION);
    exit(1);
} # prt_version
#line 512 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
sub prt_help() {
    print STDERR <<"+++EndOfHelp+++";
PROGRAM:
                        $PROGRAM - $VERSION

    
#line 367 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
#Retrieving masked regions coords from fasta files.

#line 519 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
USAGE:    
#line 364 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
#$PROGRAM [options] < fasta_file > coords_file


#line 522 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
DESCRIPTION:

    Retrieving positions for masked regions of masked sequences 
    in fasta format. Program can also output those coords in GFF
    format, where you can define several fields from command-line.


REQUIRES:

    
#line 549 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
$PROGRAM needs the following Perl modules 
installed in your system, we used those available 
from the standard Perl distribution. Those that 
are not in the standard distribution are marked 
with an '(*)', in such cases make sure that you 
already have downloaded them from CPAN 
(http://www.perl.com/CPAN) and installed.

  
#line 439 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
"Bio::Seq"
"Bio::SeqIO" - BioPerl modules to handle sequence objects. (*)
               You can download directly from CPAN or 
               from BioPerl web site at "http://bioperl.org/".
#line 558 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
  
#line 455 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
"Getopt::Long" - processing command-line options.
#line 559 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
  
#line 788 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
"Benchmark" - checking and comparing running times of code.


#line 534 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
COMMAND-LINE OPTIONS:

    
#line 563 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
A double dash on itself "--" signals end of the options
and start of file names (if present). After double dash,
you can use a single dash "-" as STDIN placeholder. 
Available options and a short description are listed here:

+ General options:

  
#line 494 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
-h, --help            Shows this help.
--version             Shows current version and exits.

#line 572 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
  
#line 592 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
-a, --all-masked       Default is looking for "N" masked sequence
                       segments. This option distinguish between 
                       upper/lower-case and assumes lower-case to 
                       define masked regions too (so called soft-masking).
-s, --soft-masked      Only outputs soft-masked sequence segments.
-m, --merge-masked     Previous options differentiate classical-masking
                       (with "N"s) from soft-masking. This option merges
                       both as if they were the same masking type
                       (it also enables soft-masking search so previous 
                       options are not required when passing this one).

#line 574 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
  
#line 611 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
-l, --large-fasta      For large genomic sequence you can enable
                       Bio::SeqIO to work with temporary files
                       avoiding memory overload (though it makes
                       the script to run slowly).

#line 576 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
+ GFF output:

  
#line 634 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
-g, --gff              Output in GFF format (default coords).

--seq-name <string>    Sets sequence GFF field (default "noname").
--source <string>      Sets source GFF field (default "masked").
-f, --feature <string> Sets feature GFF field (default "masked").
--strand <string>      Sets strand GFF field (default ".").
--group  <string>      Sets group GFF field (default none).


#line 539 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
BUGS:    Report any problem to 'jabril\@imim.es'.

AUTHOR:  $PROGRAM is under GNU-GPL (C) 2000 - Josep F. Abril

+++EndOfHelp+++
    exit(1);
} # prt_help
#line 646 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
sub main() {
    $seqin = Bio::SeqIO->new(-format => 'FASTA', -fh => \*STDIN);
    while (my $sequence = $seqin->next_seq()) {
        my ($sid,$len,$seq,@nuc,@coords,$masked_flg,$match,$msk_num);
        @coords = ();
        
#line 663 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
print STDERR "### READING FASTA............\n";
$sid  = $sequence->display_id();
$len  = $sequence->length();
$seq  = $sequence->seq();
#line 652 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
        
#line 679 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
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
#line 653 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
        
#line 702 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
$msk_num = scalar(@coords) / 2;
print STDERR "###         WRITING GFF COORDS: $msk_num masked regions found.\n";
for (my $n = 0; $n <= $#coords; $n+=2) {
    my $GFFstring = ("\%s\t" x 5).(".\t" x 3).".\n";
    printf STDOUT $GFFstring, $sid, "masked", "masked", @coords[$n..($n + 1)];
}; # for coords in @coords
#line 654 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
    }; # while 
} # main
#line 777 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
sub timing() {
    push @Timer, (new Benchmark);
    # partial time 
    $_[0] || 
        (return timestr(timediff($Timer[$#Timer],$Timer[($#Timer - 1)])));
    # total time
    return timestr(timediff($Timer[$#Timer],$Timer[0]));
} # timing
#line 861 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
sub report() { print STDERR sprintf($Messages{ shift @_ },@_) }
#line 867 "/home/ug/jabril/development/softjabril/maskedfastacoords/maskedfastacoords.nw"
sub warn() { print STDERR sprintf($Messages{ shift @_ }, @_) }
