#!/usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
#
#line 277 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
# #----------------------------------------------------------------#
# #                          fullgffcoords                         #
# #----------------------------------------------------------------#
# 
#  fullgffcoords [options] < input_files > output_files
#
#  Retrieving CDS and protein coords from GFF mapped on genomic coords.
#
#
#     Copyright (C) 2001 - Josep Francesc ABRIL FERRANDO  
#line 903 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
#
#line 1095 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
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
#line 905 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
#
#line 1089 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
# $Id: fullgffcoords.pl,v 1.1 2001-10-22 14:44:13 jabril Exp $
#line 907 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
#
use strict;
#
#line 266 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
my $PROGRAM = 'fullgffcoords.pl';
my $VERSION = '0.1_alpha';
#line 911 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
my $DATE = localtime;
my $USER = defined($ENV{USER}) ? $ENV{USER} : 'Child Process';
my $host = `hostname`;
chomp($host);
#
#line 293 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
#
# MODULES
#
#line 928 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
use Benchmark;
  
#line 936 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
my @Timer = (new Benchmark);
#line 367 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
use Getopt::Long;
Getopt::Long::Configure qw/ bundling /;
#line 297 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
#
# VARIABLES
#
#line 919 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
my ($T,$F) = (1,0); # for 'T'rue and 'F'alse
#line 1015 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
my ($n,$c); # counter and char (for &counter function)
#line 1037 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
my $line = ('#' x 80)."\n";
my $s = '### ';
my $sp = "###\n";
my $Error = "\<\<\<  ERROR  \>\>\> ";
my $Warn  = "\<\<\< WARNING \>\>\> ";
my $spl   = "\<\<\<\-\-\-\-\-\-\-\-\-\>\>\>\n";
my $spw   = "\<\<\<         \>\>\> ";
#line 1049 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
my %Messages = (
    # ERROR MESSAGES
    
#line 398 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
'UNKNOWN_CL_OPTION' =>
  $Warn."Error trapped while processing command-line:\n".(" "x16)."\%s\n",
'CMD_LINE_ERROR' =>
  $spl.$spw." Please, check your command-line options!!!\n".$Error."\n".
  $spw." ".("."x12)." Type \"$PROGRAM -h\" for help.\n".$spl,
#line 1052 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
    
#line 449 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
FILE_NO_OPEN =>
  $spl.$Warn."Cannot Open Current file \"\%s\" . Not used !!!\n".$spl,
#line 726 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
SWAPPINGCOORDS =>
  $Warn."START greater than END (\%s > \%s). SWAPPING COORDS !!!\n",
NOT_EQ_STRAND =>
  $Warn."GROUP STRAND (\%s) DOES NOT MATCH GFF-RECORD (\%s). GROUP RULES !!!\n",
#line 1053 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
    # WORKING MESSAGES
    
#line 334 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
'PROG-START'  => "$line$s\n$s Running $PROGRAM\n$s\n".
                 "$s HOST: $host\n".
                 "$s USER: $USER\n".
                 "$s DATE: $DATE\n$s\n$line$s\n" ,
'PROG-FINISH' => "$s\n$line$s\n$s $PROGRAM FINISHED\n$s\n".
                 "$s TOTAL TIME: \%s\n$line" ,
#line 654 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
SETTING_VARS => $sp."### Variable Definition Finished...\n".$sp,
#line 723 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
PARSING_GFF => $sp."### PARSING GFF FILE: \%s\n".$sp,
#line 779 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
SORTING_FEATURES => $sp."### Sorting GFF Group Features\n".$sp,
SORTING_GROUPS   => $sp."### Sorting GFF Groups\n".$sp,
#line 827 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
WRITING_FEATURES => $sp."### Writing NEW GFF Records to STDOUT...\n".$sp,
#line 1055 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
    
#line 496 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
'SHOW_VERSION' => $sp."### \%s -- Version: \%s\n".$sp,
#line 1056 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
    
#line 454 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
CHECKING_FILENAMES =>
  $sp."### Validating INPUT FILENAMES\n".$sp,
READING_FILE =>
  "###---> \"\%s\" exists, including as Input File.\n",
READING_STDIN =>
  "###---> Including GFF records from standard input.\n",  
#line 1057 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
   ); # %Messages 
#line 320 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
my ($id,$seq) = ('','');
my %CmdLineVar = (
                 
#line 575 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
GFFCOORDS  => 1,
SHOWGROUPS => 0,
#line 323 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
                 
#line 595 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
GFFVERSION => 1,
#line 324 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
                 );
#line 445 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
my @data_files = ();
#line 651 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
my ($GFF,$printGFF);
#line 720 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
my %GeneList = ();
#line 776 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
my @SortedGenes = ();
#line 301 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
#
# MAIN LOOP
#
#line 328 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
&main();

exit(0);
#line 305 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
#
# FUNCTIONS
#
#line 378 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub parse_cmdline() {
    
#line 408 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
my $cmdln_stdin = undef;
for (my $a = 0; $a <= $#ARGV; $a++) { 
    next unless $ARGV[$a] =~ /^-$/o;
    $cmdln_stdin = $a - $#ARGV;
    splice(@ARGV,$a,1);
};    

#line 381 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
    $SIG{__WARN__} = sub { &warn('UNKNOWN_CL_OPTION',$T,$_[0]) };
    GetOptions(
               
#line 569 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
"g|genomic-coords" => sub { $CmdLineVar{GFFCOORDS} = 1 },
"c|cds-coords"     => sub { $CmdLineVar{GFFCOORDS} = 2 },
"p|protein-coords" => sub { $CmdLineVar{GFFCOORDS} = 3 },
"i|include-groups" => \$CmdLineVar{SHOWGROUPS},
#line 384 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
               
#line 591 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
"1|gff-version1"    => sub { $CmdLineVar{GFFVERSION} = 1 },
"2|gff-version2"    => sub { $CmdLineVar{GFFVERSION} = 2 },
#line 385 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
               
#line 480 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
"version"   => \&prt_version, 
"h|help|?"  => \&prt_help,
#line 386 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
               ) || (&warn('CMD_LINE_ERROR',$T), exit(1));
    $SIG{__WARN__} = 'DEFAULT';

    &report("PROG-START");
    @data_files = ();
    &set_input_file($cmdln_stdin);
    @ARGV = (); # ensuring that command-line ARGVs array is empty
    
} # parse_cmdline
#line 419 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub set_input_file() {
    my $stdin_flg = $F;
    
#line 465 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
my $chk_stdin = shift @_;
my $t = scalar(@ARGV);
defined($chk_stdin) && do {
    abs($chk_stdin) > $t && ($chk_stdin = -$t);
	$chk_stdin > 0  && ($chk_stdin = 0 );
    $t += $chk_stdin;
    splice(@ARGV,$t,0,'-');
};
#line 422 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
    &report("CHECKING_FILENAMES");
  FILECHK: foreach my $test_file (@ARGV) {
        $test_file ne '-' && do {
            -e $test_file || do {
                &warn('FILE_NO_OPEN',$T,$test_file);
                next FILECHK;
            };
            &report('READING_FILE',$test_file);
            push @data_files, $test_file;
            next FILECHK;
        };
        $stdin_flg = $T;
        push @data_files, '-';
	}; # foreach
    scalar(@data_files) == 0 && do {
        push @data_files, '-';
        $stdin_flg = $T;
    };
    $stdin_flg && &report('READING_STDIN');
} # set_input_file
#line 489 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub prt_version() {
    &report('SHOW_VERSION',$PROGRAM,$VERSION);
    exit(1);
} # prt_version
#line 502 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub prt_help() {
    print STDERR <<"+++EndOfHelp+++";
PROGRAM:
                        $PROGRAM - $VERSION

    
#line 273 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
Retrieving CDS and protein coords from GFF mapped on genomic coords.

#line 509 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
USAGE:    
#line 270 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
$PROGRAM [options] < input_files > output_files


#line 512 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
DESCRIPTION:

    
#line 273 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
Retrieving CDS and protein coords from GFF mapped on genomic coords.


#line 517 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
REQUIRES:

    
#line 537 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
$PROGRAM needs the following Perl modules 
installed in your system, we used those available 
from the standard Perl distribution. Those that 
are not in the standard distribution are marked 
with an '(*)', in such cases make sure that you 
already have downloaded them from CPAN 
(http://www.perl.com/CPAN) and installed.

  
#line 372 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
"Getopt::Long" - processing command-line options.
#line 546 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
  
#line 951 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
"Benchmark" - checking and comparing running times of code.


#line 522 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
COMMAND-LINE OPTIONS:

    
#line 550 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
A double dash on itself "--" signals end of the options
and start of file names (if present). After double dash,
you can use a single dash "-" as STDIN placeholder. 
Available options and a short description are listed here:

+ General options:

  
#line 484 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
-h, --help            Shows this help.
--version             Shows current version and exits.

#line 559 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
  
#line 579 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
-g, --genomic-coords   Output GFF coords fields (columns 4th and 5th)
                       are set as genomic coords (by default).
-c, --cds-coords       Output GFF coords fields (columns 4th and 5th)
                       are set to CDS coords.
-p, --protein-coords   Output GFF coords fields (columns 4th and 5th)
                       are set to protein coords (using partial codon
                       notation as explained in ... ).
-i, --include-groups   Output a GFF record with the start/end coords
                       for each group in the input GFF file.

#line 561 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
+ GFF output:

  
#line 598 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
-2, --gff-version2
-1, --gff-version1     Output GFF version, default is 1, where version 1
                       means simple grouping field and version 2 forces
                       grouping to be tag-value pair grouping field, 
                       putting value between double quotes.


#line 527 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
BUGS:    Report any problem to 'jabril\@imim.es'.

AUTHOR:  $PROGRAM is under GNU-GPL (C) 2000 - Josep F. Abril

+++EndOfHelp+++
    exit(1);
} # prt_help
#line 608 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub main() {
    &parse_cmdline(); # PROG-START
    &set_output_subs();
    foreach my $lfile (@data_files) {
        &parse_input_files($lfile);
    }; # foreach $lfile
    &sort_by_acceptor();
    &output_extended_GFF();
    &report('PROG-FINISH',&timing($T));
} # main
#line 628 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub set_output_subs() {
    &report("SETTING_VARS");
    my $base = "\%s\t" x 8 ;
    if ($CmdLineVar{GFFVERSION} == 1) {
        $GFF = $base."\%s  \%s  \%s\n" ; # base + group coordsA coordsB
    } else {                        # base + tag "group"; coordsA; coordsB
        $GFF = $base.'gene_id "%s"; %s; %s'."\n" ;
    };
    CHECKCOORDS: {
        $CmdLineVar{GFFCOORDS} == 3 && do {
            $printGFF = \&gff_protein;
            last CHECKCOORDS;
        };
        $CmdLineVar{GFFCOORDS} == 2 && do {
            $printGFF = \&gff_cds;
            last CHECKCOORDS;
        };
        $printGFF = \&gff_genomic;
    };
} # set_output_subs
#line 660 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub gff_genomic() {
    my @data = @_;
    printf STDOUT $GFF, @data[0..8],
                        "CDS @data[9..10]", "AA @data[11..12]";
} # gff_genomic
sub gff_cds() {
    my @data = @_;
    printf STDOUT $GFF, @data[0..2,9..10,5..8],
                        "BASE @data[3..4]", "AA @data[11..12]";
} # gff_cds
sub gff_protein() {
    my @data = @_;
    printf STDOUT $GFF, @data[0..2,11..12,5..8],
                        "BASE @data[3..4]", "CDS @data[9..10]";
} # gff_protein
#line 681 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub parse_input_files() {
    my $ifile = $_[0];
    &report("PARSING_GFF",$ifile);
    open(IFILE,"< $ifile") || do {
        &warn('FILE_NO_OPEN',$T,$ifile);
        return;
    };
    ($n,$c) = (0,undef);
    while (<IFILE>) {
        my (@f,$start,$end,$group,$strand);
        $c = '*';
        
#line 970 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
next if /^\#/o;
next if /^\s*$/o;
chomp;
#line 693 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
        @f = split /\s+/o, $_, 9;
        ($group,$c) = &checkgroup($f[8]);  
        ($start,$end,$strand) = (@f[3,4],"$f[6]");
        ($start > $end) && do {
            &warn('SWAPPINGCOORDS',$T,$start,$end);
            ($start,$end) = ($end,$start);
        };
        defined($GeneList{$group}{STRAND}) || 
            ($GeneList{$group}{STRAND} = $strand);
        ($GeneList{$group}{STRAND} ne $strand) && do {
            &warn('NOT_EQ_STRAND',$T,$GeneList{$group}{STRAND},$strand);
            # does nothing, just warns.
        };
        push @{ $GeneList{$group}{RECORDS} }, [ @f[0..7] ];
        push @{ $GeneList{$group}{MIN} }, $start;
        push @{ $GeneList{$group}{MAX} }, $end;
        defined($GeneList{$group}{LENGTH}) || ($GeneList{$group}{LENGTH} = 0);
        $GeneList{$group}{LENGTH} += ($end - $start + 1);
    } continue {
        &counter(++$n,$c);
    }; # while
    &counter_end($n,$c);
    close(IFILE);
} # parse_input_files
#line 733 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub checkgroup() {
    my ($gpA,$gpB);
    my ($gpstr) = @_;
    $gpstr =~ /^([^\s]+)(?:\s+"(.+?)")?(?:.*)?/o && do {
        $gpB = defined($2) ? $2 : ''; # for GFF version 2
        $gpA = defined($1) ? $1 : $gpB;  # for GFF version 1
    };
    ($gpB ne '') && ( return $gpB, ':' );
    return $gpA, '.';
} # checkgroups
#line 749 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub sort_by_acceptor() {
    &report("SORTING_FEATURES");
    $c = 0;
    foreach my $gname (keys %GeneList) {
        my $ref = \@{ $GeneList{$gname}{RECORDS} };
        print STDERR "$gname...(".($#{$ref} + 1)."ft)".
                     (((++$c % 4) == 0) ? "\n" : "\t");
        # sorting all features by acceptor
        @{ $ref } = map  { $_->[2] }
                    sort { &sort_forward }
                    map  { [ $_->[3], $_->[4], $_ ] } @{ $ref };
        # getting group boundaries
        $GeneList{$gname}{MIN} = min(@{ $GeneList{$gname}{MIN} });
        $GeneList{$gname}{MAX} = max(@{ $GeneList{$gname}{MAX} });
        push @SortedGenes,
             [ $gname, $GeneList{$gname}{MIN}, $GeneList{$gname}{MAX},
               $GeneList{$gname}{STRAND}, $GeneList{$gname}{LENGTH} ];
    }; # foreach $gname
    ((++$c % 4) != 0) && print STDERR "\n";
    &report("SORTING_GROUPS");
    @SortedGenes = map  { $_->[2] }
                   sort { &sort_forward }
                   map  { [ $_->[1], $_->[2], $_ ] } @SortedGenes;
} # sort_by_acceptor
#line 784 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub sort_forward {
    $a->[0] <=> $b->[0]  # sorting by start
             or
    $b->[1] <=> $a->[1]; # reverse sorting by end if same start
} # sort_forward
#line 795 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub output_extended_GFF() {
    my $rpt = '%-15s%12s%12s%7s%10s'."\n";
    &report("WRITING_FEATURES");
    printf STDERR $rpt,"# GeneName","MinCoord","MaxCoord","Strand","CDSlength";
    for (my $g = 0; $g <= $#SortedGenes; $g++) {
        my ($ref,$gene,$ori,$end,$str,$len,$f,$do_it,$base);
        ($gene,$ori,$end,$str,$len) = @{ $SortedGenes[$g] };
        printf STDERR $rpt,"  $gene",$ori,$end,$str,$len;
        $ref = \@{ $GeneList{$gene}{RECORDS} };
        if ($str ne '-') {
            $do_it = \&do_it_forward;
            $base  = 0;
        } else {	 
            $do_it = \&do_it_reverse;
            $base  = $len + 1;
        };
        for ($f = 0; $f <= $#{$ref}; $f++) {
            my ($o_cds,$e_cds,$o_aa,$e_aa);
            # ($o_gn,$e_gn) = ($ref->[],$ref->[]);
            ($o_cds,$e_cds) = &$do_it($base,$ref->[$f][3],$ref->[$f][4]);
            $base = $e_cds;
            ($o_aa,$e_aa) = &to_protein($o_cds,$e_cds);
            ($str eq '-') && do { # just to display GFF-like (start<end)
                ($o_cds,$e_cds,$o_aa,$e_aa) = ($e_cds,$o_cds,$e_aa,$o_aa);
            };
            &$printGFF(@{ $ref->[$f] },$gene,$o_cds,$e_cds,$o_aa,$e_aa);
        }; # for $f
    }; # for $g
} # output_extended_GFF
#line 831 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub do_it_forward() {
    my ($bori,$gori,$gend) = @_;
    my ($cori,$cend);
    $cori = $bori + 1;
    $cend = $bori + ($gend - $gori + 1);
    return ($cori,$cend);
} # do_it_forward
#line 843 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub do_it_reverse() {
    my ($bori,$gori,$gend) = @_;
    my ($cori,$cend);
    $cori = $bori - 1;
    $cend = $cori - ($gend - $gori);
    return ($cori,$cend);
} # do_it_reverse
#line 855 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub to_protein() {
    my ($cori,$cend) = @_;
    my ($pori,$pend);
    $pori = &toprot($cori);
    $pend = &toprot($cend);
    return ($pori,$pend);
} # to_protein
sub toprot() {
    my ($val) = @_;
    my ($a,$b);
    $a = $val % 3;                # nucleotide order within codon
    ($a == 0) && ($a = 3);
    $b = int(($val - 1) / 3) + 1; # codon number
    return "$b.$a"; 
} # toprot
#line 990 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
#
sub fill_right() { $_[0].($_[2] x ($_[1] - length($_[0]))) }
sub fill_left()  { ($_[2] x ($_[1] - length($_[0]))).$_[0] }
sub fill_mid()   { 
    my $l = length($_[0]);
    my $k = int(($_[1] - $l)/2);
    ($_[2] x $k).$_[0].($_[2] x ($_[1] - ($l+$k)));
} # fill_mid
#line 1003 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
#
sub counter { # $_[0]~current_pos++ $_[1]~char
    print STDERR "$_[1]";
    (($_[0] % 50) == 0) && (print STDERR "[".&fill_left($_[0],6,"0")."]\n");
} # counter
#
sub counter_end { # $_[0]~current_pos   $_[1]~char
    (($_[0] % 50) != 0) && (print STDERR "[".&fill_left($_[0],6,"0")."]\n");
} # counter_end
#line 976 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
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
#line 940 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub timing() {
    push @Timer, (new Benchmark);
    # partial time 
    $_[0] || 
        (return timestr(timediff($Timer[$#Timer],$Timer[($#Timer - 1)])));
    # total time
    return timestr(timediff($Timer[$#Timer],$Timer[0]));
} # timing
#line 1025 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub report() { print STDERR sprintf($Messages{ shift @_ },@_) }
#line 1031 "/home/ug/jabril/development/softjabril/fullgffcoords/fullgffcoords.nw"
sub warn() { print STDERR sprintf($Messages{ shift @_ }, @_) }
