#!/usr/local/bin/perl -w
# This is perl, v5.6.1 built for i686-linux
# /usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
# 
# colors2tex.pl
#
# $Id: colors2tex.pl,v 1.2 2003-03-04 18:05:36 jabril Exp $
#
use strict;
use Getopt::Std;
#
my $USAGE = << "+++EOU+++";
################################################################################
###
#### $0
###
####\t\t@{[ (localtime)." : ".(defined($ENV{USER})?$ENV{USER}:"nouser") ]} 
###
#### USAGE:
###
###    colors2tex.pl [options]  < STDIN  > STDOUT
###
###    -h          prints this help.
###    -d          output is set to colors CMYK definition.
###    -t "color"  output is set to colors LaTeX table.
###                In this mode, a color name (defined on input)
###                is required to choose which color starts 
###                a new column in final LaTeX table.
###                "color" builds a regexp like "/^(color)$/"
###
################################################################################
+++EOU+++
#
my $tblflg = 0;
my $splitcolor = '';
my $splitcolornum = 1;
#
# MAIN
&getcmdlineopts;
&parseinput;

exit(0);
#
# SUBS
sub getcmdlineopts() {
    our($opt_d,$opt_t,$opt_h);
    getopts('dt:h');
    $opt_h && do {
        print STDERR $USAGE;
        exit(1);
    }; # $opt_h
    $opt_d && ($tblflg = 0);
    defined($opt_t) && do {
        $splitcolor = $opt_t;
        $splitcolornum = &max($splitcolornum,
                              scalar (split /\|/og, $splitcolor) );
        $tblflg = 1;
    };
} # getcmdlineopts
sub parseinput() {
    my @rec;
    &print_prologue($tblflg);
    while (<STDIN>) {
        next if /^\s*$/o;
        chomp;
        $_ =~ s/^\s*//o;
        @rec = split /\s+/og, $_;
        $rec[0] eq '#' && do {
            defined($rec[1]) || next;
            &print_color_name($tblflg, $rec[1]);
            next;
        }; # set color header
        (defined($rec[1]) && $rec[1] eq '=>') && do {
            &print_color_row($tblflg, @rec[0,5..8]);
        }; # set color line
    }; # while
    &print_trailer($tblflg);
} # parseinput
sub print_prologue() {
    my ($flg) = @_;  
    $flg || do {
        print STDOUT << '+++EOP+++';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AplotColorDefs.tex
%
% Color CMYK definition used in "gff2aplot".
%
% # $Id: colors2tex.pl,v 1.2 2003-03-04 18:05:36 jabril Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
+++EOP+++
        return;
    }; # !$flg
    print STDOUT << '+++EOP+++';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AplotColorTbl.tex
%
% Colors used in "gff2aplot": CMYK values table.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
\begin{footnotesize}
+++EOP+++
    $splitcolornum = 'c'.('@{\quad}c' x $splitcolornum);
    print STDOUT '\begin{tabular}{'.$splitcolornum.'}'."\n".
                 '  \begin{tabular}{|c|c|cccc|} \hline'."\n";
    return;
} # print_prologue
sub print_color_name() {
    my ($flg,$name) = @_;    
    $flg || do {
        print STDOUT '% '."$name\n";
        return;
    }; # !$flg
    $name =~ /^($splitcolor)$/ && &print_new_col;
    print STDOUT '% '."$name\n".'\clspc'."\n";
    return;
} # print_color_name
sub print_new_col() {
    print STDOUT << '+++EOP+++';
      \clspc
      \hline
    \end{tabular} 
   &
    \begin{tabular}{|c|c|cccc|} \hline
+++EOP+++
} # print_new_col
sub print_color_row() {
    my ($flg,$name,$cyan,$magenta,$yellow,$black) = @_; # CMYK   
    $flg || do {
        print STDOUT '\definecolor{'.$name.'}'.(" " x (20 - length($name))).
                     '{cmyk}{'."$cyan,$magenta,$yellow,$black".'}'."\n";
        return;
    }; # !$flg
    print STDOUT '\clrow{'.$name.'}'.(" " x (24 - length($name))).
                 "& $cyan & $magenta & $yellow & $black \\\\\n";
    return;
} # print_color_row
sub print_trailer() {
    my ($flg) = @_;    
    $flg || do {
        print STDOUT << '+++EOP+++';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
+++EOP+++
        return;
    }; # !$flg
    print STDOUT << '+++EOP+++';   
      \clspc
      \hline
    \end{tabular}
   \\
  \end{tabular}
\end{footnotesize}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
+++EOP+++
    return;
} # print_trailer
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
