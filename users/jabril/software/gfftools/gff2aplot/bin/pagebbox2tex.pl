#!/usr/local/bin/perl -w
# This is perl, v5.6.1 built for i686-linux
# /usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
# 
# pagebbox2tex.pl
#
# $Id: pagebbox2tex.pl,v 1.1 2002-12-24 13:07:36 jabril Exp $
#
use strict;
#
my $USAGE = "###  pagebbox2tex  <STDINPUT  >STDOUTPUT \n";
#
# MAIN
&parseinput;
exit(0);
#
# SUBS
sub parseinput() {
    my @rec;
    &print_prologue;
    while (<STDIN>) {
        next if /^\s*$/o;
        chomp;
        $_ =~ s/^\s*//o;
        $_ =~ s/["',]//og; #'"
        @rec = split /\s+/og, $_;
        $rec[0] eq '#' && do {
            print STDOUT '\hline\hline'."\n";
            print STDOUT '%                  points    -  centimeters  -     inches'."\n";
            next;
        }; # set color header
        (defined($rec[1]) && $rec[1] eq '=>') && do {
            print STDOUT (&fill_right($rec[0],12," ")).
                         (&get_sizes(@rec[4,5]))." \\\\\n";
        }; # set color line
    }; # while
    &print_trailer;    
} # parseinput
sub print_prologue() {
    print STDOUT << '+++EOP+++';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AplotPageSizeTbl.tex
%
% Page Sizes used in "gff2aplot".
%
% # $Id: pagebbox2tex.pl,v 1.1 2002-12-24 13:07:36 jabril Exp $
%
\label{sec:pagesizes}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill
\begin{table}[!ht]
\begin{center}
\setlength{\fboxsep}{2pt}
%\setlength{\arrayrulewidth}{1pt}
\fbox{
 \begin{tabular}{|c||r|r||r|r||r|r|} \hline
  \raisebox{-0.5ex}[0pt]{ PAGE } &
    \multicolumn{6}{c|}{ PAGE SIZE }\\ \cline{2-7}
  \raisebox{0.25ex}[0pt]{ FORMAT } &
    \multicolumn{2}{c||}{ (in points) } &
    \multicolumn{2}{c||}{ (in cms) } &
    \multicolumn{2}{c|}{ (in inches) } \\
+++EOP+++
} # print_prologue
sub print_trailer() {
    print STDOUT << '+++EOP+++';
  \hline % \hline
 \end{tabular}
} % fbox
\caption{\label{tbl:PageSizes}Page Sizes defined in {\prog}.}\hspace{1cm}
\vskip 1ex
\fbox{
\begin{tabular}{c@{\quad$\equiv$\quad}c}
28.35 pt & 1 cm   \\
72.00 pt & 1 inch \\
\end{tabular}
} % fbox
  %\refstepcounter{table}
  %\addcontentsline{lot}{section}{
  %   \thetable\hspace{1em}Page Sizes available at {\prog}.}
\end{center}
\end{table}
\vfill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
+++EOP+++
} # print_trailer
sub get_sizes() {
    my ($w,$h) = @_;
    my $str = '';
    $str  = ' &'.(&fill_left($w,6," ")).       ' &'.(&fill_left($h,6," "));
    $str .= ' &'.(&fill_left(&tocm($w),6," ")).' &'.(&fill_left(&tocm($h),6," "));
    $str .= ' &'.(&fill_left(&toin($w),6," ")).' &'.(&fill_left(&toin($h),6," "));
    return $str;
} # get_sizes
# to cm:   28.35 pts == 1 cm
sub tocm() { return sprintf("%.1f", ($_[0] / 28.35)); }
# to inch:    72 pts == 1 inch
sub toin() { return sprintf("%.1f", ($_[0] / 72.00)); }
#
sub fill_right() { $_[0].($_[2] x ($_[1] - length($_[0]))) }
sub fill_left()  { ($_[2] x ($_[1] - length($_[0]))).$_[0] }
sub fill_mid()   { 
    my $l = length($_[0]);
    my $k = int(($_[1] - $l)/2);
    ($_[2] x $k).$_[0].($_[2] x ($_[1] - ($l+$k)));
} # fill_mid
