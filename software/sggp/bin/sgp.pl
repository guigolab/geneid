#!/usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
#
#line 450 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
# $Id: sgp.pl,v 1.1.1.1 2001-06-02 11:24:16 jabril Exp $
#line 466 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
#
use strict;
#
#line 249 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
######################################################################
#                               SGP                                  #
######################################################################
#
#     Sinteny based Gene Prediction tool.
#
#     Copyright (C) 2000 - Josep Francesc ABRIL FERRANDO
#                                   Genis PARRA FARRE
#                                 Roderic GUIGO SERRA
#line 426 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
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

#line 227 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
#### CONSTANTS DEFINITION ####


#### VARIABLES DECLARATION ####


#### MAIN LOOP ####

  
#line 264 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
exit(0);

#line 237 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
##### SUBS #####








#line 241 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
#### EOF ####

#line 320 "/home/ug/jabril/development/sggp/SintenicGenePrediction.nw"
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
