#!/usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
#
#line 277 "/home/ug/jabril/development/softjabril/blat2gff/blat2gff.nw"
# #----------------------------------------------------------------#
# #                          blat2gff.pl                           #
# #----------------------------------------------------------------#
# 
#    blat2gff.pl [options] < inputfile > outputfile
#
#    Remember to put a short description of your script here...
# 
#     Copyright (C) 2001 - Josep Francesc ABRIL FERRANDO  
#line 488 "/home/ug/jabril/development/softjabril/blat2gff/blat2gff.nw"
#
#line 650 "/home/ug/jabril/development/softjabril/blat2gff/blat2gff.nw"
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
#line 490 "/home/ug/jabril/development/softjabril/blat2gff/blat2gff.nw"
#
#line 644 "/home/ug/jabril/development/softjabril/blat2gff/blat2gff.nw"
# $Id: blat2gff.pl,v 1.1 2001-10-22 14:44:13 jabril Exp $
#line 492 "/home/ug/jabril/development/softjabril/blat2gff/blat2gff.nw"
#
use strict;
#
#line 272 "/home/ug/jabril/development/softjabril/blat2gff/blat2gff.nw"
my $PROGRAM = 'blat2gff.pl';
my $VERSION = '0.1_alpha';
#line 496 "/home/ug/jabril/development/softjabril/blat2gff/blat2gff.nw"
my $DATE = localtime;
my $USER = defined($ENV{USER}) ? $ENV{USER} : 'Child Process';
my $host = `hostname`;
chomp($host);
#
#line 292 "/home/ug/jabril/development/softjabril/blat2gff/blat2gff.nw"
#
# MODULES
#

#
# VARIABLES
#
#line 314 "/home/ug/jabril/development/softjabril/blat2gff/blat2gff.nw"
my @frame = (3,1,2);
my %blat = (
    STRAND =>  8,
    Q_NAME =>  9,
    Q_LEN  => 10,
    Q_ORI  => 11,
    Q_END  => 12,
    T_NAME => 13,
    T_LEN  => 14,
    T_ORI  => 15,
    T_END  => 16,
    B_NUM  => 17,
    B_LEN  => 18,
    B_QPOS => 19,
    B_TPOS => 20,
    );
my %gff = ();
my $maxfields = 21;
my %group = ();
my ($group_cnt, $rec_cnt) = (0) x 2;
#line 409 "/home/ug/jabril/development/softjabril/blat2gff/blat2gff.nw"
# GetSRsAln.pl like:
#$srs->{QUERY}\t$blst_prg\tsr\t$srs->{START_Q}\t$srs->{END_Q}
#\t$srs->{SCORE}\t$srs->{STRAND_Q}\t$srs->{FRAME_Q}
#\tTarget \"$srs->{SUBJECT}\"\t$srs->{START_S}\t$srs->{END_S}
#\tE_value $srs->{E_VALUE}\tStrand $srs->{STRAND_S}
#\tFrame $srs->{FRAME_S}\t\#Projection $srs->{PROJECTION} 
my $GFFstring = ("\%s\t" x 8).'Target "%s"'.
       "\t\%s\t\%s\tE_Value .\tStrand \%s\tFrame \%s\n";
#line 300 "/home/ug/jabril/development/softjabril/blat2gff/blat2gff.nw"
#
# MAIN LOOP
#
#line 337 "/home/ug/jabril/development/softjabril/blat2gff/blat2gff.nw"
while (<STDIN>) {
    my @f = ();
    next unless /^\d/o;
    chomp;
    @f = split /\s+/og, $_;
    do { # ensure that there are 21 fields
        print STDERR "### NOT enough fields for this record:\n##### @f \n";
        next;
    } unless scalar(@f) == $maxfields;
    %gff = ();
    ++$rec_cnt;
    $f[$blat{STRAND}] =~ /(.)(.)/og && 
        ($gff{Q_STRAND} = $1, $gff{T_STRAND} = $2);
    do { # ensure that strands are +/- 
        print STDERR "### CANNOT find strand for this record:\n##### @f \n";
        next;
    } unless ($gff{Q_STRAND} =~ /[+-]/o && $gff{T_STRAND} =~ /[+-]/o);   
    $f[$blat{Q_NAME}] =~ m{/?([^/]+)$}og && ($gff{Q_NAME} = $1) ;
    $f[$blat{T_NAME}] =~ m{/?([^/]+)$}og && ($gff{T_NAME} = $1) ;
    defined($group{"$gff{T_NAME}.$gff{T_STRAND}"}) || 
        ($group{"$gff{T_NAME}.$gff{T_STRAND}"} = ++$group_cnt);
    $gff{GROUP} = $group{"$gff{T_NAME}.$gff{T_STRAND}"};
    $gff{Q_LEN} = $f[$blat{Q_LEN}]; # seq lengths are OK
    $gff{T_LEN} = $f[$blat{T_LEN}];
    $gff{Q_ORI} = $f[$blat{Q_ORI}] + 1; # HSP coords start at 0, not at 1
    $gff{Q_END} = $f[$blat{Q_END}] + 1;
    $gff{T_ORI} = $f[$blat{T_ORI}] + 1;
    $gff{T_END} = $f[$blat{T_END}] + 1;
    $gff{Q_FRAME} =
        &get_frame($gff{Q_STRAND}, $gff{Q_LEN}, $gff{Q_ORI}, $gff{Q_END});
    $gff{T_FRAME} = 
        &get_frame($gff{T_STRAND}, $gff{T_LEN}, $gff{T_ORI}, $gff{T_END});
    $gff{B_NUM} = $f[$blat{B_NUM}];
    @{ $gff{B_LEN} }  = split /,/og, $f[$blat{B_LEN}];
    @{ $gff{B_QPOS} } = split /,/og, $f[$blat{B_QPOS}];
    @{ $gff{B_TPOS} } = split /,/og, $f[$blat{B_TPOS}];
    printf STDOUT "# ".$GFFstring,
        $gff{Q_NAME}, "BLAT","$gff{Q_LEN}:$gff{T_LEN}", 
            $gff{Q_ORI}, $gff{Q_END}, 0, $gff{Q_STRAND}, $gff{Q_FRAME},
        "$gff{T_NAME}.$gff{GROUP}.$rec_cnt",
            $gff{T_ORI}, $gff{T_END}, $gff{T_STRAND}, $gff{T_FRAME};
    &loop_HSPs($rec_cnt);
}; # while read input

exit(0);
#line 304 "/home/ug/jabril/development/softjabril/blat2gff/blat2gff.nw"
#
# FUNCTIONS
#
#line 385 "/home/ug/jabril/development/softjabril/blat2gff/blat2gff.nw"
sub get_frame() {
    my ($strand,$len,$ori,$end) = @_;
    if ($strand eq '-') {
       return $frame[(($len - $end + 1) % 3)];
    } else {
       return $frame[($ori % 3)];
    }
} # get_frame
sub get_hsp() {
    my ($s,$L,$l,$n) = @_;
    my ($o,$e);
    if ($s eq '-') {
        $e = $L - $n;
        $o = $e - $l
    } else {
        $o = $n;
        $e = $o + $l;
    };
    $o++; $e++; # HSP coords start at 0, not at 1
	return ($o, $e, &get_frame($s,$L,$o,$e));
} # get_hsp
#line 419 "/home/ug/jabril/development/softjabril/blat2gff/blat2gff.nw"
sub loop_HSPs() {
    my ($num) = @_;
    my ($c,$l,$q,$t,$qo,$qe,$qf,$to,$te,$tf);
    for ($c = 0; $c < $gff{B_NUM}; $c++) {
        $l = $gff{B_LEN}[$c];
        $q = $gff{B_QPOS}[$c];
        $t = $gff{B_TPOS}[$c];
        ($qo,$qe,$qf) = &get_hsp($gff{Q_STRAND},$gff{Q_LEN},$l,$q);
        ($to,$te,$tf) = &get_hsp($gff{T_STRAND},$gff{T_LEN},$l,$t);
        printf STDOUT $GFFstring,
            $gff{Q_NAME}, "blat", "hsp", $qo, $qe, 0, $gff{Q_STRAND}, $qf,
            "$gff{T_NAME}.$gff{GROUP}.$num", $to, $te, $gff{T_STRAND}, $tf;
    }; # for
} # loop_HSPs
