#!/usr/local/bin/perl -w
# This is perl, v5.6.1 built for i686-linux
# /usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %                          GFF2APLOT                               %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
#    Converting alignments in GFF format to PostScript dotplots.
# 
#     Copyright (C) 1999-2003 - Josep Francesc ABRIL FERRANDO  
#                                       Thomas WIEHE                   
#                                      Roderic GUIGO SERRA       
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
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# $Id: gff2aplot.pl,v 1.3 2003-03-04 18:05:38 jabril Exp $
#
#
# MODULES
#
use strict;
use vars qw/
  $T $F $n $c %Messages
  %DefaultVars %CmdLineVars %CustomVars %Defaults %Vars
  %GFF_DATA %ALN_DATA %Order
  $PROGRAM $VERSION $REVISION $REVISION_DATE $LAST_UPDATE
  %VarKeys $noCV
  %TESTS %Ribbons %Shapes %GroupShapes %LineStyles %TextAlign
  %Fonts %FORMATS $formats %COLORS $colors %GRADIENTS %UNITS %USED
  $kilo $mega $giga $cm2pt $mm2pt $in2pt
  %SpecialVars %SpecialCache
  $GFF $GFF_NOGP $VECTOR $ALIGN $APLOT $APLOT_NOGP $noGFF
  $aln_COUNT $seq_COUNT
  $_counter $_prop $_element
  $_order $_elemNum $_ori $_end $_flw $_glw
  $_mnsco $_mxsco $_nori $_nend
  $_fttype $_ftprop $_ftname $_ftid
  $_ftori  $_ftend  $_ftsco  $_ftfrm
  $_ftnori $_ftnend $_ftnsco $_ftnfrm
  %StringCache
  $regexp_natural $regexp_integer $regexp_float $regexp_real
  $gradregex
  $regexp_frame
  $regexp_strand
  $regexp_group
  /;

#
$SIG{HUP}  = \&trap_signals_prog;
$SIG{ABRT} = \&trap_signals;
$SIG{INT}  = \&trap_signals;
$SIG{QUIT} = \&trap_signals;
$SIG{TERM} = \&trap_signals;
$SIG{KILL} = \&trap_signals;
$SIG{CHLD} = 'IGNORE';

#
use Getopt::Long;
Getopt::Long::Configure qw/ bundling /;
use Data::Dumper;
local $Data::Dumper::Purity   = 0;
local $Data::Dumper::Deepcopy = 1;
use Benchmark;
my @Timer = ( new Benchmark );

#
# CONSTANTS
#
( $T, $F ) = ( 1, 0 );    # for 'T'rue and 'F'alse
($PROGRAM,$VERSION,$REVISION,$REVISION_DATE,$LAST_UPDATE) =
  ( 'gff2aplot', 'v2.0',
  '$Revision: 1.3 $',                                            #'
  '$Date: 2003-03-04 18:05:38 $',                                 #'
  '$Id: gff2aplot.pl,v 1.3 2003-03-04 18:05:38 jabril Exp $',    #'
);
$REVISION      =~ s/\$//og;
$REVISION_DATE =~ s/\$//og;
$noCV    = '?';
%VarKeys = (
    L => 'LAYOUT',
    Q => 'SEQUENCE',
    S => 'SOURCE',
    T => 'STRAND',
    G => 'GROUP',
    F => 'FEATURE',
    X => 'EXTRA',
);
( $cm2pt, $mm2pt, $in2pt ) = ( 28.35, 2.835, 72 );
( $kilo, $mega, $giga ) = ( 10**3, 10**6, 10**9 );
$regexp_natural = '^\d+$';                                               #'
$regexp_integer = '^[+-]?\d+$';                                          #'
$regexp_float   = '^[+-]?(?:\d+(?:\.\d*)?|\.\d+)$';                      #'
$regexp_real    = '^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$';    #'
$gradregex      = '(gradient|rainbow)';
my $torad = 3.1415926535898 / 180;
( $GFF, $GFF_NOGP, $VECTOR, $ALIGN, $APLOT, $APLOT_NOGP, $noGFF ) =
  qw/ X x V A O o ? /;
$regexp_frame  = '^[.0123]$';     #' Accepting blast-formated frames (.123)
$regexp_strand = '^[\.\+\-]$';    #'
$regexp_group = '^(.*?)(?:"(.+?)"(?:\s+(.+))?)?$';    #'
my $Error = "\<\<\<  ERROR  \>\>\> ";
my $Warn  = "\<\<\< WARNING \>\>\> ";
my $spl   = "\<\<\<\-\-\-\-\-\-\-\-\-\>\>\>\n";
my $spw   = "\<\<\<         \>\>\> ";
my $line  = ( "#" x 80 ) . "\n";
my $sp    = "###\n";

#
# VARIABLES
#
my $Custom_path =
  defined( $ENV{GFF2APLOT_CUSTOMPATH} ) ? $ENV{GFF2APLOT_CUSTOMPATH} : '.';
my $Custom_file =
  defined( $ENV{GFF2APLOT_CUSTOMFILE} ) 
  ? $ENV{GFF2APLOT_CUSTOMFILE}
  : '.gff2aplotrc';
$Custom_path =~ s{/$}{}o;
my ( $Debug, $Verbose, $Quiet, $LogFile, $logs_filename ) =
  ( $F, $F, $F, $F, undef );
my ( @data_files, $file );
my @custom_files   = ();
my $ExtraCustomFlg = $T;
%TESTS = (
    BOOLEAN => sub {    # Value, varName, Lowercase
        my ( $_v, $_n, $_l ) = @_;
        $_l =~ /^(1|on|t(rue)|y(es))$/o  && ( $$_v = $T, return $T );
        $_l =~ /^(0|off|f(alse)|n(o))$/o && ( $$_v = $F, return $T );
        &warn( 'NOT_A_BOOLEAN', $F, $_n );
        return $F;
    },
    ALPHA => sub {    # Value, varName
        my ( $_v, $_n, undef ) = @_;
        $$_v = lc($$_v);
        $$_v =~ /^[a-z][a-z_0-9]*$/o && return $T;
        &warn( 'NOT_ALPHA', $F, $_n );
        return $F;
    },
    STRING => sub {    # Value
        my ( $_v, undef, undef ) = @_;
        defined($$_v) && do {
            $$_v =~ s{[\\]}{\\134}og;
            $$_v =~ s{[\(]}{\\050}og;
            $$_v =~ s{[\)]}{\\051}og;
        };
        return $T;
    },
    NATURAL => sub {    # Value, varName
        my ( $_v, $_n, undef ) = @_;
        $$_v =~ /$regexp_natural/o && return $T;
        &warn( 'NOT_A_NATURAL', $F, $_n );
        return $F;
    },
    INTEGER => sub {    # Value, varName
        my ( $_v, $_n, undef ) = @_;
        $$_v =~ /$regexp_integer/o && return $T;
        &warn( 'NOT_AN_INTEGER', $F, $_n );
        return $F;
    },
    FLOAT => sub {    # Value, varName
        my ( $_v, $_n, undef ) = @_;
        $$_v =~ /$regexp_float/o && return $T;
        &warn( 'NOT_A_FLOAT', $F, $_n );
        return $F;
    },
    REAL => sub {    # Value, varName
        my ( $_v, $_n, undef ) = @_;
        $$_v =~ /$regexp_real/o && return $T;
        &warn( 'NOT_A_REAL', $F, $_n );
        return $F;
    },
    COLOR => sub {    # Value, varName, Lowercase
        my ( $_v, $_n, $_l ) = @_;
        $_l =~ /^b(ack)?g(round)?(color)?/o && ( $$_v = 'bg', return $T );
        $_l =~ /^f(ore)?g(round)?(color)?/o && ( $$_v = 'fg', return $T );
        defined( $COLORS{$_l} ) && do {
            $$_v = $_l;
            defined( $USED{'COLORS'}{$_l} )
              || ( $USED{'COLORS'}{$_l} = $COLORS{$_l}->[0] );
            return $T;
        };    # defined
        &warn( 'NOT_A_COLOR', $F, $_l, $_n );
        return $F;
    },
    EXTCOLOR => sub {    # Value, varName, Lowercase
        my $maxcolornum = 6;
        my ( $_v, $_n, $_l ) = @_;
        my ( $_vt, $_lt, @_cl, @_ll, $_t );
        $_l =~ /^(.*)$gradregex$/ && do {
            $_vt = defined($1) ? $1 : '';
            $_lt = defined($2) ? $2 : '';
            defined( $GRADIENTS{$_lt} ) && do {
                defined( $USED{'GRADIENTS'}{$_lt} )
                  || ( $USED{'GRADIENTS'}{$_lt} = $_vt );
                if ( $_lt eq 'rainbow' ) {
                    $$_v = $_lt;
                }
                else {
                    $_vt =~ s/\s*$//og;
                    @_cl = split /\s+|_/og, $_vt;
                    @_ll = ();
                    for ( my $l = 0 ; ( $l <= $#_cl && $l < $maxcolornum ) ;
                        $l++ )
                    {
                        $_t = '';
                        $TESTS{COLOR}->( \$_t, $_n, $_cl[$l] )
                          && ( push @_ll, $_t );
                    }
                    scalar(@_ll) > 0 || ( push @_ll, 'fg' );
                    $$_v = join ( " ", @_ll, scalar(@_ll), $_lt );
                }
                return $T;
            };    # defined
            &warn( 'NOT_A_GRADIENT', $F, $_l, $_n );
            return $F;
        };
        return $TESTS{COLOR}->( $_v, $_n, $_l );
    },
    CMYKFLOAT => sub {    # Value, varName, Lowercase
        my ( $_v, $_n, undef ) = @_;
        ( $$_v =~ /$regexp_float/o && ( $$_v >= 0 && $$_v <= 1 ) ) && return $T;
        &warn( 'NOTCMYKFLOAT', $F, $_n, $$_v );
        return $F;
    },
    FONT => sub {    # Value, varName, Lowercase
        my ( $_v, $_n, $_l ) = @_;
        defined( $Fonts{$_l} ) && ( $$_v = $_l, return $T );
        &warn( 'NOT_A_FONT', $F, $$_v, $_n );
        return $F;
    },
    BBOX => sub {    # Value, varName, Lowercase
        my ( $_v, $_n, $_l ) = @_;
        my @it = ();
        @it = split /,/og, $_l;
        scalar(@it) != 2 && do {
            &warn( 'BBOX_NOTOK', $F, $_n, $_l, $_n );
            return $F;
        };    # not enough params
        $TESTS{PS_UNIT}->( \$it[0], $_n, $it[0] ) || return $F;
        $TESTS{PS_UNIT}->( \$it[1], $_n, $it[1] ) || return $F;
        $$_v = [ 'userdef', @it ];
        return $T;
    },
    PAGE => sub {    # Value, varName, Lowercase
        my ( $_v, $_n, $_l ) = @_;
        defined( $FORMATS{$_l} ) && do {
            $$_v = $_l;
            defined( $USED{'FORMATS'}{$_l} )
              || ( $USED{'FORMATS'}{$_l} = $FORMATS{$_l}->[0] );
            return $T;
        };
        &warn( 'NOT_A_PAGE', $F, $_l, $$_v );
        return $F;
    },
    SEQ_UNIT => sub {    # Value, varName, Lowercase
        my ( $_v, $_n, $_l ) = @_;
        my $q = undef;
        $_l =~ /^\.$|^\s*$/o && do {
            &warn( 'NOT_A_NUMBER', $F, $_n );
            return $F;
        };
        $_l =~ s/(?:\.)?
                   (?:(g(?:iga)?|m(?:ega)?|k(?:ilo)?)?
                   b(?:ase)?
                   (?:p(?:air)?)?
                   (?:s)?)$
                  //ox && ( $q = $1 );
        defined($q) || do {
            $_l =~ /[^\d]+$/o && do {
                &warn( 'NOT_A_DNAUNIT', $F, $_n, $$_v );
                return $F;
            };
            $q = '';
        };
        $TESTS{FLOAT}->( \$_l, $_n, $_l ) || return $F;
        ( $q, undef ) = split //o, $q, 2;
        $q .= 'bp';
        $$_v = [ $_l, $q ];
        return $T;
    },
    PS_UNIT => sub {    # Value, varName, Lowercase
        my ( $_v, $_n, $_l ) = @_;
        my $q = undef;
        $_l =~ /^\.$|^\s*$/o && do {
            &warn( 'NOT_A_NUMBER', $F, $_n );
            return $F;
        };
        PSLEN: {
            $_l =~ s/m(eter(s)?)?$//o && do {
                $_l =~ s/(\.)?c(enti)?$//o && ( $q = 'cm', last PSLEN );
                $_l =~ s/(\.)?m(ili)?$//o  && ( $q = 'mm', last PSLEN );
                last PSLEN;
            };
            $_l =~ s/(\.)?in(ch(es)?)?$//o && ( $q = 'in', last PSLEN );
            $_l =~ s/(\.)?(pt|point(s)?)$//o && ( $q = 'pt' );
        };    # PSLEN
        defined($q) || do {
            $_l =~ /[^\d]+$/o && do {
                &warn( 'NOT_A_PSUNIT', $F, $_n, $$_v );
                return $F;
            };
            $q = 'pt';
        };
        $TESTS{FLOAT}->( \$_l, $_n, $_l ) || return $F;
        $$_v = [ $_l, $q ];
        return $T;
    },
    RANGE => sub {    # Value, varName, Lowercase
        my ( $_v, $_n, $_l ) = @_;
        my @it = ();
        @it = split /\.\./og, $_l;
        scalar(@it) != 2 && do {
            &warn( 'RANGE_NOTOK', $F, $_n, $$_v, $_n );
            return $F;
        };    # not enough params
        $it[0] eq '*' || do {
            $TESTS{FLOAT}->( \$it[0], $_n, $it[0] ) || return $F;
        };    # $it[0] -> lower-value
        $it[1] eq '*' || do {
            $TESTS{FLOAT}->( \$it[1], $_n, $it[1] ) || return $F;
        };    # $it[1] -> upper-value
        $$_v =
          [ $it[0] eq '*' ? undef: $it[0], $it[1] eq '*' ? undef: $it[1], ];
        return $T;
    },
    PSRANGE => sub {    # Value, varName, Lowercase
        my ( $_v, $_n, $_l ) = @_;
        my @it = ();
        @it = split /\.\./og, $_l;
        scalar(@it) != 2 && do {
            &warn( 'RANGE_NOTOK', $F, $_n, $$_v, $_n );
            return $F;
        };    # not enough params
        $it[0] eq '*' || do {
            $TESTS{PS_UNIT}->( \$it[0], $_n, $it[0] ) || return $F;
        };    # $it[0] -> lower-value
        $it[1] eq '*' || do {
            $TESTS{PS_UNIT}->( \$it[1], $_n, $it[1] ) || return $F;
        };    # $it[1] -> upper-value
        $$_v =
          [ $it[0] eq '*' ? undef: [ @{ $it[0] } ],
            $it[1] eq '*' ? undef: [ @{ $it[1] } ] ];
        return $T;
    },
    SEQRANGE => sub {    # Value, varName, Lowercase
        my ( $_v, $_n, $_l ) = @_;
        my @it = ();
        @it = split /\.\./og, $_l;
        scalar(@it) != 2 && do {
            &warn( 'RANGE_NOTOK', $F, $_n, $$_v, $_n );
            return $F;
        };    # not enough params
        $it[0] eq '*' || do {
            $TESTS{SEQ_UNIT}->( \$it[0], $_n, $it[0] ) || return $F;
        };    # $it[0] -> lower-value
        $it[1] eq '*' || do {
            $TESTS{SEQ_UNIT}->( \$it[1], $_n, $it[1] ) || return $F;
        };    # $it[1] -> upper-value
        $$_v =
          [ $it[0] eq '*' ? undef: [ @{ $it[0] } ],
            $it[1] eq '*' ? undef: [ @{ $it[1] } ] ];
        return $T;
    },
    LINESTY => sub {    # Value, varName, Lowercase
        my ( $_v, $_n, $_l ) = @_;
        ( defined($$_v) && defined( $LineStyles{$_l} ) )
          && ( $$_v = $_l, return $T );
        &warn( 'LINESTY_NOTOK', $F, $$_v, $_n );
        return $F;
    },
    RIBBON => sub {    # Value, varName, Lowercase
        my ( $_v, $_n, $_l ) = @_;
        ( defined($$_v) && defined( $Ribbons{$_l} ) )
          && ( $$_v = $_l, return $T );
        &warn( 'RIBBON_NOTOK', $F, $$_v, $_n );
        return $F;
    },
    SHAPE => sub {    # Value, varName, Lowercase
        my ( $_v, $_n, $_l ) = @_;
        ( defined($$_v) && defined( $Shapes{$_l} ) )
          && ( $$_v = $_l, return $T );
        &warn( 'SHAPE_NOTOK', $F, $$_v );
        return $F;
    },
    GPSHAPE => sub {    # Value, varName, Lowercase
        my ( $_v, $_n, $_l ) = @_;
        ( defined($$_v) && defined( $GroupShapes{$_l} ) )
          && ( $$_v = $_l, return $T );
        &warn( 'GRP_SHAPE_NOTOK', $F, $$_v );
        return $F;
    },
    TXTALIGN => sub {    # Value, varName, Lowercase, flag
        my ( $_v, $_n, $_l, $a_flg ) = @_;
        ( defined($$_v)
          && defined( $TextAlign{ ( $a_flg ? 'H' : 'V' ) }{$_l} ) )
          && ( $$_v = $_l, return $T );
        &warn( 'ALIGN_NOTOK', $F, $_n, $$_v );
        return $F;
    },
    STRAND => sub {    # Value, varName, Lowercase
        my ( $_v, $_n, undef ) = @_;
        $$_v =~ /$regexp_strand/o && return $T;
        &warn( 'NOT_A_STRAND', $F, $$_v, $_n );
        return $F;
    },
);
%UNITS = (
    "pt"  => sub ($) { return $_[0] },
    "mm"  => sub ($) { return $_[0] * $mm2pt },
    "cm"  => sub ($) { return $_[0] * $cm2pt },
    "in"  => sub ($) { return $_[0] * $in2pt },
    "bp"  => sub ($) { return $_[0] },
    "kbp" => sub ($) { return $_[0] * $kilo },
    "mbp" => sub ($) { return $_[0] * $mega },
    "gbp" => sub ($) { return $_[0] * $giga },
);
%Fonts = (

    # serif
    'default'          => '/Times-Roman',
    'times'            => '/Times-Roman',
    'times-roman'      => '/Times-Roman',
    'times-italic'     => '/Times-Italic',
    'times-bold'       => '/Times-Bold',
    'times-bolditalic' => '/Times-BoldItalic',

    # sans serif
    'helvetica'             => '/Helvetica',
    'helvetica-oblique'     => '/Helvetica-Oblique',
    'helvetica-bold'        => '/Helvetica-Bold',
    'helvetica-boldoblique' => '/Helvetica-BoldOblique',

    # monospaced
    'courier'             => '/Courier',
    'courier-oblique'     => '/Courier-Oblique',
    'courier-bold'        => '/Courier-Bold',
    'courier-boldoblique' => '/Courier-BoldOblique',
);
%LineStyles = (
    's'             => 'solid',
    'solid'         => 'solid',
    'd'             => 'dotted',
    'dotted'        => 'dotted',
    'ld'            => 'ldotted',
    'long_dotted'   => 'ldotted',
    'l'             => 'dashed',
    'dashed'        => 'dashed',
    'll'            => 'ldashed',
    'long_dashed'   => 'ldashed',
    'dd'            => 'ddashed',
    'dashed_dotted' => 'ddashed',
);    # %LineStyles
%Ribbons = (    # shw lns
    'n'       => [ $F, $F ],
    'none'    => [ $F, $F ],
    'l'       => [ $F, $T ],
    'lines'   => [ $F, $T ],
    'r'       => [ $T, $F ],
    'ribbons' => [ $T, $F ],
    'b'       => [ $T, $T ],
    'both'    => [ $T, $T ],
);    # %Ribbons
%Shapes = (
    'none'                 => 'pt',
    'box'                  => 'fbox',
    'half_box'             => 'hbox',
    'diamond'              => 'fdmd',
    'sand_clock'           => 'fmdm',
    'arrow'                => 'fraw',
    'right_arrow'          => 'fraw',
    'half_arrow'           => 'hraw',
    'half_right_arrow'     => 'hraw',
    'left_arrow'           => 'flaw',
    'half_left_arrow'      => 'hlaw',
    'arrow_end'            => 'frae',
    'right_arrow_end'      => 'frae',
    'half_arrow_end'       => 'hrae',
    'half_right_arrow_end' => 'hrae',
    'left_arrow_end'       => 'flae',
    'half_left_arrow_end'  => 'hlae',
    'single'               => 'frsg',
    'right_single'         => 'frsg',
    'half_single'          => 'hrsg',
    'half_right_single'    => 'hrsg',
    'left_single'          => 'flsg',
    'half_left_single'     => 'hlsg',
    'triangle'             => 'frtn',
    'right_triangle'       => 'frtn',
    'half_triangle'        => 'hrtn',
    'half_right_triangle'  => 'hrtn',
    'left_triangle'        => 'fltn',
    'half_left_triangle'   => 'hltn',
    'triangle_up'          => 'futn',
    'triangle_down'        => 'fdtn',
    'circle'               => 'fcir',
    'arc'                  => 'farc',
    'half_arc'             => 'harc',
);    # %Shapes
%GroupShapes = (
    'none'                      => 'pt',
    'line'                      => 'glns',
    'dotted_line'               => 'glnd',
    'bracket'                   => 'gbrk',
    'curly_bracket'             => 'gbrc',
    'segment'                   => 'gsgm',
    'segment_dotted'            => 'gsgd',
    'arrow'                     => 'gfaw',
    'stroked_arrow'             => 'gsaw',
    'filled_arrow'              => 'gfaw',
    'stroked_arrow_dotted'      => 'gsad',
    'filled_arrow_dotted'       => 'gfad',
    'half_arrow'                => 'ghaw',
    'half_filled_arrow'         => 'ghaw',
    'half_filled_arrow_dotted'  => 'ghad',
    'half_stroked_arrow'        => 'ghsw',
    'half_stroked_arrow_dotted' => 'ghsd',
);    # %GroupShapes
%TextAlign = (
    H => {
        'l'      => 'lh',
        'left'   => 'lh',
        'c'      => 'ch',
        'center' => 'ch',
        'r'      => 'rh',
        'right'  => 'rh',
    },
    V => {
        't'      => 'tv',
        'top'    => 'tv',
        'c'      => 'cv',
        'center' => 'cv',
        'b'      => 'bv',
        'bottom' => 'bv',
    },
);    # %TextAlign
%SpecialVars = (
    font => sub {
        my ( $cstr, $cflg ) = @_;
        my @fld;
        $cflg || do {
            &warn( 'CANNOTDEFNOW', $F, 'CUSTOM FONT ALIAS', $cstr );
            return $F;
        };    # $cflg
        @fld = split /\s+/og, $cstr;
        ( defined( $fld[0] ) && defined( $fld[1] ) ) || do {
            &warn( 'WRONGSPECREC', $F, 'FONT', 2,
                'font  <font-alias> <PostScript-font-name>' );
            return $F;
        };    # defined
        $fld[1] =~ s{^/}{}o;
        $Fonts{ 'my_' . $fld[0] } = '/' . $fld[1];
        return $T;
    },    # font
    color => sub {
        my ( $cstr, $cflg ) = @_;
        my @fld;
        my $_n = 'Special Features [CMYK COLOR] ';
        $cflg || do {
            &warn( 'CANNOTDEFNOW', $F, 'CUSTOM CMYK COLOR', $cstr );
            return $F;
        };    # $cflg
        @fld = split /\s+/og, $cstr;
        scalar(@fld) < 5 && do {
            &warn( 'WRONGSPECREC', $F, 'CMYK COLOR', 5,
                'color  <color-alias> <cyan> <magenta> <yellow> <black>' );
            return $F;
        };    # defined
              # CMYKFLOAT checking color amount between 0..1
        $TESTS{CMYKFLOAT}->( \$fld[1], $_n . '(cyan)', $fld[1] ) || return $F;
        $TESTS{CMYKFLOAT}->( \$fld[2], $_n . '(magenta)', $fld[2] )
          || return $F;
        $TESTS{CMYKFLOAT}->( \$fld[3], $_n . '(yellow)', $fld[3] ) || return $F;
        $TESTS{CMYKFLOAT}->( \$fld[4], $_n . '(black)',  $fld[4] ) || return $F;
        $COLORS{ ( 'my_' . $fld[0] ) } = [ ++$colors, @fld[ 1 .. 5 ] ];
        return $T;
    },    # color
    box => sub {
        my ( $cstr, undef ) = @_;
        my ( @fld, @tld );
        my $_n = 'Special Features [BOX] ';
        @fld = split /\s+/og, $cstr, 5;
        scalar(@fld) < 4 && do {
            &warn( 'WRONGSPECREC', $F, 'BOX', 4,
                'box  <Xori> <Xend> <Yori> <Yend> [ L:<line_color> ] '
                . '[ S:<line_style> ] [ W:<line_width> ] [ F:<fill_color> ] [ K:<score> ]'
            );
            return $F;
        };    # defined
        $TESTS{SEQ_UNIT}->( \$fld[0], $_n . '(Xori)', lc( $fld[0] ) )
          || return $F;
        $TESTS{SEQ_UNIT}->( \$fld[1], $_n . '(Xend)', lc( $fld[1] ) )
          || return $F;
        $TESTS{SEQ_UNIT}->( \$fld[2], $_n . '(Yori)', lc( $fld[2] ) )
          || return $F;
        $TESTS{SEQ_UNIT}->( \$fld[3], $_n . '(Yend)', lc( $fld[3] ) )
          || return $F;
        @tld = &get_spfeat_props( $fld[4], $F, $T, $_n );
        ( $tld[2] =~ /$gradregex/o ) && ( $tld[2] = "$tld[11] $tld[2]" );
        ( $tld[5] =~ /$gradregex/o ) && ( $tld[5] = "$tld[11] $tld[5]" );
        push @{ $SpecialCache{PRE} },
          [
            &get_units( \@{ $fld[0] } ),
            &get_units( \@{ $fld[1] } ),
            &get_units( \@{ $fld[2] } ),
            &get_units( \@{ $fld[3] } ),
            join (
                " ", @tld[ 2, 4, 3, 5 ],
                &get_units( \@{ $fld[0] } ), &get_units( \@{ $fld[1] } ),
                &get_units( \@{ $fld[2] } ), &get_units( \@{ $fld[3] } )
            ) . " xbox"
        ];
        return $T;
    },    # box
    circle => sub {
        my ( $cstr, undef ) = @_;
        my ( @fld, @tld );
        my $_n = 'Special Features [CIRCLE] ';
        @fld = split /\s+/og, $cstr, 5;
        scalar(@fld) < 4 && do {
            &warn( 'WRONGSPECREC', $F, 'CIRCLE', 4,
                'circle  <Cx> <Cy> <Rx> <Ry> '
                . '[ A:<angle> ] [ L:<line_color> ] [ S:<line_style> ] '
                . '[ W:<line_width> ] [ F:<fill_color> ] [ K:<score> ]' );
            return $F;
        };    # defined
        $TESTS{SEQ_UNIT}->( \$fld[0], $_n . '(Cx)', lc( $fld[0] ) )
          || return $F;
        $TESTS{SEQ_UNIT}->( \$fld[1], $_n . '(Cy)', lc( $fld[1] ) )
          || return $F;
        $TESTS{SEQ_UNIT}->( \$fld[2], $_n . '(Rx)', lc( $fld[2] ) )
          || return $F;
        $TESTS{SEQ_UNIT}->( \$fld[3], $_n . '(Ry)', lc( $fld[3] ) )
          || return $F;
        @tld = &get_spfeat_props( $fld[4], $F, $T, $_n );
        ( $tld[2] =~ /$gradregex/o ) && ( $tld[2] = "$tld[11] $tld[2]" );
        ( $tld[5] =~ /$gradregex/o ) && ( $tld[5] = "$tld[11] $tld[5]" );
        push @{ $SpecialCache{PRE} },
          [
            &get_units( \@{ $fld[0] } ),
            &get_units( \@{ $fld[0] } ),
            &get_units( \@{ $fld[1] } ),
            &get_units( \@{ $fld[1] } ),
            join (
                " ", @tld[ 2, 4, 3, 5 ],
                &get_units( \@{ $fld[2] } ), &get_units( \@{ $fld[3] } ),
                $tld[1],                     &get_units( \@{ $fld[0] } ),
                &get_units( \@{ $fld[1] } )
            ) . " xcir"
        ];
        return $T;
    },    # circle
    arrow => sub {
        my ( $cstr, undef ) = @_;
        my ( @fld, @tld );
        my $_n = 'Special Features [ARROW] ';
        @fld = split /\s+/og, $cstr, 4;
        scalar(@fld) < 3 && do {
            &warn( 'WRONGSPECREC', $F, 'ARROW', 3,
                'arrow  <Cx> <Cy> <length> [ A:<angle> ] '
                . '[ L:<line_color> ] [ S:<line_style> ] [ W:<line_width> ] '
                . '[ N:<font_name> ] [ Z:<font_size> ] [ H:<horizontal_align> ] '
                . '[ V:<vertical_align> ] [ T:<text_color> ] [ "<string>" ]' );
            return $F;
        };    # defined
        $TESTS{SEQ_UNIT}->( \$fld[0], $_n . '(Cx)', lc( $fld[0] ) )
          || return $F;
        $TESTS{SEQ_UNIT}->( \$fld[1], $_n . '(Cy)', lc( $fld[1] ) )
          || return $F;
        $TESTS{PS_UNIT}->( \$fld[2], $_n . '(length)', lc( $fld[2] ) )
          || return $F;
        @tld = &get_spfeat_props( $fld[3], $T, $T, $_n );
        ( $tld[2] =~ /$gradregex/o ) && ( $tld[2] = "$tld[11] $tld[2]" );
        $tld[0] ne '' || ( $tld[0] = '( )' );
        push @{ $SpecialCache{POST} },
          [
            &get_units( \@{ $fld[0] } ),
            &get_units( \@{ $fld[0] } ),
            &get_units( \@{ $fld[1] } ),
            &get_units( \@{ $fld[1] } ),
            '[ ' . join (
                " ",                         @tld[ 0, 6 .. 10, 2, 4 ],
                &get_units( \@{ $fld[2] } ), $tld[1],
                &get_units( \@{ $fld[0] } ), &get_units( \@{ $fld[1] } )
            ) . " xarw"
        ];
        return $T;
    },    # arrow
    text => sub {
        my ( $cstr, undef ) = @_;
        my ( @fld, @tld );
        my $_n = 'Special Features [TEXT] ';
        @fld = split /\s+/og, $cstr, 3;
        scalar(@fld) < 2 && do {
            &warn( 'WRONGSPECREC', $F, 'TEXT', 3,
                'text  <Cx> <Cy> "<string>" [ A:<angle> ] '
                . '[ N:<font_name> ] [ Z:<font_size> ] [ H:<horizontal_align> ] '
                . '[ V:<vertical_align> ] [ T:<text_color> ]' );
            return $F;
        };    # defined
        $TESTS{SEQ_UNIT}->( \$fld[0], $_n . '(Cx)', lc( $fld[0] ) )
          || return $F;
        $TESTS{SEQ_UNIT}->( \$fld[1], $_n . '(Cy)', lc( $fld[1] ) )
          || return $F;
        @tld = &get_spfeat_props( $fld[2], $T, $F, $_n );
        $tld[0] ne '' && do {
            push @{ $SpecialCache{POST} },
              [ &get_units( \@{ $fld[0] } ), &get_units( \@{ $fld[0] } ),
                &get_units( \@{ $fld[1] } ), &get_units( \@{ $fld[1] } ),
                '[ '
                . join ( " ", @tld[ 1, 0, 6 .. 10 ],
                    &get_units( \@{ $fld[0] } ), &get_units( \@{ $fld[1] } ) )
                . " xtxt" ];
            return $T;
        };
        &warn( 'EMPTYSTRING', $F, $_n );
        return $F;
    },    #
);
my ( $seqname, $source, $feature, $start, $end, $score, $strand, $frame )
  ;    # GFF temporary vars
my (
    $seqname_1, $seqname_2, $source_1, $source_2, $feature_1, $feature_2,
    $start_1,   $start_2,   $end_1,    $end_2,    $score_1,   $score_2,
    $strand_1,  $strand_2,  $frame_1,  $frame_2
);    # APLOT temporary vars
my ( $tag, $group, $group_id, $label, $group_gff_counter,
    $group_aplot_counter );    # GROUPING temporary vars
( $_counter, $_prop, $_element ) = ( 0 .. 2 );
( $_order, $_elemNum, $_ori, $_end,  $_mnsco,
  $_mxsco, $_flw,     $_glw, $_nori, $_nend )
  = ( 0 .. 9 );
( $_fttype, $_ftprop, $_ftname, $_ftid,   $_ftori,  $_ftend,
   $_ftsco, $_ftfrm,  $_ftnori, $_ftnend, $_ftnsco, $_ftnfrm )
  = ( 0 .. 11 );
my @vect_ary;
my %vect_type = (
    'sco'   => 1,
    'pos'   => 2,
    'box'   => 3,
    'chr'   => 0,
    'ERROR' => -1
);
my %StrandNum;
my $ttt = 0;

foreach my $w (qw/ + ++ +- +. .+ . .. .- -. -+ -- - /) {
    $StrandNum{$w} = $ttt++;
};    # foreach
$colors = 0;
%COLORS = (    # [ ColorNUMBER, qw/ CYAN MAGENTA YELLOW BLACK / ]
        # black+grey+white
    black         => [ ++$colors, qw/ 0.00 0.00 0.00 1.00 / ],
    verydeepgrey  => [ ++$colors, qw/ 0.00 0.00 0.00 0.90 / ],
    deepgrey      => [ ++$colors, qw/ 0.00 0.00 0.00 0.80 / ],
    verydarkgrey  => [ ++$colors, qw/ 0.00 0.00 0.00 0.70 / ],
    darkgrey      => [ ++$colors, qw/ 0.00 0.00 0.00 0.60 / ],
    grey          => [ ++$colors, qw/ 0.00 0.00 0.00 0.50 / ],
    lightgrey     => [ ++$colors, qw/ 0.00 0.00 0.00 0.40 / ],
    verylightgrey => [ ++$colors, qw/ 0.00 0.00 0.00 0.30 / ],
    palegrey      => [ ++$colors, qw/ 0.00 0.00 0.00 0.20 / ],
    verypalegrey  => [ ++$colors, qw/ 0.00 0.00 0.00 0.10 / ],
    white         => [ ++$colors, qw/ 0.00 0.00 0.00 0.00 / ],

    # pink                            
    verydeeppink  => [ ++$colors, qw/ 0.00 0.60 0.15 0.35 / ],
    deeppink      => [ ++$colors, qw/ 0.00 0.75 0.20 0.15 / ],
    verydarkpink  => [ ++$colors, qw/ 0.00 0.90 0.25 0.00 / ],
    darkpink      => [ ++$colors, qw/ 0.00 0.75 0.30 0.00 / ],
    pink          => [ ++$colors, qw/ 0.00 0.60 0.25 0.00 / ],
    lightpink     => [ ++$colors, qw/ 0.00 0.45 0.20 0.00 / ],
    verylightpink => [ ++$colors, qw/ 0.00 0.30 0.15 0.00 / ],
    palepink      => [ ++$colors, qw/ 0.00 0.20 0.10 0.00 / ],
    verypalepink  => [ ++$colors, qw/ 0.00 0.10 0.05 0.00 / ],

    # magenta                                 
    verydeepmagenta  => [ ++$colors, qw/ 0.00 0.60 0.00 0.40 / ],
    deepmagenta      => [ ++$colors, qw/ 0.00 0.75 0.00 0.20 / ],
    verydarkmagenta  => [ ++$colors, qw/ 0.00 0.90 0.00 0.00 / ],
    darkmagenta      => [ ++$colors, qw/ 0.00 0.75 0.00 0.00 / ],
    magenta          => [ ++$colors, qw/ 0.00 0.60 0.00 0.00 / ],
    lightmagenta     => [ ++$colors, qw/ 0.00 0.45 0.00 0.00 / ],
    verylightmagenta => [ ++$colors, qw/ 0.00 0.30 0.00 0.00 / ],
    palemagenta      => [ ++$colors, qw/ 0.00 0.20 0.00 0.00 / ],
    verypalemagenta  => [ ++$colors, qw/ 0.00 0.10 0.00 0.00 / ],

    # violet                                  
    verydeepviolet  => [ ++$colors, qw/ 0.30 0.60 0.00 0.40 / ],
    deepviolet      => [ ++$colors, qw/ 0.35 0.75 0.00 0.20 / ],
    verydarkviolet  => [ ++$colors, qw/ 0.45 0.90 0.00 0.00 / ],
    darkviolet      => [ ++$colors, qw/ 0.35 0.75 0.00 0.00 / ],
    violet          => [ ++$colors, qw/ 0.30 0.60 0.00 0.00 / ],
    lightviolet     => [ ++$colors, qw/ 0.22 0.45 0.00 0.00 / ],
    verylightviolet => [ ++$colors, qw/ 0.15 0.30 0.00 0.00 / ],
    paleviolet      => [ ++$colors, qw/ 0.10 0.20 0.00 0.00 / ],
    verypaleviolet  => [ ++$colors, qw/ 0.05 0.10 0.00 0.00 / ],

    # blue                            
    verydeepblue  => [ ++$colors, qw/ 0.60 0.60 0.00 0.40 / ],
    deepblue      => [ ++$colors, qw/ 0.75 0.75 0.00 0.20 / ],
    verydarkblue  => [ ++$colors, qw/ 0.90 0.90 0.00 0.00 / ],
    darkblue      => [ ++$colors, qw/ 0.75 0.75 0.00 0.00 / ],
    blue          => [ ++$colors, qw/ 0.60 0.60 0.00 0.00 / ],
    lightblue     => [ ++$colors, qw/ 0.45 0.45 0.00 0.00 / ],
    verylightblue => [ ++$colors, qw/ 0.30 0.30 0.00 0.00 / ],
    paleblue      => [ ++$colors, qw/ 0.20 0.20 0.00 0.00 / ],
    verypaleblue  => [ ++$colors, qw/ 0.10 0.10 0.00 0.00 / ],

    # skyblue                                 
    verydeepskyblue  => [ ++$colors, qw/ 0.60 0.40 0.00 0.40 / ],
    deepskyblue      => [ ++$colors, qw/ 0.75 0.50 0.00 0.20 / ],
    verydarkskyblue  => [ ++$colors, qw/ 0.90 0.60 0.00 0.00 / ],
    darkskyblue      => [ ++$colors, qw/ 0.75 0.50 0.00 0.00 / ],
    skyblue          => [ ++$colors, qw/ 0.60 0.40 0.00 0.00 / ],
    lightskyblue     => [ ++$colors, qw/ 0.45 0.30 0.00 0.00 / ],
    verylightskyblue => [ ++$colors, qw/ 0.30 0.22 0.00 0.00 / ],
    paleskyblue      => [ ++$colors, qw/ 0.20 0.14 0.00 0.00 / ],
    verypaleskyblue  => [ ++$colors, qw/ 0.10 0.07 0.00 0.00 / ],

    # cyan                            
    verydeepcyan  => [ ++$colors, qw/ 0.60 0.00 0.00 0.40 / ],
    deepcyan      => [ ++$colors, qw/ 0.75 0.00 0.00 0.20 / ],
    verydarkcyan  => [ ++$colors, qw/ 0.90 0.00 0.00 0.00 / ],
    darkcyan      => [ ++$colors, qw/ 0.75 0.00 0.00 0.00 / ],
    cyan          => [ ++$colors, qw/ 0.60 0.00 0.00 0.00 / ],
    lightcyan     => [ ++$colors, qw/ 0.45 0.00 0.00 0.00 / ],
    verylightcyan => [ ++$colors, qw/ 0.30 0.00 0.00 0.00 / ],
    palecyan      => [ ++$colors, qw/ 0.20 0.00 0.00 0.00 / ],
    verypalecyan  => [ ++$colors, qw/ 0.10 0.00 0.00 0.00 / ],

    # seagreen                        
    verydeepseagreen  => [ ++$colors, qw/ 0.60 0.00 0.15 0.40 / ],
    deepseagreen      => [ ++$colors, qw/ 0.75 0.00 0.20 0.20 / ],
    verydarkseagreen  => [ ++$colors, qw/ 0.90 0.00 0.25 0.00 / ],
    darkseagreen      => [ ++$colors, qw/ 0.75 0.00 0.30 0.00 / ],
    seagreen          => [ ++$colors, qw/ 0.60 0.00 0.25 0.00 / ],
    lightseagreen     => [ ++$colors, qw/ 0.45 0.00 0.20 0.00 / ],
    verylightseagreen => [ ++$colors, qw/ 0.30 0.00 0.15 0.00 / ],
    paleseagreen      => [ ++$colors, qw/ 0.20 0.00 0.10 0.00 / ],
    verypaleseagreen  => [ ++$colors, qw/ 0.10 0.00 0.05 0.00 / ],

    # green                           
    verydeepgreen  => [ ++$colors, qw/ 0.60 0.00 0.60 0.40 / ],
    deepgreen      => [ ++$colors, qw/ 0.75 0.00 0.75 0.20 / ],
    verydarkgreen  => [ ++$colors, qw/ 0.90 0.00 0.90 0.00 / ],
    darkgreen      => [ ++$colors, qw/ 0.75 0.00 0.75 0.00 / ],
    green          => [ ++$colors, qw/ 0.60 0.00 0.60 0.00 / ],
    lightgreen     => [ ++$colors, qw/ 0.45 0.00 0.45 0.00 / ],
    verylightgreen => [ ++$colors, qw/ 0.30 0.00 0.30 0.00 / ],
    palegreen      => [ ++$colors, qw/ 0.20 0.00 0.20 0.00 / ],
    verypalegreen  => [ ++$colors, qw/ 0.10 0.00 0.10 0.00 / ],

    # limegreen                       
    verydeeplimegreen  => [ ++$colors, qw/ 0.15 0.00 0.60 0.40 / ],
    deeplimegreen      => [ ++$colors, qw/ 0.20 0.00 0.75 0.20 / ],
    verydarklimegreen  => [ ++$colors, qw/ 0.25 0.00 0.90 0.00 / ],
    darklimegreen      => [ ++$colors, qw/ 0.30 0.00 0.75 0.00 / ],
    limegreen          => [ ++$colors, qw/ 0.25 0.00 0.60 0.00 / ],
    lightlimegreen     => [ ++$colors, qw/ 0.20 0.00 0.45 0.00 / ],
    verylightlimegreen => [ ++$colors, qw/ 0.15 0.00 0.30 0.00 / ],
    palelimegreen      => [ ++$colors, qw/ 0.10 0.00 0.20 0.00 / ],
    verypalelimegreen  => [ ++$colors, qw/ 0.05 0.00 0.10 0.00 / ],

    # yellow                                  
    verydeepyellow  => [ ++$colors, qw/ 0.00 0.00 0.60 0.40 / ],
    deepyellow      => [ ++$colors, qw/ 0.00 0.00 0.75 0.20 / ],
    verydarkyellow  => [ ++$colors, qw/ 0.00 0.00 0.90 0.00 / ],
    darkyellow      => [ ++$colors, qw/ 0.00 0.00 0.75 0.00 / ],
    yellow          => [ ++$colors, qw/ 0.00 0.00 0.60 0.00 / ],
    lightyellow     => [ ++$colors, qw/ 0.00 0.00 0.45 0.00 / ],
    verylightyellow => [ ++$colors, qw/ 0.00 0.00 0.30 0.00 / ],
    paleyellow      => [ ++$colors, qw/ 0.00 0.00 0.20 0.00 / ],
    verypaleyellow  => [ ++$colors, qw/ 0.00 0.00 0.10 0.00 / ],

    # orange                                  
    verydeeporange  => [ ++$colors, qw/ 0.00 0.30 0.80 0.40 / ],
    deeporange      => [ ++$colors, qw/ 0.00 0.32 0.85 0.20 / ],
    verydarkorange  => [ ++$colors, qw/ 0.00 0.35 0.90 0.00 / ],
    darkorange      => [ ++$colors, qw/ 0.00 0.30 0.75 0.00 / ],
    orange          => [ ++$colors, qw/ 0.00 0.25 0.60 0.00 / ],
    lightorange     => [ ++$colors, qw/ 0.00 0.20 0.45 0.00 / ],
    verylightorange => [ ++$colors, qw/ 0.00 0.15 0.30 0.00 / ],
    paleorange      => [ ++$colors, qw/ 0.00 0.10 0.20 0.00 / ],
    verypaleorange  => [ ++$colors, qw/ 0.00 0.05 0.10 0.00 / ],

    # red                                     
    verydeepred  => [ ++$colors, qw/ 0.00 0.60 0.60 0.40 / ],
    deepred      => [ ++$colors, qw/ 0.00 0.75 0.75 0.20 / ],
    verydarkred  => [ ++$colors, qw/ 0.00 0.90 0.90 0.00 / ],
    darkred      => [ ++$colors, qw/ 0.00 0.75 0.75 0.00 / ],
    red          => [ ++$colors, qw/ 0.00 0.60 0.60 0.00 / ],
    lightred     => [ ++$colors, qw/ 0.00 0.45 0.45 0.00 / ],
    verylightred => [ ++$colors, qw/ 0.00 0.30 0.30 0.00 / ],
    palered      => [ ++$colors, qw/ 0.00 0.20 0.20 0.00 / ],
    verypalered  => [ ++$colors, qw/ 0.00 0.10 0.10 0.00 / ],

    # brown                           
    verydeepbrown  => [ ++$colors, qw/ 0.35 0.85 0.90 0.35 / ],
    deepbrown      => [ ++$colors, qw/ 0.30 0.75 0.90 0.25 / ],
    verydarkbrown  => [ ++$colors, qw/ 0.25 0.70 0.90 0.15 / ],
    darkbrown      => [ ++$colors, qw/ 0.20 0.60 0.75 0.05 / ],
    brown          => [ ++$colors, qw/ 0.17 0.50 0.65 0.00 / ],
    lightbrown     => [ ++$colors, qw/ 0.15 0.40 0.55 0.00 / ],
    verylightbrown => [ ++$colors, qw/ 0.12 0.28 0.45 0.00 / ],
    palebrown      => [ ++$colors, qw/ 0.10 0.18 0.30 0.00 / ],
    verypalebrown  => [ ++$colors, qw/ 0.05 0.10 0.20 0.00 / ],
);    # %COLORS
%GRADIENTS = (
    'gradient' =>
'/gradient { dup 1 eq { bg 9 4 roll 1 add } if dup dup 1 sub exch 4 mul 3 add 1 roll 4 mul array astore exch dup dup 0 eq exch 1 eq or { 3 -1 roll mul 4 mul cvi 4 getinterval aload pop } { 3 -1 roll 2 copy 5 2 roll mul floor dup 4 1 roll 4 mul cvi 8 getinterval aload pop 9 dict begin /_ku X /_yu X /_mu X /_cu X /_kd X /_yd X /_md X /_cd X 1 exch div 3 -1 roll exch div exch sub /_sc X _cd _cu _cd sub _sc mul add _md _mu _md sub _sc mul add _yd _yu _yd sub _sc mul add _kd _ku _kd sub _sc mul add end } ie } B',
    'rainbow' =>
'/rainbow { dup (-) ne { dup 0.1 le { 0.2 div 0.5 add 1 0 0 } { dup 0.4 le { 1 1 3 -1 roll 0.1 sub 0.3 div sub 0 0 } { dup 0.6 le { 1 0 3 -1 roll 0.4 sub 0.2 div 0 } { dup 0.8 le { 1 exch 0.6 sub 0.2 div sub 0 1 0 } { 0 exch 0.8 sub 0.2 div 1 0 } ie } ie } ie } ie } { pop } ie } B',
);    # %GRADIENTS
$formats = 0;
%FORMATS = (    # [ FormatNUMBER, X(short edge), Y(long edge) ]
        # A Series
    a0  => [ ++$formats, 2384, 3370 ],
    a1  => [ ++$formats, 1684, 2384 ],
    a2  => [ ++$formats, 1190, 1684 ],
    a3  => [ ++$formats, 842,  1190 ],
    a4  => [ ++$formats, 595,  842 ],
    a5  => [ ++$formats, 420,  595 ],
    a6  => [ ++$formats, 297,  420 ],
    a7  => [ ++$formats, 210,  297 ],
    a8  => [ ++$formats, 148,  210 ],
    a9  => [ ++$formats, 105,  148 ],
    a10 => [ ++$formats, 73,   105 ],

    # B Series
    b0  => [ ++$formats, 2920, 4127 ],
    b1  => [ ++$formats, 2064, 2920 ],
    b2  => [ ++$formats, 1460, 2064 ],
    b3  => [ ++$formats, 1032, 1460 ],
    b4  => [ ++$formats, 729,  1032 ],
    b5  => [ ++$formats, 516,  729 ],
    b6  => [ ++$formats, 363,  516 ],
    b7  => [ ++$formats, 258,  363 ],
    b8  => [ ++$formats, 181,  258 ],
    b9  => [ ++$formats, 127,  181 ],
    b10 => [ ++$formats, 91,   127 ],

    # USA Formats
    executive => [ ++$formats, 540,  720 ],
    folio     => [ ++$formats, 612,  936 ],
    legal     => [ ++$formats, 612,  1008 ],
    letter    => [ ++$formats, 612,  792 ],
    quarto    => [ ++$formats, 610,  780 ],
    statement => [ ++$formats, 396,  612 ],
    '10x14'   => [ ++$formats, 720,  1008 ],
    ledger    => [ ++$formats, 1224, 792 ],
    tabloid   => [ ++$formats, 792,  1224 ],
);    # %FORMATS
      # Program status strings.
%Messages = (

    # ERROR MESSAGES
    FILE_NO_OPEN => $spl . $Warn
    . "Cannot Open Current file \"\%s\" . Not used !!!\n"
    . $spl,
    USER_HALT => $spl . $Warn
    . "$PROGRAM has been stopped by user !!!\n" . $spl . $Warn
    . "---------- Exiting NOW !!! ----------\n"
    . $spl,
    PROCESS_HALT => $spl . $Warn
    . "------- $PROGRAM is down !!! -------\n" . $spl . $Warn
    . "---------- Exiting NOW !!! ----------\n"
    . $spl,
    UNKNOWN_CL_OPTION => $Warn
    . "Error trapped while processing command-line:\n"
    . ( " " x 16 ) . "\%s\n",
    CMD_LINE_ERROR => $spl . $spw
    . " Please, check your command-line options!!!\n" . $Error . "\n" . $spw
    . " "
    . ( "." x 12 )
    . " Type \"$PROGRAM -h\" for help.\n"
    . $spl,
    SECTION_NOT_DEF => $Warn
    . "You probably forgot a section header, unable to parse this record.\n",
    VAR_NOT_DEFINED     => $Warn . "\%s variable not defined: \"\%s\" .\n",
    VARTYPE_NOT_DEFINED => $Error
    . "Variable type \"\%s\" not defined,\n" . $spw
    . "  could not check value for \"\%s\".\n",
    BAD_REGEXP => $Warn
    . "Ill-formed regular expression found in custom file:\n" . $spw
    . " ---> \%s <---\n",
    NOT_A_BOOLEAN => $Warn
    . "\"\%s\" variable requires a boolean value:\n" . $spw
    . "     (ON/OFF, 1/0, TRUE/FALSE, YES/NO)\n",
    NOT_ALPHA => $Warn
    . "\"\%s\" variable parameter must fit \"[a-zA-Z][a-zA-Z_0-9]*\"...\n",
    NOT_A_NATURAL  => $Warn . "\"\%s\" variable is not a positive integer.\n",
    NOT_AN_INTEGER => $Warn . "\"\%s\" variable is NOT an integer.\n",
    NOT_A_FLOAT    => $Warn . "\"\%s\" variable is NOT a decimal number.\n",
    NOT_A_REAL     => $Warn
    . "\"\%s\" variable is NOT a real number (with exponent).\n",
    NOT_A_COLOR => $Warn . "\"\%s\" color not defined in \"\%s\" variable.\n",
    NOT_A_GRADIENT => $Warn
    . "\"\%s\" gradient color wrongly defined in \"\%s\" variable.\n",
    NOTCMYKFLOAT => $Warn
    . "\"\%s\": value must be a float between 0 and 1, \"\%s\" not valid.\n",
    NOT_A_FONT => $Warn . "Sorry, \"\%s\" font is not defined for \"\%s\".\n",
    BBOX_NOTOK => $Warn
    . "CANNOT understand \"\%s\" set as \"\%s\" ...\n" . $spw
    . " Format must be: \"\%s=<width>,<height>\" \n",
    NOT_A_PAGE => $Warn
    . "\"\%s\" page-size is not defined for \"\%s\" variable.\n",
    NOT_A_NUMBER => $Warn
    . "\"\%s\" variable needs a number (with or without units).\n",
    NOT_A_DNAUNIT => $Warn
    . "\"\%s\" variable requires nucleotide units.\n" . $spw
    . " \"\%s\" is not valid (units must be in Gb, Mb, Kb, or bases).\n",
    NOT_A_PSUNIT => $Warn
    . "\"\%s\" variable requires PostScript units.\n" . $spw
    . " \"\%s\" is not valid (units must be in points, cm, mm, or inches).\n",
    RANGE_NOTOK => $Warn
    . "CANNOT understand \"\%s\" set as \"\%s\" ...\n" . $spw
    . " Format must be: \"\%s=<lower-value>..<upper-value>\" \n",
    LINESTY_NOTOK => $Warn
    . "LINE style \"\%s\" is not defined, setting to \"solid\"...\n",
    RIBBON_NOTOK => $Warn
    . "RIBBON style \"\%s\" is not defined, using \"\%s=none\" instead...\n",
    SHAPE_NOTOK => $Warn
    . "SHAPE \"\%s\" is not defined, using default shape...\n",
    GRP_SHAPE_NOTOK => $Warn
    . "GROUP SHAPE \"\%s\" is not defined, using defaults...\n",
    SPECIAL_NOT_DEF => $Warn
    . "Current special feature \"\%s\", NOT recognized, NOT set...\n",
    CANNOTDEFNOW => $Warn
    . "\%s \"\%s\", must be defined at the very beginning \n" . $spl
    . "of the first loaded custom file (within an EXTRA section, \"\# X \#\")...\n",
    WRONGSPECREC => $Warn
    . "A \%s special feature record requires \n" . $spl
    . "  at least \%s fields, and must follow this format:\n" . $spl
    . "  \"\%s\"\n",
    EMPTYSTRING => $Warn . "\"\%s\", contains an EMPTY STRING, skipping...\n",
    NOT_ENOUGH_FIELDS => $Warn
    . "Not enough fields in file \"\%s\", line \%s :\n\t\%s\n",
    ORI_GREATER_END => $Warn
    . "Start greater than end \"\%s > \%s\" in file \"\%s\" line \"\%s\".\n",
    STRAND_MISMATCH => $Warn
    . " Strand mismatch definition \"\%s\" in file \"\%s\" line \"\%s\".\n",
    FRAME_MISMATCH => $Warn
    . " Frame mismatch definition \"\%s\" in file \"\%s\" line \"\%s\".\n",
    NOT_A_STRAND => $Warn
    . "\"\%s\" STRAND not recognized on \"\%s\"... "
    . "Must be [+], [.] or [-]\n",

    # WORKING MESSAGES
    CHECKING_FILENAMES => $sp . "### Validating INPUT FILENAMES\n" . $sp,
    READING_FILE       => "###---> \"\%s\" exists, including as Input File.\n",
    READING_STDIN => "###---> Including GFF records from standard input.\n",
    CHECKING_CMDLN_OPTS => $sp . "### Checking COMMAND-LINE Options\n" . $sp,
    CHECKING_CMDLN_VARS => $sp
    . "### Checking Custom Variables SET by COMMAND-LINE\n"
    . $sp,
    CHECKING_CUSTOM_NAMES => $sp . "### Validating CUSTOM FILENAMES\n" . $sp,
    READING_FROM_PATH     =>
    "###---> Custom File NOT FOUND in local path: \"\%s\"\n"
    . "###     Trying to find in \"GFF2APLOT_CUSTOMPATH\": \"\%s\"\n",
    READING_CUSTOM_FILE =>
    "###---> \"\%s\" exists, including as Custom File.\n",
    NO_CUSTOM_FILES =>
    "###---> NO CUSTOM FILES found. Using program DEFAULTS.\n",
    CFILE_NOTOPEN => "###---> \"\%s\" NOT FOUND, skipping!!!\n",
    SETPAGESIZE   => "###---> Page size set to \"\%s\" ( \%s )...\n",
    SETPAGEAXES   =>
    "###---> Setting plot AXES: \n###     \%s\n###     \%s\n###     \%s\n",
    PLOTLIMITS => "###---> SEQUENCE BOUNDARIES selected for this plot:\n"
    . "###      \%s\n###      \%s\n",
    ONZOOM => "###---> SEQUENCE ZOOM that has been choosen:\n"
    . "###      \%s\n###      \%s\n",
    NOZOOM     => "###---> SEQUENCE ZOOM is not enabled for this plot...\n",
    ONZOOMAREA => "###---> ZOOM AREA that is going to be highlighted:\n"
    . "###      \%s\n###      \%s\n",
    NOZOOMAREA => "###---> NO ZOOM AREA was selected for this plot... \n",
    PLOTSCORES => "###---> SCORE BOUNDARIES selected for this plot:\n"
    . "###      \%s\n###      \%s\n",
    TICKDONE => "###---> TICKMARKS were set: \n"
    . "###      Nucleotide-Scale: \%s\n###      Score-Scale:      \%s\n",
    SCALEDONE => "###---> X/Y scale settings are: \n###     \%s\n"
    . "###     Are X/Y axes going to have same length ?  \%s \n",
    PSHEADER   => "###---> PostScript Header DONE...\n",
    PSCOLORS   => "###---> PostScript Colors Table SET...\n",
    PSFORMATS  => "###---> PostScript Page-formats Table SET...\n",
    PSVARS     => "###---> PostScript Variables SET...\n",
    PSCODE     => "###---> PostScript Functions SET...\n",
    PSPLOT     => "###---> Writting PostScript Page:\n",
    PSPLOTDONE => "###---> PostScript Page FINISHED...\n",
    BLOCKRB    => "###       + drawing \%s ribbons...\n",
    BLOCKRBL   => "###       + drawing \%s ribbons outline...\n",
    BLOCKXY    => "###       + drawing \%s annotation (\%s)...\n",
    BLOCKA     => "###       + drawing alignment features (\%s)...\n",
    BLOCKFX    => "###       + drawing special features from custom files...\n",
    BLOCKTODO  => "###     + \%s section :\n",
    BLOCKDONE  => "###       + \%s section DONE...\n",
    NOGFFDATA  => "###     ... \%s GFF elements NOT FOUND !!!\n",
    NOAPLOTDATA      => "###     ... Alignment elements NOT FOUND !!!\n",
    SHOW_VERSION     => $sp . "### \%s -- \%s\n" . $sp,
    READ_CUSTOM_FILE => $sp
    . "### Reading Customization Parameters from \"\%s\"\n"
    . $sp,
    NO_CUSTOM_FOUND => $sp
    . "### NO CUSTOM FILES found: Using program DEFAULTS.\n"
    . $sp,
    READ_GFF_FILE => $sp . "### Reading GFF records from \"\%s\"\n" . $sp,
    SORT_GFF      => "\%sSorting \%s\n",
    SORT_SEQ      => "\%sSequence: \%s\n",
    SORT_SRC      => "\%sSource: \%s\n",
    SORT_STR      => "\%sStrand: \%s\n",
    SORT_GRP      => "\%sGroup: \%s\n",
    SORT_FTR      => "\%s       Sorted \%s elements.\n",
    SORT_GPN      => "\%s\`Sorted \%s groups.\n",
);    # %Messages
my $total_time = 0;
my $DATE       = localtime;
my $USER       = defined( $ENV{USER} ) ? $ENV{USER} : 'Child Process';
my $PVER = sprintf( "v%vd", $^V );
%{ $USED{'COLORS'} } = (
    'white'         => $COLORS{'white'}->[0],
    'black'         => $COLORS{'black'}->[0],
    'lightred'      => $COLORS{'lightred'}->[0],
    'verylightgrey' => $COLORS{'verylightgrey'}->[0],
    'grey'          => $COLORS{'grey'}->[0],
    'red'           => $COLORS{'red'}->[0],
);
%{ $USED{'GRADIENTS'} } = ();
%{ $USED{'FORMATS'} }   = (
    'a4' => $FORMATS{'a4'}->[0],
);

#
# MAIN PROGRAM LOOP
#

# &set_default_vars;

%CmdLineVars = ();    # Reseting Command-Line OPTIONS
&parse_command_line;

%SpecialCache = (     # Reseting Special features cache
    PRE  => [],
    POST => [],
);
%CustomVars = ();    # Reseting Customization OPTIONS
&parse_custom_files;

%Vars = ();
&merge_custom_vars;

%GFF_DATA = %ALN_DATA = ();    # Reseting DATA
&parse_GFF_files;

&map_vars_data;

&sort_elements;

&set_page_vars;

&make_plot;

$total_time = &timing($T);
&header( "$PROGRAM HAS FINISHED", "[ $total_time ]" );

&close_logfile();
exit(0);

#
# MAIN FUNCTIONS
#
sub set_default_vars() {
    %DefaultVars = (
        LAYOUT => {    ## '# L #'
                # PAGE definition
            page_bbox => { TYPE => 'BBOX', VALUE => undef },
            page_size => { TYPE => 'PAGE', VALUE => 'a4' },

            # Page MARGINS
            margin_left   => { TYPE => 'PS_UNIT', VALUE => [ 1, 'cm' ] },
            margin_right  => { TYPE => 'PS_UNIT', VALUE => [ 1, 'cm' ] },
            margin_top    => { TYPE => 'PS_UNIT', VALUE => [ 1, 'cm' ] },
            margin_bottom => { TYPE => 'PS_UNIT', VALUE => [ 1, 'cm' ] },

            # Page COLORS                     
            background_color => { TYPE => 'COLOR', VALUE => 'white' },
            foreground_color => { TYPE => 'COLOR', VALUE => 'black' },

            # GLOBAL Labels                   
            show_title => { TYPE => 'BOOLEAN', VALUE => $T },
            title      => { TYPE => 'STRING',  VALUE => undef },             # T
            title_font => { TYPE => 'FONT',    VALUE => 'helvetica-bold' },
            title_fontsize => { TYPE => 'PS_UNIT', VALUE => [ 24, 'pt' ] },
            show_subtitle  => { TYPE => 'BOOLEAN', VALUE => $T },
            subtitle       => { TYPE => 'STRING',  VALUE => undef },        # ST
            subtitle_font  => { TYPE => 'FONT',    VALUE => 'helvetica' },
            subtitle_fontsize => { TYPE => 'PS_UNIT', VALUE => [ 16,   'pt' ] },
            show_x_label      => { TYPE => 'BOOLEAN', VALUE => $T },
            x_label           => { TYPE => 'STRING',  VALUE => undef },
            x_label_font => { TYPE => 'FONT', VALUE => 'helvetica-bold' },
            x_label_fontsize => { TYPE => 'PS_UNIT', VALUE => [ 12, 'pt' ] },
            show_y_label     => { TYPE => 'BOOLEAN', VALUE => $T },
            y_label          => { TYPE => 'STRING',  VALUE => undef },
            y_label_font => { TYPE => 'FONT', VALUE => 'helvetica-bold' },
            y_label_fontsize => { TYPE => 'PS_UNIT', VALUE => [ 12, 'pt' ] },
            show_percent_box_label => { TYPE => 'BOOLEAN', VALUE => $T },
            percent_box_label      => { TYPE => 'STRING',  VALUE => undef },
            percent_box_label_font =>
            { TYPE => 'FONT', VALUE => 'helvetica-bold' },
            percent_box_label_fontsize =>
            { TYPE => 'PS_UNIT', VALUE => [ 12, 'pt' ] },
            show_percent_box_sublabel => { TYPE => 'BOOLEAN', VALUE => $T },
            percent_box_sublabel      => { TYPE => 'STRING',  VALUE => undef },
            percent_box_sublabel_font =>
            { TYPE => 'FONT', VALUE => 'helvetica' },
            percent_box_sublabel_fontsize =>
            { TYPE => 'PS_UNIT', VALUE => [ 9, 'pt' ] },
            show_extra_box_label => { TYPE => 'BOOLEAN', VALUE => $T },
            extra_box_label      => { TYPE => 'STRING',  VALUE => undef },
            extra_box_label_font =>
            { TYPE => 'FONT', VALUE => 'helvetica-bold' },
            extra_box_label_fontsize =>
            { TYPE => 'PS_UNIT', VALUE => [ 10, 'pt' ] },
            show_extra_box_sublabel => { TYPE => 'BOOLEAN', VALUE => $T },
            extra_box_sublabel      => { TYPE => 'STRING',  VALUE => undef },
            extra_box_sublabel_font => { TYPE => 'FONT', VALUE => 'helvetica' },
            extra_box_sublabel_fontsize =>
            { TYPE => 'PS_UNIT', VALUE => [ 9, 'pt' ] },
            x_sequence_coords     => { TYPE => 'SEQRANGE', VALUE => undef },
            x_sequence_start      => { TYPE => 'SEQ_UNIT', VALUE => undef },
            x_sequence_end        => { TYPE => 'SEQ_UNIT', VALUE => undef },
            y_sequence_coords     => { TYPE => 'SEQRANGE', VALUE => undef },
            y_sequence_start      => { TYPE => 'SEQ_UNIT', VALUE => undef },
            y_sequence_end        => { TYPE => 'SEQ_UNIT', VALUE => undef },
            x_sequence_zoom       => { TYPE => 'SEQRANGE', VALUE => undef },
            x_sequence_zoom_start => { TYPE => 'SEQ_UNIT', VALUE => undef },
            x_sequence_zoom_end   => { TYPE => 'SEQ_UNIT', VALUE => undef },
            y_sequence_zoom       => { TYPE => 'SEQRANGE', VALUE => undef },
            y_sequence_zoom_start => { TYPE => 'SEQ_UNIT', VALUE => undef },
            y_sequence_zoom_end   => { TYPE => 'SEQ_UNIT', VALUE => undef },
            zoom                  => { TYPE => 'BOOLEAN',  VALUE => $F },
            zoom_area             => { TYPE => 'BOOLEAN',  VALUE => $F },
            zoom_marks            => { TYPE => 'BOOLEAN',  VALUE => $F },
            zoom_area_mark_width => { TYPE => 'PS_UNIT', VALUE => [ 2, 'pt' ] },
            zoom_area_mark_style => { TYPE => 'LINESTY', VALUE => 'solid' },
            zoom_area_mark_color => { TYPE => 'COLOR',   VALUE => 'lightred' },
            zoom_area_fill_color => { TYPE => 'COLOR',   VALUE => undef },
            alignment_name       => { TYPE => 'STRING',  VALUE => undef },
            x_sequence_name      => { TYPE => 'STRING',  VALUE => undef },
            y_sequence_name      => { TYPE => 'STRING',  VALUE => undef },
            aplot_xy_same_length => { TYPE => 'BOOLEAN', VALUE => $T },
            aplot_xy_scale       => { TYPE => 'FLOAT',   VALUE => undef },
            alignment_scale_width => { TYPE => 'BOOLEAN', VALUE => $F },
            alignment_scale_color => { TYPE => 'BOOLEAN', VALUE => $F },
            show_ribbons          => { TYPE => 'BOOLEAN', VALUE => undef },
            ribbon_color_merge    => { TYPE => 'BOOLEAN', VALUE => $F },
            color_merge_factor    => { TYPE => 'FLOAT',   VALUE => 0.5 },
            ribbon_style          => { TYPE => 'RIBBON',  VALUE => undef },
            ribbon_color          => { TYPE => 'COLOR',   VALUE => undef },
            show_grid             => { TYPE => 'BOOLEAN', VALUE => $F },
            show_percent_box_grid => { TYPE => 'BOOLEAN', VALUE => $T },
            show_percent_box      => { TYPE => 'BOOLEAN', VALUE => $F },
            show_extra_box        => { TYPE => 'BOOLEAN', VALUE => $F },
            aplot_box_bgcolor     => { TYPE => 'COLOR',   VALUE => 'bg' },
            percent_box_bgcolor   => { TYPE => 'COLOR',   VALUE => 'bg' },
            extra_box_bgcolor     => { TYPE => 'COLOR',   VALUE => 'bg' },
            percent_box_height => { TYPE => 'PS_UNIT', VALUE => [ -1, 'pt' ] },
            extra_box_height    => { TYPE => 'PS_UNIT', VALUE => [ -1, 'pt' ] },
            show_tickmark_label => { TYPE => 'BOOLEAN', VALUE => $T },
            show_only_bottom_ticks    => { TYPE => 'BOOLEAN',  VALUE => $F },
            show_aplot_x_ticks        => { TYPE => 'BOOLEAN',  VALUE => $T },
            show_aplot_y_ticks        => { TYPE => 'BOOLEAN',  VALUE => $T },
            show_percent_x_ticks      => { TYPE => 'BOOLEAN',  VALUE => $T },
            show_percent_y_ticks      => { TYPE => 'BOOLEAN',  VALUE => $T },
            show_extrabox_x_ticks     => { TYPE => 'BOOLEAN',  VALUE => $T },
            show_extrabox_y_ticks     => { TYPE => 'BOOLEAN',  VALUE => $T },
            aplot_major_tickmark      => { TYPE => 'INTEGER',  VALUE => undef },
            aplot_minor_tickmark      => { TYPE => 'INTEGER',  VALUE => 5 },
            aplot_score_range         => { TYPE => 'RANGE',    VALUE => undef },
            percent_major_tickmark    => { TYPE => 'INTEGER',  VALUE => undef },
            percent_minor_tickmark    => { TYPE => 'INTEGER',  VALUE => 4 },
            percent_box_score_range   => { TYPE => 'RANGE',    VALUE => undef },
            extra_major_tickmark      => { TYPE => 'INTEGER',  VALUE => 2 },
            extra_minor_tickmark      => { TYPE => 'INTEGER',  VALUE => 5 },
            extra_box_score_range     => { TYPE => 'RANGE',    VALUE => undef },
            major_tickmark_nucleotide => { TYPE => 'SEQ_UNIT', VALUE => undef },
            minor_tickmark_nucleotide => { TYPE => 'SEQ_UNIT', VALUE => undef },
            major_tickmark_score      => { TYPE => 'FLOAT',    VALUE => undef },
            minor_tickmark_score      => { TYPE => 'FLOAT',    VALUE => undef },
            show_feature_label        => { TYPE => 'BOOLEAN',  VALUE => $F },
            feature_label_font => { TYPE => 'FONT', VALUE => 'helvetica' },
            feature_label_fontsize =>
            { TYPE => 'PS_UNIT', VALUE => [ 6, 'pt' ] },
            feature_x_label_length => { TYPE => 'INTEGER', VALUE => undef },
            feature_x_label_angle  => { TYPE => 'INTEGER', VALUE => 0 },
            feature_x_label_rotate => { TYPE => 'BOOLEAN', VALUE => $F },
            feature_y_label_length => { TYPE => 'INTEGER', VALUE => undef },
            feature_y_label_angle  => { TYPE => 'INTEGER', VALUE => 0 },
            feature_y_label_rotate => { TYPE => 'BOOLEAN', VALUE => $T },
            show_group_label       => { TYPE => 'BOOLEAN', VALUE => $T },
            group_label_font     => { TYPE => 'FONT',    VALUE => 'helvetica' },
            group_label_fontsize => { TYPE => 'PS_UNIT', VALUE => [ 8, 'pt' ] },
            group_x_label_length => { TYPE => 'INTEGER', VALUE => undef },
            group_x_label_angle  => { TYPE => 'INTEGER', VALUE => 0 },
            group_x_label_rotate => { TYPE => 'BOOLEAN', VALUE => $F },
            group_y_label_length => { TYPE => 'INTEGER', VALUE => undef },
            group_y_label_angle  => { TYPE => 'INTEGER', VALUE => 0 },
            group_y_label_rotate => { TYPE => 'BOOLEAN', VALUE => $T },
            show_ps_warnings     => { TYPE => 'BOOLEAN', VALUE => $T },
            hide_credits         => { TYPE => 'BOOLEAN', VALUE => $F },
            align_tag            => { TYPE => 'ALPHA',   VALUE => 'target' },
            vector_tag           => { TYPE => 'ALPHA',   VALUE => 'vector' },
            label_tag            => { TYPE => 'ALPHA',   VALUE => 'id' },
            start_tag            => { TYPE => 'ALPHA',   VALUE => 'start' },
            end_tag              => { TYPE => 'ALPHA',   VALUE => 'end' },
            score_tag            => { TYPE => 'ALPHA',   VALUE => 'e_value' },
            strand_tag           => { TYPE => 'ALPHA',   VALUE => 'strand' },
            frame_tag            => { TYPE => 'ALPHA',   VALUE => 'frame' },
            scotype_tag          => { TYPE => 'ALPHA',   VALUE => 'vectype' },
            window_tag           => { TYPE => 'ALPHA',   VALUE => 'window' },
            step_tag             => { TYPE => 'ALPHA',   VALUE => 'step' },
            scoary_tag           => { TYPE => 'ALPHA',   VALUE => 'scores' },
        },
        SEQUENCE => {    ## '# Q #'
            feature_color => { TYPE => 'EXTCOLOR', VALUE => undef },
            ribbon_color  => { TYPE => 'COLOR',    VALUE => undef },
        },
        SOURCE => {    ## '# S #'
            hide          => { TYPE => 'BOOLEAN',  VALUE => $F },
            source_layer  => { TYPE => 'INTEGER',  VALUE => undef },
            feature_color => { TYPE => 'EXTCOLOR', VALUE => undef },
            ribbon_color  => { TYPE => 'COLOR',    VALUE => undef },
        },
        STRAND => {    ## '# T #'
            hide          => { TYPE => 'BOOLEAN',  VALUE => $F },
            strand_layer  => { TYPE => 'INTEGER',  VALUE => undef },
            feature_color => { TYPE => 'EXTCOLOR', VALUE => undef },
            ribbon_color  => { TYPE => 'COLOR',    VALUE => undef },
        },
        GROUP => {    ## '# G #'
            hide              => { TYPE => 'BOOLEAN', VALUE => $F },
            group_color       => { TYPE => 'COLOR',   VALUE => 'fg' },
            group_shape       => { TYPE => 'GPSHAPE', VALUE => 'bracket' },
            show_group_limits => { TYPE => 'BOOLEAN', VALUE => $F },
            group_label       => { TYPE => 'STRING',  VALUE => undef },
            show_group_label  => { TYPE => 'BOOLEAN', VALUE => $T },

# show_group_rule            => { TYPE => 'BOOLEAN', VALUE => $T     },
# show_group_arrow           => { TYPE => 'BOOLEAN', VALUE => $T     },
# feature_arrows_color       => { TYPE => 'COLOR'  , VALUE => 'fg'   },
            Show_JOINS       => { TYPE => 'BOOLEAN',  VALUE => $T },
            Join_Lines_COLOR => { TYPE => 'COLOR',    VALUE => 'fg' },
            feature_color    => { TYPE => 'EXTCOLOR', VALUE => undef },
            ribbon_color     => { TYPE => 'COLOR',    VALUE => undef },
        },
        FEATURE => {    ## '# F #'
            group_layer     => { TYPE => 'INTEGER',  VALUE => undef },
            hide            => { TYPE => 'BOOLEAN',  VALUE => $F },
            feature_color   => { TYPE => 'EXTCOLOR', VALUE => 'lightred' },
            alignment_color => { TYPE => 'EXTCOLOR', VALUE => 'fg' },
            ribbon_color    => { TYPE => 'COLOR',    VALUE => 'verylightgrey' },
            feature_shape   => { TYPE => 'SHAPE',    VALUE => 'box' },
            feature_label   => { TYPE => 'STRING',   VALUE => undef },
            show_feature_label    => { TYPE => 'BOOLEAN', VALUE => $F },
            show_ribbons          => { TYPE => 'BOOLEAN', VALUE => $F },
            ribbon_style          => { TYPE => 'RIBBON',  VALUE => 'none' },
            feature_scale_height  => { TYPE => 'BOOLEAN', VALUE => $T },
            feature_scale_color   => { TYPE => 'BOOLEAN', VALUE => $F },
            alignment_scale_width => { TYPE => 'BOOLEAN', VALUE => $F },
            alignment_scale_color => { TYPE => 'BOOLEAN', VALUE => $F },
            feature_layer         => { TYPE => 'INTEGER', VALUE => undef },

# Show_HalfHeightBOX         => { TYPE => 'BOOLEAN', VALUE => $T     },
# HalfSizeBox_BGCOLOR        => { TYPE => 'COLOR'  , VALUE => 'DEFAULT' },
# Show_FullHeightBOX         => { TYPE => 'BOOLEAN', VALUE => $T     },
# FullSizeBox_BGCOLOR        => { TYPE => 'COLOR'  , VALUE => 'DEFAULT' },
# Show_BOX_LABEL             => { TYPE => 'BOOLEAN', VALUE => $T     },
# Show_UserDef_BOX_LABEL     => { TYPE => 'BOOLEAN', VALUE => $T     },
# Show_RIBBON                => { TYPE => 'BOOLEAN', VALUE => $T     },
# Ribbon_BGCOLOR             => { TYPE => 'COLOR'  , VALUE => 'DEFAULT' },
            Show_GFF              => { TYPE => 'BOOLEAN', VALUE => $F },
            Show_GFF_ReverseOrder => { TYPE => 'BOOLEAN', VALUE => $F },
            Show_FUNCTION         => { TYPE => 'BOOLEAN', VALUE => $F },

# APlotLine_GroupScore       => { TYPE => 'BOOLEAN', VALUE => $F     },
# APlotLine_ScaleWidth       => { TYPE => 'BOOLEAN', VALUE => $F     },
# APlotLine_ScaleGrey        => { TYPE => 'BOOLEAN', VALUE => $F     },
            Show_SELECTION_BOX   => { TYPE => 'BOOLEAN', VALUE => $T },
            SelectionBox_BGCOLOR => { TYPE => 'COLOR',   VALUE => 'grey' },
            Function_COLOR       => { TYPE => 'COLOR',   VALUE => 'red' },
        },
    );    # %DefaultVars
    %Defaults = ();
    foreach my $sct ( keys %DefaultVars ) {
        foreach my $vnm ( keys %{ $DefaultVars{$sct} } ) {
            $Defaults{$sct}{$vnm} = $DefaultVars{$sct}{$vnm}{'VALUE'};
        };    # foreach $vnm
    };    # foreach $sct
    print
      LOGFILE ( Data::Dumper->Dump( [ \%DefaultVars ], [qw( *DefaultVars )] ) )
      if ( $LogFile && $Debug );
}    # set_default_vars

sub parse_command_line() {
    my $cmdln_stdin = undef;
    for ( my $a = 0 ; $a <= $#ARGV ; $a++ ) {
        next unless $ARGV[$a] =~ /^-$/o;
        $cmdln_stdin = $a - $#ARGV;
        splice( @ARGV, $a, 1 );
    }

    $SIG{__WARN__} = sub { &warn( 'UNKNOWN_CL_OPTION', $T, $_[0] ) };
    GetOptions(
        "v|verbose"         => \$Verbose,
        "V|logs-filename=s" => \$logs_filename,    # Print_Report -> LogFile
        "q|quiet"           => \$Quiet,
        "P|page-bbox"          => \$CmdLineVars{LAYOUT}{"page_bbox"},
        "p|page-size"          => \$CmdLineVars{LAYOUT}{"page_size"},
        "margin-left=s"        => \$CmdLineVars{LAYOUT}{"margin_left"},
        "margin-right=s"       => \$CmdLineVars{LAYOUT}{"margin_right"},
        "margin-top=s"         => \$CmdLineVars{LAYOUT}{"margin_top"},
        "margin-bottom=s"      => \$CmdLineVars{LAYOUT}{"margin_bottom"},
        "B|background-color=s" =>
        \$CmdLineVars{LAYOUT}{"background_color"},    # BACKGROUND_COLOR
        "F|foreground-color=s" =>
        \$CmdLineVars{LAYOUT}{"foreground_color"},    # FOREGROUND_COLOR
        "T|title=s"    => \$CmdLineVars{LAYOUT}{"title"},       # TITLE
        "t|subtitle=s" => \$CmdLineVars{LAYOUT}{"subtitle"},    # SUBTITLE
        "X|x-label=s"  => \$CmdLineVars{LAYOUT}{"x_label"},
        "Y|y-label=s"  => \$CmdLineVars{LAYOUT}{"y_label"},
        "L|percent-box-label=s" => \$CmdLineVars{LAYOUT}{"percent_box_label"},
        "l|extra-box-label=s"   => \$CmdLineVars{LAYOUT}{"extra_box_label"},
        "x|x-sequence-coords=s" => \$CmdLineVars{LAYOUT}{"x_sequence_coords"},
        "S|start-x-sequence=i"  => \$CmdLineVars{LAYOUT}
        {"x_sequence_start"},    # SEQUENCE1_ORIGIN # Zoom_SEQUENCE1_ORIGIN
        "E|end-x-sequence=i" => \$CmdLineVars{LAYOUT}
        {"x_sequence_end"},      # SEQUENCE1_END    # Zoom_SEQUENCE1_END
        "y|y-sequence-coords=s" => \$CmdLineVars{LAYOUT}{"y_sequence_coords"},
        "s|start-y-sequence=i"  => \$CmdLineVars{LAYOUT}
        {"y_sequence_start"},    # SEQUENCE2_ORIGIN # Zoom_SEQUENCE2_ORIGIN
        "e|end-y-sequence=i" => \$CmdLineVars{LAYOUT}
        {"y_sequence_end"},      # SEQUENCE2_END    # Zoom_SEQUENCE2_END
        "x-sequence-zoom=s"  => \$CmdLineVars{LAYOUT}{"x_sequence_zoom"},
        "y-sequence-zoom=s"  => \$CmdLineVars{LAYOUT}{"y_sequence_zoom"},
        "Z|zoom"             => \$CmdLineVars{LAYOUT}{"zoom"},
        "z|zoom-area"        => \$CmdLineVars{LAYOUT}{"zoom_area"},
        "A|alignment-name=s" =>
        \$CmdLineVars{LAYOUT}{"alignment_name"},    # Align_NAME
        "N|x-sequence-name=s" =>
        \$CmdLineVars{LAYOUT}{"x_sequence_name"},    # X-Sequence_NAME
        "n|y-sequence-name=s" =>
        \$CmdLineVars{LAYOUT}{"y_sequence_name"},    # Y-Sequence_NAME
        "r|aplot-xy-noteq" =>
        sub { $CmdLineVars{LAYOUT}{"aplot_xy_same_length"} = $F },
        "R|aplot-xy-scale=f" => \$CmdLineVars{LAYOUT}{"aplot_xy_scale"},
        "W|aln-scale-width"  => \$CmdLineVars{LAYOUT}
        {"alignment_scale_width"},  # APlotLine_ScaleWidth; APlotLine_GroupScore
        "w|aln-scale-color" => \$CmdLineVars{LAYOUT}
        {"alignment_scale_color"},   # APlotLine_ScaleGrey; APlotLine_GroupScore
        "K|show-ribbons=s" => sub {
            $CmdLineVars{LAYOUT}{"show_ribbons"} = $T;
            $CmdLineVars{LAYOUT}{"ribbon_style"} = $_[1];
        },
        "G|show-grid"        => \$CmdLineVars{LAYOUT}{"show_grid"},
        "g|hide-grid"        => sub { $CmdLineVars{LAYOUT}{"show_grid"} = $F },
        "I|show-percent-box" =>
        \$CmdLineVars{LAYOUT}{"show_percent_box"},    # Display_PERCENT-BOX
        "i|hide-percent-box" =>
        sub { $CmdLineVars{LAYOUT}{"show_percent_box"} = $F }
        ,    # Display_PERCENT-BOX
        "O|show-extra-box" =>
        \$CmdLineVars{LAYOUT}{"show_extra_box"},    # Display_EXTRA-BOX
        "o|hide-extra-box" =>
        sub { $CmdLineVars{LAYOUT}{"show_extra_box"} = $F }, # Display_EXTRA-BOX
        "D|aplot-box-color=s" =>
        \$CmdLineVars{LAYOUT}{"aplot_box_bgcolor"},          # APlotBox_BqGCOLOR
        "d|percent-box-color=s" =>
        \$CmdLineVars{LAYOUT}{"percent_box_bgcolor"},    # PercentBox_BGCOLOR
        "b|extra-box-color=s" =>
        \$CmdLineVars{LAYOUT}{"extra_box_bgcolor"},      # ExtraBox_BGCOLOR
        "nopswarnings" => sub { $CmdLineVars{LAYOUT}{"show_ps_warnings"} = $F },
        "a|hide-credits" =>
        \$CmdLineVars{LAYOUT}{"hide_credits"},           # Show_Credits
        "debug" => \$Debug,                              # Dumps Vars -> LogFile
        "layout-var=s@"        => \@{ $CmdLineVars{VARS}{LAYOUT} },
        "sequence-var=s@"      => \@{ $CmdLineVars{VARS}{SEQUENCE} },
        "source-var=s@"        => \@{ $CmdLineVars{VARS}{SOURCE} },
        "strand-var=s@"        => \@{ $CmdLineVars{VARS}{STRAND} },
        "group-var=s@"         => \@{ $CmdLineVars{VARS}{GROUP} },
        "feature-var=s@"       => \@{ $CmdLineVars{VARS}{FEATURE} },
        "C|custom-filename=s@" => \@custom_files,
        "version"              => \&prt_version,
        "h|help|?"             => \&prt_help,
    ) || ( &warn( 'CMD_LINE_ERROR', $T ), exit(1) );
    $SIG{__WARN__} = 'DEFAULT';

    CHKLOG:
    ( defined($logs_filename) ) && do {
        open( LOGFILE, "> " . $logs_filename )
          || ( &warn( 'FILE_NO_OPEN', $T, $logs_filename ), last CHKLOG );
        $LogFile = 1;
    };

    &header(
        '',            "RUNNING $PROGRAM", '', "User: $USER",
        "Date: $DATE", "Perl: $PVER"
    );

    &header("SETTING DEFAULTS");
    %DefaultVars = ();
    &set_default_vars;

    &header("CHECKING COMMAND-LINE OPTIONS");
    @data_files = ();
    &set_input_file($cmdln_stdin);
    @ARGV = ();    # ensuring that command-line ARGVs array is empty

    &set_custom_files();

    &check_command_line_vars();

    &footer("COMMAND-LINE CHECKED");
}    # parse_command_line

sub prt_version() {
    my $comment = $Messages{'SHOW_VERSION'};
    $comment = sprintf( $comment, $PROGRAM, $VERSION );
    &prt_to_stderr($comment);
    exit(1);
}    # prt_version

sub set_input_file() {
    my $stdin_flg = $F;
    my $chk_stdin = shift @_;
    my $t         = scalar(@ARGV);
    defined($chk_stdin) && do {
        abs($chk_stdin) > $t && ( $chk_stdin = -$t );
        $chk_stdin > 0 && ( $chk_stdin = 0 );
        $t += $chk_stdin;
        splice( @ARGV, $t, 0, '-' );
    };
    &report("CHECKING_FILENAMES");
    FILECHK: foreach my $test_file (@ARGV) {
        $test_file ne '-' && do {
            -e $test_file || do {
                &warn( 'FILE_NO_OPEN', $T, $test_file );
                next FILECHK;
            };
            &report( 'READING_FILE', $test_file );
            push @data_files, $test_file;
            next FILECHK;
        };
        $stdin_flg = $T;
        push @data_files, '-';
    };    # foreach
    scalar(@data_files) == 0 && do {
        push @data_files, '-';
        $stdin_flg = $T;
    };
    $stdin_flg && &report('READING_STDIN');
}    # set_input_file

sub check_command_line_vars() {
    &report("CHECKING_CMDLN_OPTS");
    &command_line_for_layout();
    defined( $CmdLineVars{VARS} ) && &report("CHECKING_CMDLN_VARS");
    &command_line_for_vars();
    &footer("COMMAND-LINE OPTIONS CHECKED");
}    # check_command_line_vars

sub command_line_for_layout() {
    ( $n, $c ) = ( 0, '' );
    foreach my $cvar ( keys %{ $CmdLineVars{LAYOUT} } ) {
        defined( $CmdLineVars{LAYOUT}{$cvar} ) && do {
            my @cary = ( 'LAYOUT', $cvar, $CmdLineVars{LAYOUT}{$cvar} );
            $c =
              &varscheck( $F, 'LAYOUT', \@cary, \%CmdLineVars ) ? '.' : $noCV;
            &counter( ++$n, $c );
        };    # defined($CmdLineVars{LAYOUT}{$cvar})
    };    # foreach $cvar
    &counter_end( $n, $c );
    print LOGFILE ( Data::Dumper->Dump( [ \%{ $CmdLineVars{LAYOUT} } ],
          [qw( *{$CmdLineVars{LAYOUT}} )] ) )
      if ( $LogFile && $Debug );
}    # command_line_for_layout

sub command_line_for_vars() {
    ( $n, $c ) = ( 0, '' );
    foreach my $calias (qw/ L Q S T G F /) {
        my ( $cvar, $cflg );
        $cvar = $VarKeys{$calias};
        $cflg = ( $calias ne 'L' ) ? $T : $F;
        defined( $CmdLineVars{VARS}{$cvar} ) && do {
            foreach my $ccv ( @{ $CmdLineVars{VARS}{$cvar} } ) {
                my @clin = ();
                TWOTHREE: {
                    $cflg && do {
                        $ccv =~ /^(.*?):{2}(.*?)={1}(.*?)$/o
                          && ( @clin = ( $1, $2, $3 ) );
                        last TWOTHREE;
                    };
                    $ccv =~ /^(.*?)={1}(.*?)$/o
                      && ( @clin = ( $cvar, $1, $2 ) );
                };    # TWOTHREE
                $c =
                  &varscheck( $cflg, $cvar, \@clin, \%CmdLineVars ) 
                  ? $calias
                  : $noCV;
                &counter( ++$n, $c );
            };    # foreach $ccv
        };    # $CmdLineVars{VARS}{$cvar}
    };    # foreach $calias
    &counter_end( $n, $c );
    print LOGFILE ( Data::Dumper->Dump( [ \%{ $CmdLineVars{VARS} } ],
          [qw( *{$CmdLineVars{VARS}} )] ) )
      if ( $LogFile && $Debug );
}    # command_line_for_vars

sub set_custom_files() {
    unshift @custom_files, $Custom_file;
    my @files = ();
    &report("CHECKING_CUSTOM_NAMES");
    MLOOP: foreach my $test_file (@custom_files) {
        FILECHK: {
            -e $test_file && last FILECHK;
            ( $test_file =~ m{/}og || $Custom_path eq '.' ) || do {
                my $tmpfl = $test_file;
                $test_file = "$Custom_path/$test_file";
                &report( 'READING_FROM_PATH', $tmpfl, $test_file );
                -e $test_file && last FILECHK;
            };
            scalar(@custom_files) == 1 && do {
                &report('NO_CUSTOM_FILES');
                last MLOOP;
            };
            &report( 'CFILE_NOTOPEN', $test_file );
            next MLOOP;
        };    # FILECHK
        &report( 'READING_CUSTOM_FILE', $test_file );
        push @files, $test_file;
    };    # MLOOP: foreach
    @custom_files = @files;
}    # set_custom_files

sub parse_custom_files() {
    &header("READING CUSTOM FILES");
    MAIN: {
        scalar(@custom_files) == 0 && do {
            &report( 'NO_CUSTOM_FOUND', $file );
            last MAIN;
        };
        LOAD: foreach $file (@custom_files) {
            open( THIS, "< $file" )
              || ( &warn( 'FILE_NO_OPEN', $T, $file ), next LOAD );
            &report( 'READ_CUSTOM_FILE', $file );
            ( $n, $c ) = ( 0, undef );
            my ( @line, $main, $_c, $_v, $v_flag );
            while (<THIS>) {
                /^\#/o && do {
                    /^\# ([XLQSTGF]) \#/o && do {
                        $_c = $1;
                        $c  = '*';
                        $_c ne 'X' && ( $ExtraCustomFlg = $F );
                        $v_flag = ( $_c ne 'L' ) ? $T : $F;
                        $_v = $VarKeys{$_c};
                        next;
                    };
                    $c = '.';
                    next;
                };
                ( $c = '.', next ) if /^\s*$/o;
                chomp;
                ( $main, undef ) = split /\b\s+\#/o;
                $_c eq 'X' && do {
                    $c =
                      &parse_special_customs( $main, $ExtraCustomFlg ) 
                      ? $_c
                      : $noCV;
                    next;
                };
                TWOTHREE: {
                    $v_flag && do {
                        $main =~ /^(.*?):{2}(.*?)={1}(.*?)$/o
                          && ( @line = ( $1, $2, $3 ) );
                        last TWOTHREE;
                    };
                    $main =~ /^(.*?)={1}(.*?)$/o && ( @line = ( $_v, $1, $2 ) );
                };    # TWOTHREE
                $c =
                  &varscheck( $v_flag, $_v, \@line, \%CustomVars ) 
                  ? $_c
                  : $noCV;
            }
            continue {
                &counter( ++$n, $c );
            };    # WHILE
            &counter_end( $n, $c );
            close(THIS);
        };    # LOAD
    };    # MAIN
    print LOGFILE (
      Data::Dumper->Dump(
          [ \%SpecialCache, \%CustomVars ],
          [qw( *SpecialCache   *CustomVars )]
      ) )
      if ( $LogFile && $Debug );
    &footer("CUSTOM FILES LOADED");
}    # parse_custom_files

sub varscheck() {
    my ( $flag, $class, $rec, $varec ) = @_;
    defined( $DefaultVars{$class} )
      || ( &warn( 'SECTION_NOT_DEF', $F ), return $F );
    my $_var = \%{ $DefaultVars{$class} };
    defined( $$rec[2] ) || return $F;
    defined( $_var->{ $$rec[1] } )
      || ( &warn( 'VAR_NOT_DEFINED', $F, $class, $$rec[1] ), return $F );
    &checkvarvalues( $_var->{ $$rec[1] }{'TYPE'}, \$$rec[2], $$rec[1] )
      || return $F;
    $flag && do {
        my @tmpary = ();
        defined( @{ $varec->{$class}{ $$rec[1] } } ) || do {
            @{ $varec->{$class}{ $$rec[1] } } = ();
        };
        @tmpary = &find_regexp( $$rec[0] );
        $tmpary[3] = lc( $tmpary[3] ) if ( $class eq 'FEATURE' );
        ( shift @tmpary ) || return $F;
        push @{ $varec->{$class}{ $$rec[1] } }, ( @tmpary, $$rec[2] );
        return $T;
    };
    $varec->{$class}{ $$rec[1] } = $$rec[2];
    return $T;
}    # varscheck

sub checkvarvalues() {
    my ( $_test, $_val, $_var ) = @_;
    my $_t = lc($$_val);
    defined( $TESTS{$_test} ) || do {
        &warn( 'VARTYPE_NOT_DEFINED', $T, $_test, $_var );
        return $F;
    };    # !defined($TESTS->{$_test})
    return $TESTS{$_test}->( $_val, $_var, $_t );
}    # checkvarvalues

sub find_regexp() {
    my $string = $_[0];
    my ( $isOK_flg, $not_flg, $id_flg, $tmpstr, $tmpid, $escflg );
    $isOK_flg = $T;
    $not_flg  = $F;
    $id_flg   = undef;
    $string =~ s{^!}{}o && ( $not_flg = $T );    # not_regexp is true
    $string =~ s{(\\@)$}{@@}o;    # scaping trailing '@' (if element ends with)
    ( $string, $escflg ) = &escape_input($string);
    ( $tmpstr, $tmpid ) = ( undef, undef );
    ( reverse($string) =~ m{^([^\/@]*?)(?:@){1}(.*)$}o ) && do {
        $tmpstr = reverse($2);
        $tmpid  = reverse($1);
    };    # reverse($string)
    ( defined($tmpid) && $tmpid ne "" ) && ( $id_flg = $tmpid );
    ( defined($tmpstr) && $tmpstr ne "" ) || do {
        $string eq '@' && ( $isOK_flg = $F );
        $tmpstr = $string;
    };    # (defined($tmpstr) && $tmpstr ne "")
    REGEXPS: {
        $tmpstr eq '*' && do {
            $string = '^.*$';    #'
            last REGEXPS;
        };    # $string eq '*'
        $tmpstr =~ m{^/(.*)/$}o && do {
            ( $string, $isOK_flg ) = &eval_regexp($1);
            last REGEXPS;
        };    # $tmpstr is a regexp
              # ($tmpstr =~ m{^/.*[^/]$}o || $tmpstr =~ m{^[^/].*/$}o) && do {
              #     $isOK_flg = $F;
              #     last REGEXPS;
              # }; # $tmpstr is a bad defined regexp
        $string = '^' . ( $escflg ? $tmpstr : quotemeta($tmpstr) ) . '$';    #'
    };    # REGEXPS
    $isOK_flg && do {
        unless ( eval { "" =~ m{$string}; $T } ) {
            $isOK_flg = $F;
        };    # check if final regexp string is OK
    };    # $isOK_flg
    $isOK_flg || do {
        &warn( 'BAD_REGEXP', $F, "$string" );
        $string = "";
    };    # NOT $isOK_flg
    return ( $isOK_flg, $not_flg, $id_flg, $string );
}    # find_regexp

sub eval_regexp() {
    my $str = $_[0];
    my $flag;
    eval { "" =~ m{$str}; $flag = $T; } || ( $flag = $F );
    return ( $str, $flag );
}    # eval_regexp

sub escape_input() {
    my $var = $_[0];
    my $c;
    $c = ( $var =~ s{([;,<>&!\{\}`'"])}{\\$1}g ) || 0;    #"'`
    return $var, ( $c > 0 ) ? $T : $F;
}    # escape_input

sub parse_special_customs() {
    my ( $string, $sflg ) = @_;
    my ( @psc, $tpsc );
    $string =~ s/^\s+//o;
    @psc = split /\s+/o, $string, 2;
    $tpsc = lc( $psc[0] );
    defined( $SpecialVars{$tpsc} ) || do {
        &warn( 'SPECIAL_NOT_DEF', $F, $tpsc );
        return $F;
    };
    return $SpecialVars{$tpsc}->( $psc[1], $sflg );
}    # parse_special_customs

sub get_spfeat_props() {
    my ( $cstr, $cflg, $gflg, $vna ) = @_;
    my ( @kao, @rr, $rr );
    my ( $a, $l, $s, $w, $f, $h, $v, $t, $fn, $fz, $fk, $str ) = (undef) x 12;
    $cflg && do {
        $cstr =~ s/\"(.+?)\"//o && ( $str = $1 );
    };    # $cflg
    @kao = split /\s+/og, $cstr;
    scalar(@kao) && do {
        foreach my $oc (@kao) {
            $oc =~ /^\s*$/o && next;     # empty fields
            $oc =~ /^A:(.*)$/o && do {
                $rr = $1;
                $TESTS{FLOAT}->( \$rr, $vna . '(angle)', $rr )
                  && ( $a = $rr );
                next;
            };
            $oc =~ /^S:(.*)$/o && do {
                $rr = $1;
                $TESTS{LINESTY}->( \$rr, $vna . '(line_style)', lc($rr) )
                  && ( $s = $rr );
                next;
            };
            $oc =~ /^W:(.*)$/o && do {
                $rr = $1;
                $TESTS{PS_UNIT}->( \$rr, $vna . '(line_width)', lc($rr) )
                  && ( $w = $rr );
                next;
            };
            if ($gflg) {
                $oc =~ /^L:(.*)$/o && do {
                    $rr = $1;
                    $TESTS{EXTCOLOR}->( \$rr, $vna . '(line_color)', lc($rr) )
                      && ( $l = $rr );
                    next;
                };
                $oc =~ /^F:(.*)$/o && do {
                    $rr = $1;
                    $TESTS{EXTCOLOR}->( \$rr, $vna . '(fill_color)', lc($rr) )
                      && ( $f = $rr );
                    next;
                };
                $oc =~ /^K:(.*)$/o && do {
                    $rr = $1;
                    $TESTS{CMYKFLOAT}->( \$rr, $vna . '(score)', $rr )
                      && ( $fk = $rr );
                    next;
                };
            }
            ;    # if $gflg
            if ($cflg) {
                $oc =~ /^H:(.*)$/o && do {
                    $rr = $1;
                    $TESTS{TXTALIGN}->( \$rr, $vna . '(horizontal_align)',
                      lc($rr), $T ) && ( $h = $rr );
                    next;
                };
                $oc =~ /^V:(.*)$/o && do {
                    $rr = $1;
                    $TESTS{TXTALIGN}->( \$rr, $vna . '(vertical_align)',
                      lc($rr), $F ) && ( $v = $rr );
                    next;
                };
                $oc =~ /^T:(.*)$/o && do {
                    $rr = $1;
                    $TESTS{COLOR}->( \$rr, $vna . '(text_color)', lc($rr) )
                      && ( $t = $rr );
                    next;
                };
                $oc =~ /^N:(.*)$/o && do {
                    $rr = $1;
                    $TESTS{FONT}->( \$rr, $vna . '(font_name)', lc($rr) )
                      && ( $fn = $rr );
                    next;
                };
                $oc =~ /^Z:(.*)$/o && do {
                    $rr = $1;
                    $TESTS{PS_UNIT}->( \$rr, $vna . '(font_size)', lc($rr) )
                      && ( $fz = &get_units( \@{$rr} ) );
                    next;
                };
            }
            ;    # if $cflg
        };    # foreach
    };    # scalar(@kao)
    $a = defined($a) ? $a : 0;
    @rr =
      ( '', "$a", defined($l) ? "$l" : 'fg',
        defined($s) ? "($LineStyles{$s})" : '(solid)',
        defined($w) ? "@{$w}"             : 'Clw',
        defined($f) ? "$f " . &tobool($T) : &tobool($F),
        defined($v) ? "($TextAlign{'V'}{$v})"
        : ( sin( $a * $torad ) >= 0 ? '(bv)' : '(tv)' ),
        defined($h) ? "($TextAlign{'H'}{$h})"
        : ( cos( $a * $torad ) >= 0 ? '(lh)' : '(rh)' ),
        defined($fz) ? "$fz" : '8 pt',
        defined($fn) ? $Fonts{$fn} : $Fonts{'helvetica'},
        defined($t) ? "$t " . &tobool($T) : &tobool($F),
        defined($fk) ? "$fk" : '1', );
    $rr[0] = ( defined($str) ) ? &tostring($str) : '';
    return @rr;
}    # get_spfeat_props

sub parse_GFF_files() {
    &header("PARSING INPUT GFF RECORDS");
    LOAD: foreach $file (@data_files) {
        open( THIS, "< $file" )
          || ( &warn( 'FILE_NO_OPEN', $T, $file ), next LOAD );
        $file eq '-' && ( $file = 'STANDARD INPUT' );
        &report( 'READ_GFF_FILE', $file );
        ( $n, $c ) = ( 0, undef );
        while (<THIS>) {
            my ( @line, $main );
            ( $c = '.', next ) if /^\#/o;
            ( $c = '.', next ) if /^\s*$/o;
            chomp;

            # $c = $noGFF;
            ( $main, undef ) = split /\s+\#/o;
            @line = split /\s+/o, $main, 9;
            scalar(@line) < 8
              && (
                &warn( 'NOT_ENOUGH_FIELDS', $F, $file, $n,
                    join ( " ", @line ) ), next );
            $c = &GFF_format( &fieldscheck( \@line ) );
        }
        continue {
            &counter( ++$n, $c );
        };    # WHILE
        &counter_end( $n, $c );
        close(THIS);
    };    # LOAD
    print LOGFILE ( Data::Dumper->Dump( [ \%GFF_DATA, \%ALN_DATA ],
          [qw( *GFF_DATA   *ALN_DATA )] ) )
      if ( $LogFile && $Debug );
    &footer("DATA LOADED");
}    # sub parse_GFF

sub GFF_format() {
    my $gff = $_[0];

    # return "x" if $GFF == $version1;
    return $GFF        if $gff eq $GFF;         # $version2
    return $GFF_NOGP   if $gff eq $GFF_NOGP;    # $version2 (ungrouped)
    return $VECTOR     if $gff eq $VECTOR;      # VECTOR: GFFv2 particular case
    return $ALIGN      if $gff eq $ALIGN;       # ALIGN: GFFv2 particular case
    return $APLOT      if $gff eq $APLOT;       # Old aplot format (with colons)
    return $APLOT_NOGP if $gff eq $APLOT_NOGP;  # Old aplot format (ungrouped)
    return $noGFF;
}    # GFF_format

sub fieldscheck() {
    my ($list) = @_;
    my ( $seqname, $start, $end ) = @$list[ 0, 3, 4 ];

    # ($list->[0],$list->[3],$list->[4]);
    ( &fcolon($seqname) && &fcolon($start) && &fcolon($end) ) && do {
        return &load_aplot($list);
    };
    return &load_gff($list);
}    # fieldscheck
sub fcolon() { return ( $_[0] =~ /^.+:.+$/o ? $T : $F ) }

sub check_coords() {    # ((ori,end)_1,...,(ori,end)_n)
    my ( $aflg, @ary ) = @_;
    for ( my $j = 0 ; $j <= $#ary ; $j += 2 ) {
        $ary[$j] > $ary[ $j + 1 ] && do {
            &warn( 'ORI_GREATER_END', $F, $ary[$j], $ary[ $j + 1 ], $file,
                $n + 1 );
            return $F;
        };    # $ary[$j] > $ary[$j+1]
    };    # for
    return $T;
}    # check_coords

sub check_score() {    # (scoref_1,...,scoref_n)
    foreach my $sco (@_) {
        ( defined($$sco) && $$sco eq '.' ) && ( $$sco = 0 );
    };    # foreach
}    # check_score

sub check_strand() {    # (str_1,...,srt_n)
    foreach my $str (@_) {
        $str !~ /$regexp_strand/o && do {
            &warn( 'STRAND_MISMATCH', $F, $str, $file, $n + 1 );
            return $F;
        };
    };    # foreach
    return $T;
}    # check_strand

sub check_frame() {    # (frm_1,...,frm_n)
    foreach my $frm (@_) {
        $frm !~ /$regexp_frame/o && do {
            &warn( 'FRAME_MISMATCH', $F, $frm, $file, $n + 1 );
            return $F;
        };
    };    # foreach
    return $T;
}    # check_frame

sub load_gff() {    # if errors found > return $noGFF
    my ($list) = @_;
    my $w_gff;
    ( $seqname, $source, $feature, $start, $end, $score, $strand, $frame ) =
      @$list[ 0, 1, 2, 3, 4, 5, 6, 7 ];
    $w_gff = &load_grouping( $T, $list->[8] );
    &check_gff_fields($w_gff) || ( $w_gff = $noGFF );
    return $w_gff;
}    # load_gff

sub check_gff_fields() {
    &check_coords( $T, $start, $end ) || ( return $F );
    &check_strand($strand) || ( return $F );
    &check_score( \$score );
    &check_frame($frame) || ( return $F );
    &add_gff_record( $_[0] );
    return $T;
}    # check_gff_fields

sub load_aplot() {    # if errors found > return $noGFF
    my ($list) = @_;
    my $w_gff;
    ( $seqname_1, $seqname_2, $source_1, $source_2, $feature_1,
      $feature_2, $start_1,   $start_2,  $end_1,    $end_2,
      $strand_1,  $strand_2,  $frame_1,  $frame_2 )
      = &remove_colon( @$list[ 0, 1, 2, 3, 4, 6, 7 ] );
    defined($strand_2) || ( $strand_2 = $strand_1 );
    defined($frame_2)  || ( $frame_2  = $frame_1 );
    ( $score_1, $score_2 ) = ( $list->[5], undef );
    $w_gff = &load_grouping( $F, $list->[8] );
    &check_aplot_fields($w_gff) || ( $w_gff = $noGFF );
    return $w_gff;
}    # load_aplot

sub remove_colon() {
    my @ary_out = ();
    my ( $a, $b ) = ( undef, undef );
    foreach my $fld (@_) {
        ( $a, $b ) = split /:/o, $fld, 2;
        push @ary_out,
          ( ( defined($a) ? $a : '.' ), ( defined($b) ? $b : undef ) );
    }
    return @ary_out;
}    # remove_colon

sub check_aplot_fields() {
    &check_coords( $F, $start_1, $end_1, $start_2, $end_2 ) || ( return $F );
    &check_strand( $strand_1, $strand_2 ) || ( return $F );
    &check_score( \$score_1, \$score_2 );
    &check_frame( $frame_1, $frame_2 ) || ( return $F );
    &add_gff_record( $_[0] );
    return $T;
}    # check_aplot_fields

sub add_gff_record() {
    my $_gff = $_[0];
    my ( $VarName, $Counter, $Type, $isaln, $myfunc, $t );
    $isaln = ( $_gff =~ /$ALIGN|$APLOT|$APLOT_NOGP/o ) ? $T : $F;
    ISALN: {
        $isaln && do {
            $myfunc = \&load_aln_var;
            ( $VarName, $Counter, $Type ) =
              ( \%ALN_DATA, \$aln_COUNT, 'SEQUENCE' );
            $seqname = join ( ":", $seqname_1, $seqname_2 );
            $source =
              defined($source_2) 
              ? join ( ":", $source_1, $source_2 )
              : $source_1;
            $strand = defined($strand_2) ? $strand_1 . $strand_2 : $strand_1;
            $feature =
              defined($feature_2) 
              ? join ( ":", $feature_1, $feature_2 )
              : $feature_1;
            last ISALN;
        };    # $isaln
        $myfunc = \&load_gff_var;
        ( $VarName, $Counter, $Type ) = ( \%GFF_DATA, \$seq_COUNT, 'SEQUENCE' );
    };    # ISALN
    &load_var( $seqname, $VarName, $Counter, $Type, $myfunc );
    ( $VarName, $Counter, $Type ) = (
      \%{ $VarName->{$seqname}[$_element] },
      \$$VarName{$seqname}[$_counter][$_elemNum],
      'SOURCE'
    );
    &load_var( $source, $VarName, $Counter, $Type, $myfunc );
    ( $VarName, $Counter, $Type ) = (
      \%{ $VarName->{$source}[$_element] },
      \$$VarName{$source}[$_counter][$_elemNum],
      'STRAND'
    );
    &load_var( $strand, $VarName, $Counter, $Type, $myfunc );
    ( $VarName, $Counter, $Type ) = (
      \%{ $VarName->{$strand}[$_element] },
      \$$VarName{$strand}[$_counter][$_elemNum],
      'GROUP'
    );
    &load_var( $group, $VarName, $Counter, $Type, $myfunc );
    $feature =
      lc($feature);    # GFF feature (3rd field) [MUST NOT be case-sensitive]
    LDALN: {
        $isaln && do {
            push @{ $VarName->{$group}[$_element] }, [
                'A',         # Type == ALIGNMENT
                {},          # Properties hash, now empty
                $feature,    # GFF feature [AGAIN: MUST NOT be case-sensitive]
                $label,      # Record ID if exist, order# otherwise
                $start_1, $end_1, $score_1, $frame_1,
                $start_2, $end_2, $score_2, $frame_2,
            ];
            last LDALN;
        };    # $isaln
        my $first_fld = ( $_gff eq $VECTOR ) ? 'V' : 'G';
        push @{ $VarName->{$group}[$_element] }, [
            $first_fld,    # Type == plain GFF or vector
            {},            # Properties hash, now empty
            $feature,      # GFF feature [AGAIN: MUST NOT be case-sensitive]
            $label,        # Record ID if exist, order# otherwise
            $start, $end, $score, $frame,
        ];
        $_gff eq $VECTOR && do {
            @{ $VarName->{$group}[$_element][8] } = [@vect_ary];
        };
    };    # LDALN
    $t = ++$VarName->{$group}[$_counter][$_elemNum];
    &set_var_defaults( 'FEATURE',
        \%{ $VarName->{$group}[$_element][ ( $t - 1 ) ][$_prop] } );
    return;
}    # add_gff_record

sub load_var() {
    my ( $_value, $_var, $_cnt, $_type, $_frf ) = @_;
    defined( $$_var{$_value} ) || do {
        $$_var{$_value}[$_counter] = &$_frf( ++$$_cnt );
        &set_var_defaults( $_type, \%{ $$_var{$_value}[$_prop] } );
    };
    return;
}    # load_var
sub load_gff_var() { return [ $_[0], (0) x 7 ]; }
sub load_aln_var() { return [ $_[0], (0) x 9 ]; }

sub set_var_defaults() {
    my ( $sect, $varhash ) = @_;

    # $$varhash = \%{$Defaults{$sect}};
    return;
}    # set_var_defaults

sub load_grouping() {
    my ( $_type, $attributes ) = @_;
    my ( $grp_string, $grp_counter,  $grp_GP, $grp_NOGP,  $grp_tag );
    my ( $grp_flag,   $group_string, @tt,     @new_group, @grouping_list );
    GFF_CHOICE: {
        $_type && do {
            $grp_string  = "$seqname\_$source\_$strand";
            $grp_counter = ++$group_gff_counter;
            $grp_GP      = $GFF;
            $grp_NOGP    = $GFF_NOGP;
            $grp_tag     = '';
            last GFF_CHOICE;
        };
        $grp_string  = "$seqname_1\_$seqname_2\_$source_1\_$strand_1$strand_2";
        $grp_counter = ++$group_aplot_counter;
        $grp_GP      = $APLOT;
        $grp_NOGP    = $APLOT_NOGP;
        $grp_tag = $Vars{LAYOUT}{align_tag};  # %SOURCE is a temporary hash name
    };
    $label    = '';
    $group_id = $grp_counter;
    defined($attributes) || do {
        $group = "$grp_string\_$group_id";
        return $grp_NOGP;
    };
    @grouping_list = split /\s*;\s*/og, $attributes;
    $grp_flag = 0;
    $group_string = shift @grouping_list;
    ( $group_string =~ /$regexp_group/o ) && do {
        @new_group =
          ( defined($1) ? $1 : '', defined($2) ? $2 : '',
            defined($3) ? $3 : '' );
    };
    $new_group[0] =~ s/\s*$//o;
    ( $new_group[0] eq '' ) && do {    # type 2 attributes
        $grp_flag = 1;
        $new_group[0] = $grp_tag;
    };
    ( $new_group[1] eq '' ) && do {    # type 1 attributes
        $grp_flag = 1;
        $new_group[1] = $new_group[0];
        $new_group[0] = $grp_tag;
    };
    ( $tag, $group ) = ( lc( $new_group[0] ), $new_group[1] );

    # Here looking for colon field separator in aplot GFF-like grouping
    ( $grp_flag && $group =~ /^(.*?):(.*?)$/o )
      && ( ( $group, $label ) = ( $1, $2 ) );
    $_type && do {    # GFF grouping
        @tt =
          ( ( ( defined( $new_group[2] ) && $new_group[2] ne '' )
            ? $new_group[2]
            : '@@NO@@NE@@' ), @grouping_list );
        $tag =~ /^$Vars{LAYOUT}{align_tag}$/ && do {
            return &load_GFF_align( \@tt );
        };
        $tag =~ /^$Vars{LAYOUT}{vector_tag}$/ && do {
            return &load_GFF_vector( \@tt );
        };
    };
    for ( my $element = 0 ; $element <= $#grouping_list ; $element++ ) {
        $grouping_list[$element] =~ /$regexp_group/o && do {
            @new_group =
              ( defined($1) ? $1 : '', defined($2) ? $2 : '',
                defined($3) ? $3 : '' );
        };
        $new_group[0] =~ s/\s*$//og;
        lc( $new_group[0] ) =~ /^$Vars{LAYOUT}{label_tag}$/ && do {
            ( $new_group[1] eq '' ) || ( $label = $new_group[1] );

            #            $label = $new_group[1];
            #            $label eq "" && do {
#                (undef,$label,undef) = split /\s+/og, $new_group[0];
            #            };
        };
    };    # for $element
    return $grp_GP;
}    # load_grouping

sub load_GFF_align() {
    my ($strf) = @_;
    my ( $aln_ori, $aln_end, $aln_sco, $aln_str, $aln_frm, $aln_lbl ) =
      (undef) x 6;
    my ( @lst, $fst );
    ($fst) = shift @$strf;
    $fst ne '@@NO@@NE@@' && do {
        @lst = split /\s+/og, $fst;
        ( ( defined( $lst[0] ) && $lst[0] =~ /$regexp_float/o )
          && ( defined( $lst[1] ) && $lst[1] =~ /$regexp_float/o ) )
          && do {
            ( $aln_ori, $aln_end ) = @lst[ 0, 1 ];
            ( defined( $lst[2] ) && $lst[2] =~ /$regexp_strand/o ) && do {
                $aln_str = $lst[2];
                ( defined( $lst[3] ) && $lst[3] =~ /$regexp_frame/o ) && do {
                    $aln_frm = $lst[3];
                };    # $aln_frm
            };    # $aln_str
        };    # numbers?($lst[0] && $lst[1])
    };    # defined($strf[0])
    for ( my $gh = 0 ; $gh <= $#{$strf} ; $gh++ ) {
        my ( @q, $t, @new_id );
        @q = split /\s+/og, $strf->[$gh];
        $t = lc( $q[0] );

        # print STDERR join(' : ',$strf->[$gh], join('***',@q).'**|', 
        #                   $t, $Vars{LAYOUT}{label_tag})."\n";
        $t =~ /^($Vars{LAYOUT}{label_tag})$/ && do {
            ( $strf->[$gh] =~ /$regexp_group/o ) && do {
                @new_id =
                  ( defined($1) ? $1 : '', defined($2) ? $2 : '',
                    defined($3) ? $3 : '' );
            };
            ( $new_id[1] eq '' ) || ( $aln_lbl = $new_id[1] );
            next;
        };    # id_tag
        $t =~ /^$Vars{LAYOUT}{start_tag}$/ && do {
            ( defined( $q[1] ) && $q[1] =~ /$regexp_float/o )
              && ( $aln_ori = $q[1] );
            next;
        };    # start_tag
        $t =~ /^$Vars{LAYOUT}{end_tag}$/ && do {
            ( defined( $q[1] ) && $q[1] =~ /$regexp_float/o )
              && ( $aln_end = $q[1] );
            next;
        };    # end_tag
        $t =~ /^$Vars{LAYOUT}{strand_tag}$/ && do {
            ( defined( $q[1] ) && $q[1] =~ /$regexp_strand/o )
              && ( $aln_str = $q[1] );
            next;
        };    # strand_tag
        $t =~ /^$Vars{LAYOUT}{frame_tag}$/ && do {
            ( defined( $q[1] ) && $q[1] =~ /$regexp_frame/o )
              && ( $aln_frm = $q[1] );
        };    # frame_tag
        $t =~ /^$Vars{LAYOUT}{score_tag}$/ && do {
            ( defined( $q[1] ) && $q[1] =~ /$regexp_real/o )
              && ( $aln_sco = $q[1] );
        };    # frame_tag
    };    # for $gh
    ( defined($aln_ori) && defined($aln_end) ) || return $GFF;
    ( $seqname_1, $seqname_2 ) = ( $seqname, $group );
    ( $source_1,  $feature_1 ) = ( $source,  $feature );
    ( $source_2,  $feature_2 ) = (undef) x 2;
    ( $start_1,   $end_1 )     = ( $start,   $end );
    ( $start_2,   $end_2 )     = ( $aln_ori, $aln_end );
    ( $score_1,   $score_2 )   = ( $score,   $aln_sco );
    ( $strand_1,  $frame_1 )   = ( $strand,  $frame );
    $strand_2 = defined($aln_str) ? $aln_str : $strand_1;
    $frame_2  = defined($aln_frm) ? $aln_frm : $frame_1;
    $label    = defined($aln_lbl) ? $aln_lbl : $label;
    return $ALIGN;
}    # load_GFF_align

sub load_GFF_vector() {
    my ($strf) = @_;
    my ( $v_type, $v_window, $v_step, $vct_lbl ) =
      ( $vect_type{'sco'}, (undef) x 3 );
    my ( @lst, $fst );
    @vect_ary = ();
    ($fst) = shift @$strf;
    $fst ne '@@NO@@NE@@' && do {
        @lst = split /\s+/og, $fst;
        ( defined( $lst[0] ) && defined( $vect_type{ lc( $lst[0] ) } ) ) && do {
            $v_type = $vect_type{ lc( $lst[0] ) };
            ( defined( $lst[1] ) && $lst[1] =~ /$regexp_float/o ) && do {
                $v_window = $lst[1];
                ( defined( $lst[2] ) && $lst[2] =~ /$regexp_float/o ) && do {
                    $v_step = $lst[2];
                };    # $v_step
            };    # $v_window
        };    # $v_type
    };    # defined($strf[0])
    for ( my $gh = 0 ; $gh <= $#{$strf} ; $gh++ ) {
        my ( @q, $t, @new_id );
        @q = split /\s+/og, $strf->[$gh];
        $t = lc( $q[0] );

        # print STDERR join(' : ',$strf->[$gh], join('***',@q).'**|', 
        #                   $t, $Vars{LAYOUT}{label_tag})."\n";
        $t =~ /^($Vars{LAYOUT}{label_tag})$/ && do {
            ( $strf->[$gh] =~ /$regexp_group/o ) && do {
                @new_id =
                  ( defined($1) ? $1 : '', defined($2) ? $2 : '',
                    defined($3) ? $3 : '' );
            };
            ( $new_id[1] eq '' ) || ( $vct_lbl = $new_id[1] );
            next;
        };    # id_tag
        $t =~ /^($Vars{LAYOUT}{scotype_tag})$/ && do {
            defined( $vect_type{ lc( $strf->[$gh] ) } )
              && ( $v_type = $vect_type{ lc( $strf->[$gh] ) } );
            next;
        };    # scotype_tag
        $t =~ /^($Vars{LAYOUT}{window_tag})$/ && do {
            ( defined( $q[1] ) && $q[1] =~ /$regexp_float/o )
              && ( $v_window = $q[1] );
            next;
        };    # window_tag
        $t =~ /^($Vars{LAYOUT}{step_tag})$/ && do {
            ( defined( $q[1] ) && $q[1] =~ /$regexp_float/o )
              && ( $v_step = $q[1] );
            next;
        };    # step_tag
        $t =~ /^$Vars{LAYOUT}{scoary_tag}$/ && do {
            @vect_ary = split /\s+/og, $q[1];
        };    # frame_tag
    };    # for $gh
    scalar(@vect_ary) > 0 || return $GFF;
    ( $v_type == $vect_type{'sco'} && !defined($v_window) )
      && ( $v_type = $vect_type{'ERROR'} );
    defined($v_step) || ( $v_step = $v_window );
    @vect_ary = ( $v_type, $v_window, $v_step, @vect_ary );
    $label = defined($vct_lbl) ? $vct_lbl : $label;
    return $VECTOR;
}    # load_GFF_vector

sub sort_elements() {
    &header("SORTING ELEMENTS BY ACCEPTOR (START)");
    %Order = ();

    # sorting %GFF_DATA contents 
    scalar(%GFF_DATA) && do {
        &report( 'SORT_GFF', '*- ', 'ANNOTATION DATA' );
        &sort_elements_loop( \%GFF_DATA, 'GFF' );
    };    # scalar(%GFF_DATA) > 0
          # sorting %ALN_DATA contents
    scalar(%ALN_DATA) && do {
        &report( 'SORT_GFF', '*- ', 'ALIGNMENT DATA' );
        &sort_elements_loop( \%ALN_DATA, 'ALN' );
    };    # scalar(%ALN_DATA) > 0
    print LOGFILE (
      Data::Dumper->Dump(
          [ \%Order, \%GFF_DATA, \%ALN_DATA ],
          [qw( *Order   *GFF_DATA   *ALN_DATA )]
      ) )
      if ( $LogFile && $Debug );
    &footer("ELEMENTS SORTED");
}    # sort_elements

sub sort_elements_loop() {
    my ( $s_ref, $ktr ) = @_;
    my (
        $v_max, $v_min, $w_max, $w_min, $s_max,
        $s_min, $l_max, $L_max, $aflg
    );
    my $sq_ord = \@{ $Order{$ktr} };

    # my @tt_coords = ();
    $aflg = $ktr eq 'ALN' ? $T : $F;
    @{$sq_ord} = ();
    foreach my $s_seq ( keys %{$s_ref} ) {
        &report( 'SORT_SEQ', '|  *- ', $s_seq );
        my @sc_coords = ();
        push @{$sq_ord}, [ $s_seq, $s_ref->{$s_seq}[$_counter][$_order], () ];
        my $sc_ord = \@{ $sq_ord->[ $#{$sq_ord} ][2] };
        my $ss_ref = \%{ $s_ref->{$s_seq}[$_element] };
        foreach my $s_src ( keys %{$ss_ref} ) {
            &report( 'SORT_SRC', ( ( '|  ' x 2 ) . '*- ' ), $s_src );
            my @sr_coords = ();
            push @{$sc_ord},
              [ $s_src, $ss_ref->{$s_src}[$_counter][$_order], () ];
            my $sr_ord = \@{ $sc_ord->[ $#{$sc_ord} ][2] };
            my $sc_ref = \%{ $ss_ref->{$s_src}[$_element] };
            foreach my $s_str ( keys %{$sc_ref} ) {
                &report( 'SORT_STR', ( ( '|  ' x 3 ) . '*- ' ), $s_str );
                my @ft_coords = ();
                push @{$sr_ord}, [
                    $s_str,

                    # $sc_ref->{$s_str}[$_counter][$_order],
                    ()
                ];    # if uncomment '$sc_ref' set next to [2] instead of [1].
                my $st_ord    = \@{ $sr_ord->[ $#{$sr_ord} ][1] };
                my $sr_ref    = \%{ $sc_ref->{$s_str}[$_element] };
                my $sortfunct =
                  ( $s_str eq '-' ) ? \&sort_reverse : \&sort_forward;
                my ( @layer_ary, $s_elem );
                @layer_ary = (
                    \%{ $sc_ref->{$s_str}[$_prop] },
                    \%{ $ss_ref->{$s_src}[$_prop] },
                    \%{ $s_ref->{$s_seq}[$_prop] }
                );
                foreach my $s_grp ( keys %{$sr_ref} ) {
                    &report( 'SORT_GRP', ( ( '|  ' x 4 ) . '*- ' ), $s_grp );
                    my $sg_ref = \@{ $sr_ref->{$s_grp}[$_element] };
                    $s_elem = $sr_ref->{$s_grp}[$_counter][$_elemNum];
                    $s_elem > 1 && do {
                        @{$sg_ref} = map { $_->[3] } sort { &$sortfunct } map {
                            [
                                &get_layer(
                                    'feature_layer',    # {feature_layer}
                                    \%{ $_->[$_ftprop] },
                                    \%{ $sr_ref->{$s_grp}[$_prop] },
                                    @layer_ary
                                ),
                                $_->[$_ftori],
                                $_->[$_ftend],
                                $_
                            ];
                        } @{$sg_ref};    # maps layer,start,end,arrayelement
                    };    # $s_elem > 1
                     # @ft_coords = ( map { $_->[$_ftori], $_->[$_ftend] } @{ $sg_ref } );
                    $v_min = &min( map { $_->[$_ftori] } @{$sg_ref} );
                    $v_max = &max( map { $_->[$_ftend] } @{$sg_ref} );
                    $sr_ref->{$s_grp}[$_counter][$_ori] = $v_min;
                    $sr_ref->{$s_grp}[$_counter][$_end] = $v_max;
                    @ft_coords = ( map { $_->[$_ftsco] } @{$sg_ref} );
                    $s_min = &min(@ft_coords);
                    $s_max = &max(@ft_coords);
                    $sr_ref->{$s_grp}[$_counter][$_mnsco] = $s_min;
                    $sr_ref->{$s_grp}[$_counter][$_mxsco] = $s_max;

# @ft_coords = ( map { &max(length($_->[$_ftname]),length($_->[$_ftid])) }
                    #                   @{ $sg_ref } );
                    $l_max = &max(
                        map {
                            &max(
                                length( $_->[$_ftname] ),
                                length( $_->[$_ftid] )
                            );
                        } @{$sg_ref}
                    );
                    $sr_ref->{$s_grp}[$_counter][$_flw] = $l_max;
                    $L_max = length($s_grp);

                    #
                    $aflg && do {

# @ft_coords = ( map { $_->[$_ftnori], $_->[$_ftnend] } @{ $sg_ref } );
                        $w_min = &min( map { $_->[$_ftnori] } @{$sg_ref} );
                        $w_max = &max( map { $_->[$_ftnend] } @{$sg_ref} );
                        $sr_ref->{$s_grp}[$_counter][$_nori] = $w_min;
                        $sr_ref->{$s_grp}[$_counter][$_nend] = $w_max;
                        push @{$st_ord},
                          [
                            $s_grp, $v_min, $v_max, $s_min, $s_max,
                            $l_max, $L_max, $w_min, $w_max
                        ];
                    };
                    $aflg || do {
                        push @{$st_ord},
                          [
                            $s_grp, $v_min, $v_max, $s_min,
                            $s_max, $l_max, $L_max
                        ];
                    };

                    #
                    &report( 'SORT_FTR', ( '|  ' x 5 ), $s_elem );
                };    # foreach $s_grp
                shift @layer_ary;

                # @ft_coords = ( map { $_->[1], $_->[2] } @{ $st_ord } );
                $v_min = &min( map { $_->[1] } @{$st_ord} );
                $v_max = &max( map { $_->[2] } @{$st_ord} );
                $sc_ref->{$s_str}[$_counter][$_ori] = $v_min;
                $sc_ref->{$s_str}[$_counter][$_end] = $v_max;

                # @ft_coords = ( map { $_->[3], $_->[4] } @{ $st_ord } );
                $s_min = &min( map { $_->[3] } @{$st_ord} );
                $s_max = &max( map { $_->[4] } @{$st_ord} );
                $sc_ref->{$s_str}[$_counter][$_mnsco] = $s_min;
                $sc_ref->{$s_str}[$_counter][$_mxsco] = $s_max;

                # @ft_coords = ( map { $_->[5] } @{ $st_ord } );
                $l_max = &max( map { $_->[5] } @{$st_ord} );
                $L_max = &max( map { $_->[6] } @{$st_ord} );
                $sc_ref->{$s_str}[$_counter][$_flw] = $l_max;
                $sc_ref->{$s_str}[$_counter][$_glw] = $L_max;
                $aflg && do {

                    # @ft_coords = ( map { $_->[6], $_->[7] } @{ $st_ord } );
                    $w_min = &min( map { $_->[7] } @{$st_ord} );
                    $w_max = &max( map { $_->[8] } @{$st_ord} );
                    $sc_ref->{$s_str}[$_counter][$_nori] = $w_min;
                    $sc_ref->{$s_str}[$_counter][$_nend] = $w_max;
                    push @sr_coords,
                      [
                        $v_min, $v_max, $s_min, $s_max,
                        $l_max, $L_max, $w_min, $w_max
                    ];
                };

                # push @sr_coords, $v_min, $v_max;
                $aflg || do {
                    push @sr_coords,
                      [ $v_min, $v_max, $s_min, $s_max, $l_max, $L_max ];
                };

                #
                @{$st_ord} = map { $_->[3] } sort { &$sortfunct } map {
                    [
                        &get_layer(
                            'group_layer',    # {group_layer}
                            \%{ $sr_ref->{ $_->[0] }[$_prop] },
                            @layer_ary
                        ),
                        $_->[1],
                        $_->[2],
                        $_->[0]
                    ];
                } @{$st_ord};
                $s_elem = $sc_ref->{$s_str}[$_counter][$_elemNum];
                &report( 'SORT_GPN', ( '|  ' x 4 ), $s_elem );
            };    # foreach $s_str
            @{$sr_ord} = map { $_->[2] } sort {
                $a->[0] <=> $b->[0]         # sort by layer
                  || $a->[1] <=> $b->[1]    # sort by "strand number"
              } map {
                [
                    &get_layer(
                        'strand_layer',    # {strand_layer}
                        \%{ $sc_ref->{ $_->[0] }[$_prop] },
                        \%{ $ss_ref->{$s_src}[$_prop] },
                        \%{ $s_ref->{$s_seq}[$_prop] }
                    ),
                    $StrandNum{ $_->[0] },
                    $_
                ];
            } @{$sr_ord};

            # @tt_coords = ( map { $_->[0], $_->[1] } @sr_coords );
            $v_min = &min( map { $_->[0] } @sr_coords );
            $v_max = &max( map { $_->[1] } @sr_coords );
            $ss_ref->{$s_src}[$_counter][$_ori] = $v_min;
            $ss_ref->{$s_src}[$_counter][$_end] = $v_max;

            # @tt_coords = ( map { $_->[2], $_->[3] } @sr_coords );
            $s_min = &min( map { $_->[2] } @sr_coords );
            $s_max = &max( map { $_->[3] } @sr_coords );
            $ss_ref->{$s_src}[$_counter][$_mnsco] = $s_min;
            $ss_ref->{$s_src}[$_counter][$_mxsco] = $s_max;
            $l_max = &max( map { $_->[4] } @sr_coords );
            $L_max = &max( map { $_->[5] } @sr_coords );
            $ss_ref->{$s_src}[$_counter][$_flw] = $l_max;
            $ss_ref->{$s_src}[$_counter][$_glw] = $L_max;
            $aflg && do {

                # @tt_coords = ( map { $_->[5], $_->[6] } @sr_coords );
                $w_min = &min( map { $_->[6] } @sr_coords );
                $w_max = &max( map { $_->[7] } @sr_coords );
                $ss_ref->{$s_src}[$_counter][$_nori] = $w_min;
                $ss_ref->{$s_src}[$_counter][$_nend] = $w_max;
                push @sc_coords,
                  [
                    $v_min, $v_max, $s_min, $s_max,
                    $l_max, $L_max, $w_min, $w_max
                ];
            };

            # push @sc_coords, $v_min, $v_max;
            $aflg || do {
                push @sc_coords,
                  [ $v_min, $v_max, $s_min, $s_max, $l_max, $L_max ];
            };
        };    # foreach $s_src
        &sort_by_inputorderlayer( $sc_ord, $ss_ref,
            \%{ $s_ref->{$s_seq}[$_prop] } );

        # @tt_coords = ( map { $_->[0], $_->[1] } @sc_coords );
        $s_ref->{$s_seq}[$_counter][$_ori] = &min( map { $_->[0] } @sc_coords );
        $s_ref->{$s_seq}[$_counter][$_end] = &max( map { $_->[1] } @sc_coords );

        # @tt_coords = ( map { $_->[2], $_->[3] } @sc_coords );
        $s_ref->{$s_seq}[$_counter][$_mnsco] =
          &min( map { $_->[2] } @sc_coords );
        $s_ref->{$s_seq}[$_counter][$_mxsco] =
          &max( map { $_->[3] } @sc_coords );
        $s_ref->{$s_seq}[$_counter][$_flw] = &max( map { $_->[4] } @sc_coords );
        $s_ref->{$s_seq}[$_counter][$_glw] = &max( map { $_->[5] } @sc_coords );
        $aflg && do {

            # @tt_coords = ( map { $_->[5], $_->[6] } @sc_coords );
            $s_ref->{$s_seq}[$_counter][$_nori] =
              &min( map { $_->[6] } @sc_coords );
            $s_ref->{$s_seq}[$_counter][$_nend] =
              &max( map { $_->[7] } @sc_coords );
        };
    };    # foreach $s_seq
    &sort_by_inputorder($sq_ord);
    return;
}    #sort_elements_loop

sub get_layer() {
    my ( $name, @refs ) = @_;
    my ( $out, $v );
    $out = 0;    #
                 # print STDERR ("#"x40)."\n";
    foreach $v (qw/ FEATURE GROUP STRAND SOURCE SEQUENCE /) {

        # print STDERR "### GET LAYER: $v ($name) : ".
#              (defined($Defaults{$v}{$name})?$Defaults{$v}{$name}:"undef")."\n";
        defined( $Defaults{$v}{$name} ) && ( $out = $Defaults{$v}{$name} );
    };    # foreach
    foreach $v (@refs) {

        # print STDERR "### GET LAYER: $v ($name) : ".
        #              (defined($v->{$name})?$v->{$name}:"undef")."\n";
        defined( $v->{$name} ) && ( $out = $v->{$name} );
    };    # foreach
          # print STDERR "### GET LAYER: $name = $out \n";
    return $out;
}    # get_layer

sub sort_by_inputorder() {
    my $ref = $_[0];
    @{$ref} =
      map { [ $_->[0], $_->[2] ] }
      sort { $a->[1] <=> $b->[1] }
      map { [ $_->[0], $_->[1], $_->[2] ] } @{$ref};
}    # sort_by_inputorder

sub sort_by_inputorderlayer() {
    my ( $ref, $soref, $sqref ) = @_;
    @{$ref} = map { [ $_->[0], $_->[3] ] } sort {
        $a->[1] <=> $b->[1]         # {source_layer}
          || $a->[2] <=> $b->[2]    # input order
      } map {
        [
            $_->[0],
            &get_layer(
                'source_layer', \%{ $soref->{ $_->[0] }[$_prop] },
                $sqref
            ),
            $_->[1],
            $_->[2]
        ];
    } @{$ref};
}    # sort_by_inputorderlayer

sub sort_forward {
    $a->[0] <=> $b->[0]          # sorting by layer
      or $a->[1] <=> $b->[1]     # sorting by start
      or $b->[2] <=> $a->[2];    # reverse sorting by end if same start
}    # sort_forward
     #

sub sort_reverse {
    $a->[0] <=> $b->[0]          # sorting by layer
      or $b->[2] <=> $a->[2]     # reverse sorting by end
      or $a->[1] <=> $b->[1];    # sorting by start if same end
}    # sort_reverse

sub merge_custom_vars() {
    &header("MAPPING CUSTOMIZATION INPUTS TO MAIN VARS");
    foreach my $var_name ( keys %{ $Defaults{LAYOUT} } ) {
        defined( $CmdLineVars{LAYOUT}{$var_name} ) && do {
            $Vars{LAYOUT}{$var_name} = $CmdLineVars{LAYOUT}{$var_name};
            next;
        };
        defined( $CustomVars{LAYOUT}{$var_name} ) && do {
            $Vars{LAYOUT}{$var_name} = $CustomVars{LAYOUT}{$var_name};
            next;
        };
        $Vars{LAYOUT}{$var_name} = $Defaults{LAYOUT}{$var_name};
    };    # foreach $var_name
    foreach my $_sec ( keys %Defaults ) {
        $_sec eq 'LAYOUT' && next;    # skip layout variables
        foreach my $_var ( keys %{ $Defaults{$_sec} } ) {
            defined( $CustomVars{$_sec}{$_var} ) && do {
                push @{ $Vars{$_sec}{$_var} }, @{ $CustomVars{$_sec}{$_var} };
            };    # defined($CustomVars{$_sec}{$_var})
            defined( $CmdLineVars{$_sec}{$_var} ) && do {
                push @{ $Vars{$_sec}{$_var} }, @{ $CmdLineVars{$_sec}{$_var} };
            };    # defined($CustomVars{$_sec}{$_var})
        };    # foreach $vnm
    };    # foreach $sct
    print LOGFILE ( Data::Dumper->Dump( [ \%Vars ], [qw( *Vars )] ) )
      if ( $LogFile && $Debug );
    &footer("VALUES SET for MAIN VARS");
}    # merge_custom_vars

sub map_vars_data() {
    &header("SETTING CUSTOM VALUES TO GFF ELEMENTS");
    foreach my $v_sec ( keys %Vars ) {
        $v_sec eq 'LAYOUT' && next;    # skip layout variables
        foreach my $v_var ( keys %{ $Vars{$v_sec} } ) {
            my @v_values = @{ $Vars{$v_sec}{$v_var} };
            for ( my $foo = 0 ; $foo < $#v_values ; $foo += 4 ) {
                my @tl = @v_values[ $foo .. ( $foo + 3 ) ];
                &map_vars_to_GFF( \%GFF_DATA, $v_sec, $v_var, @tl );
                &map_vars_to_GFF( \%ALN_DATA, $v_sec, $v_var, @tl );
            };    # for $foo
        };    # foreach $vnm
    };    # foreach $sct
    print LOGFILE ( Data::Dumper->Dump( [ \%GFF_DATA, \%ALN_DATA ],
          [qw( *GFF_DATA   *ALN_DATA )] ) )
      if ( $LogFile && $Debug );
    &footer("VALUES SET for GFF ELEMENTS");
}    # map_vars_data

sub map_vars_to_GFF() {
    my (
        $mainref, $v_sec,   $v_var,     $neg_flg,
        $id_flg,  $reg_exp, $the_value, $name_test
      )
      = @_;
    $name_test = $neg_flg ? \&match_regexp_neg : \&match_regexp;
    foreach my $seq ( keys %{$mainref} ) {
        my $seq_ref = \@{ $mainref->{$seq} };
        $v_sec eq 'SEQUENCE' && do {
            &does_feat_match(
                $name_test, $seq,   $reg_exp, $seq_ref,
                $_prop,     $v_sec, $v_var,   $the_value
            );
            next;
        };
        foreach my $src ( keys %{ $seq_ref->[$_element] } ) {
            my $src_ref = \@{ $seq_ref->[$_element]{$src} };
            $v_sec eq 'SOURCE' && do {
                &does_feat_match(
                    $name_test, $src,   $reg_exp, $src_ref,
                    $_prop,     $v_sec, $v_var,   $the_value
                );
                next;
            };
            foreach my $str ( keys %{ $src_ref->[$_element] } ) {
                my $str_ref = \@{ $src_ref->[$_element]{$str} };
                $v_sec eq 'STRAND' && do {
                    &does_feat_match(
                        $name_test, $str,   $reg_exp, $str_ref,
                        $_prop,     $v_sec, $v_var,   $the_value
                    );
                    next;
                };
                foreach my $grp ( keys %{ $str_ref->[$_element] } ) {
                    my $grp_ref = \@{ $str_ref->[$_element]{$grp} };
                    $v_sec eq 'GROUP' && do {
                        &does_feat_match(
                            $name_test, $grp,   $reg_exp, $grp_ref,
                            $_prop,     $v_sec, $v_var,   $the_value
                        );
                        next;
                    };
                    foreach my $feat ( 0 .. $#{ $grp_ref->[$_element] } ) {
                        my ( $feat_ref, $ft_name, $ft_id, $ft_vars );
                        $feat_ref = \@{ $grp_ref->[$_element][$feat] };
                        ( $ft_name, $ft_id ) =
                          ( $feat_ref->[$_ftname], $feat_ref->[$_ftid] );

                        # check ID
                        defined($id_flg) && do {
                            $ft_id =~ /^$id_flg$/o && do {
                                &does_feat_match(
                                    $name_test, $ft_name, $reg_exp,
                                    $feat_ref,  $_ftprop, $v_sec,
                                    $v_var,     $the_value
                                );
                            };    # $ft_id =~ /^$id_flg$/o
                            next;
                        };    # $id_flg ne $NULL
                              # check regexp
                        &does_feat_match(
                            $name_test, $ft_name, $reg_exp, $feat_ref,
                            $_ftprop,   $v_sec,   $v_var,   $the_value
                        );
                    };    # foreach my $feat
                };    # foreach my $grp
            };    # foreach my $str
        };    # foreach my $src
    };    # foreach my $seq
    return;
}    # map_vars_to_GFF

sub does_feat_match() {
    my ( $thetest, $name, $rexp, $gffref, $prop, $sct, $var, $value ) = @_;
    &$thetest( $name, $rexp ) && do {
        ( ref( $gffref->[$prop] ) eq 'REF' )
          && &set_all_defaults( $gffref, $sct );
        $gffref->[$prop]{$var} = $value;
    };
}    # does_feat_match

sub match_regexp() {
    $_[0] =~ /$_[1]/ && return $T;
    return $F;
}    # match_regexp

sub match_regexp_neg() {
    $_[0] !~ /$_[1]/ && return $T;
    return $F;
}    # match_regexp_neg

sub set_all_defaults() {
    my ( $the_hash, $the_sec ) = @_;
    $the_hash->[$_prop] = ();
    foreach my $vnm ( keys %{ $Defaults{$the_sec} } ) {
        $the_hash->[$_prop]{$vnm} = \$Defaults{$the_sec}{$vnm};
    };    # foreach $nm
}    # set_all_defaults

sub set_page_vars() {
    &header("SETTING PAGE LAYOUT");
    my $var = \%{ $Vars{LAYOUT} };
    &set_page_size($var);
    &set_page_axes($var);
    &set_seq_boundaries($var);
    &set_score_boundaries($var);
    &set_tickmark_vars($var);
    &find_prop_ratio($var);
    &set_misc_vars($var);
    &set_page_labels($var);
    &footer("PAGE LAYOUT SET for CURRENT PLOT");
}    # set_page_vars

sub set_page_size() {
    my ($vrf) = @_;
    PSIZES: {
        defined( $vrf->{page_bbox} ) && do {
            $vrf->{_page_width}  = &get_units( \@{ $vrf->{page_bbox}[1] } );
            $vrf->{_page_height} = &get_units( \@{ $vrf->{page_bbox}[2] } );
            $FORMATS{ $vrf->{page_bbox}[0] } =
              [ ++$formats, $vrf->{_page_width}, $vrf->{_page_height} ];
            $vrf->{page_size} = $vrf->{page_bbox}[0];
            last PSIZES;
        };
        ( $vrf->{_page_width}, $vrf->{_page_height} ) =
          @{ $FORMATS{ $vrf->{page_size} } }[ 1, 2 ];
    };    # PSIZES
    &report(
        'SETPAGESIZE', $vrf->{page_size},
        "$vrf->{_page_width} x $vrf->{_page_height}"
    );
    return;
}    # set_page_size

sub set_page_margins() {
    my ($vrf) = @_;
    $vrf->{_page_margins} = [
        &get_units( \@{ $vrf->{margin_left}[0] } ),
        &get_units( \@{ $vrf->{margin_right}[0] } ),
        &get_units( \@{ $vrf->{margin_top}[0] } ),
        &get_units( \@{ $vrf->{margin_bottom}[0] } )
    ];
    $vrf->{_page_clip_width} =
      $vrf->{_page_width} -
      ( $vrf->{_page_margins}[0] + $vrf->{_page_margins}[1] );
    $vrf->{_page_clip_height} =
      $vrf->{_page_height} -
      ( $vrf->{_page_margins}[2] + $vrf->{_page_margins}[3] );
    return;
}    # set_page_margins

sub set_page_axes() {
    my ($vrf) = @_;
    $vrf->{_x_sequence_name} =
      defined( $vrf->{x_sequence_name} ) ? $vrf->{x_sequence_name} : '.+';
    $vrf->{_y_sequence_name} =
      defined( $vrf->{y_sequence_name} ) ? $vrf->{y_sequence_name} : '.+';
    if ( defined( $vrf->{alignment_name} ) ) {
        ( $vrf->{_x_sequence_name}, $vrf->{_y_sequence_name} ) = split /:/o,
          $vrf->{alignment_name}, 2;
    }
    elsif ( scalar( @{ $Order{ALN} } ) > 0 ) {
        ( $vrf->{_x_sequence_name}, $vrf->{_y_sequence_name} ) = split /:/o,
          $Order{ALN}[0][0], 2;
    }
    ;    # defined($vrf->{alignment_name})
    $vrf->{_alignment_name} =
      [ $vrf->{_x_sequence_name}, $vrf->{_y_sequence_name} ];

    #
# $Debug && print STDERR "###OOPS(1)### $vrf->{_x_sequence_name} x $vrf->{_y_sequence_name} = @{$vrf->{_alignment_name}}\n";
    #
    ( $vrf->{_plot_aln}, $vrf->{_swap_aln} ) =
      &search_alnname( \@{ $vrf->{_alignment_name} }, \@{ $Order{ALN} } );
    $vrf->{_swap_aln} && do {
        ( $vrf->{x_sequence_name}, $vrf->{y_sequence_name} ) =
          ( $vrf->{y_sequence_name}, $vrf->{x_sequence_name} );
    };    # $vrf->{_swap_aln}
    $vrf->{_plot_aln} && do {
        $vrf->{_x_sequence_name} =
          defined( $vrf->{x_sequence_name} ) 
          ? $vrf->{x_sequence_name}
          : $vrf->{_alignment_name}[0];
        $vrf->{_y_sequence_name} =
          defined( $vrf->{y_sequence_name} ) 
          ? $vrf->{y_sequence_name}
          : $vrf->{_alignment_name}[1];
    };    # $vrf->{_plot_aln}
          #
     # $Debug && print STDERR "###OOPS(2)### $vrf->{_x_sequence_name} x $vrf->{_y_sequence_name} = @{$vrf->{_alignment_name}}\n";
     #
    ( $vrf->{_plot_x_sequence}, $vrf->{_x_sequence_name} ) =
      &search_seqname( $vrf->{_x_sequence_name}, undef, \@{ $Order{GFF} } );
    ( $vrf->{_plot_y_sequence}, $vrf->{_y_sequence_name} ) = &search_seqname(
        $vrf->{_y_sequence_name}, $vrf->{_x_sequence_name},
        \@{ $Order{GFF} }
    );
    ( $vrf->{_report_x_sequence_name}, $vrf->{_report_y_sequence_name} ) =
      ( $vrf->{_x_sequence_name}, $vrf->{_y_sequence_name} );
    $vrf->{_alignment_name} =
      join ( ":",
        $vrf->{_alignment_name}[0] ne '.+' ? $vrf->{_alignment_name}[0]
        : $vrf->{_x_sequence_name},
        $vrf->{_alignment_name}[1] ne '.+' ? $vrf->{_alignment_name}[1]
        : $vrf->{_y_sequence_name} );
    $vrf->{_report_alignment_name} = $vrf->{_alignment_name};

    #
# $Debug && print STDERR "###OOPS(3)### $vrf->{_x_sequence_name} x $vrf->{_y_sequence_name} = $vrf->{_alignment_name}\n";
    #
    # ($vrf->{alignment_name} = $vrf->{_alignment_name}) =~ s/$;/ : /o;
    my ( $yf, $nf ) =
      ( 'was found. PLOTTING IT...', 'NOT found. NOT PLOTTED...' );
    &report( 'SETPAGEAXES',
        ( "ALIGNMENT \"$vrf->{_report_alignment_name}\" "
        . ( $vrf->{_plot_aln} ? $yf : $nf )
        . ( $vrf->{_swap_aln}
        ? "\n###     The alignment found makes program"
        . " to swap X-Y annotations..."
        : "" ) ),
        ( "X-SEQUENCE \"$vrf->{_report_x_sequence_name}\" "
        . ( $vrf->{_plot_x_sequence} ? $yf : $nf ) ),
        ( "Y-SEQUENCE \"$vrf->{_report_y_sequence_name}\" "
        . ( $vrf->{_plot_y_sequence} ? $yf : $nf ) ) );
    return;
}    # set_page_axes

sub search_alnname() {
    my ( $_ra, $_ro ) = @_;
    my ( $sq_a, $sq_b, $naln, $nm_a, $nm_b );

    #
    ( $sq_a, $sq_b ) = @$_ra;
    foreach $naln (@$_ro) {
        ( $nm_a, $nm_b ) = split /:/o, $naln->[0], 2;
        ( "$nm_a" eq "$sq_a" && "$nm_b" eq "$sq_b" ) && do {
            @$_ra = ( $nm_a, $nm_b );
            return ( $T, $F );
        };
    };    # foreach $naln
          # swapping
    foreach $naln (@$_ro) {
        ( $nm_a, $nm_b ) = split /:/o, $naln->[0], 2;
        ( "$nm_a" eq "$sq_b" && "$nm_b" eq "$sq_a" ) && do {
            @$_ra = ( $nm_a, $nm_b );
            return ( $T, $T );
        };
    };    # foreach $naln
          #
    return ( $F, $F );
}    # search_alnname

sub search_seqname() {
    my ( $sq_a, $sq_b, $_ro ) = @_;
    my ( $a_found, $naln );
    foreach $naln (@$_ro) {
        ( $sq_a eq '.+' && defined($sq_b) && "$naln->[0]" eq "$sq_b" ) && next;
        ( $sq_a eq '.+' || "$naln->[0]" eq "$sq_a" ) && do {
            $sq_a = $naln->[0];
            return ( $T, $sq_a );
        };
    };    # foreach $naln
    return ( $F, undef );
}    # search_seqname

sub get_first_coord() {
    my ( $ary, $cct ) = @_;
    my ( $thelast, $cur, $tht );
    $thelast = $#{$ary};
    for ( $cur = 0 ; $cur < $thelast ; $cur++ ) {
        defined( $ary->[$cur] ) && do {
            ( $tht, $ary->[$cur] ) = ( $ary->[$cur], undef );
            return &get_units( \@$tht );
        };    # defined($ary[$cur])
    };    # for $cur
    $$cct++;
    return &get_units( \@{ $ary->[$thelast] } );
}    # get_first_coord

sub set_seq_boundaries() {
    my ($vrf) = @_;
    my ( $a, $x, $y, $rrf, %tcrds, $cc );
    %tcrds = (
        xs => [],
        xe => [],
        ys => [],
        ye => [],
    );
    ( $a, $x, $y ) = (
      $vrf->{_alignment_name}, $vrf->{_x_sequence_name},
      $vrf->{_y_sequence_name}
    );
    $vrf->{_plot_aln} && do {
        $rrf = \@{ $ALN_DATA{$a}[$_counter] };
        if ( $vrf->{_swap_aln} ) {
            push @{ $tcrds{xs} }, $rrf->[$_nori];
            push @{ $tcrds{xe} }, $rrf->[$_nend];
            push @{ $tcrds{ys} }, $rrf->[$_ori];
            push @{ $tcrds{ye} }, $rrf->[$_end];
        }
        else {
            push @{ $tcrds{xs} }, $rrf->[$_ori];
            push @{ $tcrds{xe} }, $rrf->[$_end];
            push @{ $tcrds{ys} }, $rrf->[$_nori];
            push @{ $tcrds{ye} }, $rrf->[$_nend];
        }
        ;    # if $vrf->{_swap_aln}
    };
    $vrf->{_plot_x_sequence} && do {
        $rrf = \@{ $GFF_DATA{$x}[$_counter] };
        push @{ $tcrds{xs} }, $rrf->[$_ori];
        push @{ $tcrds{xe} }, $rrf->[$_end];
    };
    $vrf->{_plot_y_sequence} && do {
        $rrf = \@{ $GFF_DATA{$y}[$_counter] };
        push @{ $tcrds{ys} }, $rrf->[$_ori];
        push @{ $tcrds{ye} }, $rrf->[$_end];
    };

    # print LOGFILE (Data::Dumper->Dump([ \%tcrds ], [ qw( *tcrds ) ]))
    #     if ($LogFile && $Debug);
    #
    $vrf->{_x_start} =
      [ ( scalar( @{ $tcrds{xs} } ) > 0 ? &min( @{ $tcrds{xs} } ) : 0 ), 'bp' ];
    $vrf->{_x_end} = [
        ( scalar( @{ $tcrds{xe} } ) > 0 ? &max( @{ $tcrds{xe} } ) : 1000 ),
        'bp'
    ];
    $vrf->{_y_start} =
      [ ( scalar( @{ $tcrds{ys} } ) > 0 ? &min( @{ $tcrds{ys} } ) : 0 ), 'bp' ];
    $vrf->{_y_end} = [
        ( scalar( @{ $tcrds{ye} } ) > 0 ? &max( @{ $tcrds{ye} } ) : 1000 ),
        'bp'
    ];

    #
    %tcrds = (
        xs => [
            defined( $vrf->{x_sequence_zoom} ) ? $vrf->{x_sequence_zoom}[0]
            : undef,
            defined( $vrf->{x_sequence_zoom_start} )
            ? $vrf->{x_sequence_zoom_start}
            : undef,
            defined( $vrf->{x_sequence_coords} ) ? $vrf->{x_sequence_coords}[0]
            : undef,
            defined( $vrf->{x_sequence_start} ) ? $vrf->{x_sequence_start}
            : undef, $vrf->{_x_start} ],
        xe => [
            defined( $vrf->{x_sequence_zoom} ) ? $vrf->{x_sequence_zoom}[1]
            : undef,
            defined( $vrf->{x_sequence_zoom_end} ) ? $vrf->{x_sequence_zoom_end}
            : undef,
            defined( $vrf->{x_sequence_coords} ) ? $vrf->{x_sequence_coords}[1]
            : undef,
            defined( $vrf->{x_sequence_end} ) ? $vrf->{x_sequence_end} : undef,
            $vrf->{_x_end} ],
        ys => [
            defined( $vrf->{y_sequence_zoom} ) ? $vrf->{y_sequence_zoom}[0]
            : undef,
            defined( $vrf->{y_sequence_zoom_start} )
            ? $vrf->{y_sequence_zoom_start}
            : undef,
            defined( $vrf->{y_sequence_coords} ) ? $vrf->{y_sequence_coords}[0]
            : undef,
            defined( $vrf->{y_sequence_start} ) ? $vrf->{y_sequence_start}
            : undef, $vrf->{_y_start} ],
        ye => [
            defined( $vrf->{y_sequence_zoom} ) ? $vrf->{y_sequence_zoom}[1]
            : undef,
            defined( $vrf->{y_sequence_zoom_end} ) ? $vrf->{y_sequence_zoom_end}
            : undef,
            defined( $vrf->{y_sequence_coords} ) ? $vrf->{y_sequence_coords}[1]
            : undef,
            defined( $vrf->{y_sequence_end} ) ? $vrf->{y_sequence_end} : undef,
            $vrf->{_y_end} ],
    );    # %tcrds
          # print LOGFILE (Data::Dumper->Dump([ \%tcrds ], [ qw( *tcrds ) ]))
          #     if ($LogFile && $Debug);
          #
    $vrf->{_x_start} =
      &get_units(
        defined( $vrf->{x_sequence_coords} )
        ? \@{ $vrf->{x_sequence_coords}[0] }
        : \@{ $vrf->{_x_start} } );
    $vrf->{_x_end} =
      &get_units(
        defined( $vrf->{x_sequence_coords} )
        ? \@{ $vrf->{x_sequence_coords}[1] }
        : \@{ $vrf->{_x_end} } );
    $vrf->{_y_start} =
      &get_units(
        defined( $vrf->{y_sequence_coords} )
        ? \@{ $vrf->{y_sequence_coords}[0] }
        : \@{ $vrf->{_y_start} } );
    $vrf->{_y_end} =
      &get_units(
        defined( $vrf->{y_sequence_coords} )
        ? \@{ $vrf->{y_sequence_coords}[1] }
        : \@{ $vrf->{_y_end} } );
    $vrf->{zoom_area} && do {
        $cc = 0;
        $vrf->{_x_zoomarea_start} = &get_first_coord( \@{ $tcrds{xs} }, \$cc );
        $vrf->{_x_zoomarea_end}   = &get_first_coord( \@{ $tcrds{xe} }, \$cc );
        $vrf->{_y_zoomarea_start} = &get_first_coord( \@{ $tcrds{ys} }, \$cc );
        $vrf->{_y_zoomarea_end}   = &get_first_coord( \@{ $tcrds{ye} }, \$cc );
        $cc == 4 && ( $vrf->{zoom_area} = $F, $vrf->{zoom} = $F );
    };    # zoom_area
    $vrf->{zoom} && do {
        $cc = 0;
        $vrf->{_x_zoom_start} = &get_first_coord( \@{ $tcrds{xs} }, \$cc );
        $vrf->{_x_zoom_end}   = &get_first_coord( \@{ $tcrds{xe} }, \$cc );
        $vrf->{_y_zoom_start} = &get_first_coord( \@{ $tcrds{ys} }, \$cc );
        $vrf->{_y_zoom_end}   = &get_first_coord( \@{ $tcrds{ye} }, \$cc );
        $cc == 4 && ( $vrf->{zoom} = $F );
        $vrf->{zoom} && do {
            $vrf->{_x_start} = $vrf->{_x_zoom_start};
            $vrf->{_x_end}   = $vrf->{_x_zoom_end};
            $vrf->{_y_start} = $vrf->{_y_zoom_start};
            $vrf->{_y_end}   = $vrf->{_y_zoom_end};
        };    # zoom...
    };    # zoom
    $vrf->{zoom} || ( $vrf->{zoom_marks} = $F );    # !zoom
    &report(
        'PLOTLIMITS',
        "X : $vrf->{_x_start} to $vrf->{_x_end}",
        "Y : $vrf->{_y_start} to $vrf->{_y_end}"
    );
    &report(
        ( $vrf->{zoom}
        ? (
              'ONZOOM',
              "X : $vrf->{_x_zoom_start} to $vrf->{_x_zoom_end}",
              "Y : $vrf->{_y_zoom_start} to $vrf->{_y_zoom_end}"
        )
        : 'NOZOOM' )
    );
    &report(
        ( $vrf->{zoom_area}
        ? (
              'ONZOOMAREA',
              "X : $vrf->{_x_zoomarea_start} to $vrf->{_x_zoomarea_end}",
              "Y : $vrf->{_y_zoomarea_start} to $vrf->{_y_zoomarea_end}"
        )
        : 'NOZOOMAREA' )
    );
    return;
}    # set_seq_boundaries

sub set_score_boundaries() {
    my ($vrf) = @_;
    my ( $mnsco, $mxsco );
    ( $vrf->{_aln_min_score}, $vrf->{_aln_max_score} ) = ( 0, 100 );
    $vrf->{_plot_aln} && do {
        $vrf->{_aln_min_score} =
          $ALN_DATA{ $vrf->{_alignment_name} }[$_counter][$_mnsco];
        $vrf->{_aln_max_score} =
          $ALN_DATA{ $vrf->{_alignment_name} }[$_counter][$_mxsco];
    };    # _plot_aln
    defined( $vrf->{aplot_score_range} ) && do {
        ( $mnsco, $mxsco ) = @{ $vrf->{aplot_score_range} };
        defined($mnsco) && ( $vrf->{_aln_min_score} = $mnsco );
        defined($mxsco) && ( $vrf->{_aln_max_score} = $mxsco );
    };    # defined($vrf->{aplot_score_range})
    ( $vrf->{_pbox_min_score}, $vrf->{_pbox_max_score} ) =
      ( $vrf->{_aln_min_score}, $vrf->{_aln_max_score} );
    ( $vrf->{show_percent_box} && defined( $vrf->{percent_box_score_range} ) )
      && do {
        ( $mnsco, $mxsco ) = @{ $vrf->{percent_box_score_range} };
        defined($mnsco) && ( $vrf->{_pbox_min_score} = $mnsco );
        defined($mxsco) && ( $vrf->{_pbox_max_score} = $mxsco );
    };    # show_percent_box
    &report(
        'PLOTSCORES',
        "ALIGNMENT:   $vrf->{_aln_min_score} to $vrf->{_aln_max_score}",
        "PERCENT box: $vrf->{_pbox_min_score} to $vrf->{_pbox_max_score}"
    );
}    # set_score_boundaries

sub set_tickmark_vars() {
    my ($vrf) = @_;
    my ( @tary, $stplen );
    $vrf->{show_only_bottom_ticks} && do {
        WHATICK: {
            $vrf->{show_extra_box} && ( @tary = ( $F, $F, $T ), last WHATICK );
            $vrf->{show_percent_box}
              && ( @tary = ( $F, $T, $F ), last WHATICK );
            @tary = ( $T, $F, $F );
        };    # WHATICK
        ( $vrf->{show_aplot_x_ticks}, $vrf->{show_percent_x_ticks},
          $vrf->{show_extrabox_x_ticks} )
          = @tary;
    };
    $vrf->{_x_len} = $vrf->{_x_end} - $vrf->{_x_start} + 1;
    $vrf->{_y_len} = $vrf->{_y_end} - $vrf->{_y_start} + 1;
    if ( defined( $vrf->{aplot_major_tickmark} )
        && !defined( $vrf->{major_tickmark_nucleotide} ) )
    {
        $stplen =
          10**( int( &logdec( &min( $vrf->{_x_len}, $vrf->{_y_len} ) ) ) + 1 ) /
          $vrf->{aplot_major_tickmark};
    }
    else {
        $stplen =
          10**( int( &logdec( &min( $vrf->{_x_len}, $vrf->{_y_len} ) ) - 0.175 )
        );
    }
    defined( $vrf->{major_tickmark_nucleotide} ) || do {
        $vrf->{major_tickmark_nucleotide} = $stplen;
    };    # defined($vrf->{major_tickmark_nucleotide})
    defined( $vrf->{minor_tickmark_nucleotide} ) || do {
        $vrf->{minor_tickmark_nucleotide} =
          $vrf->{major_tickmark_nucleotide} / $vrf->{aplot_minor_tickmark};
    };    # defined($vrf->{minor_tickmark_nucleotide})
    $vrf->{_pbox_len} = $vrf->{_pbox_max_score} - $vrf->{_pbox_min_score} + 1;
    if ( defined( $vrf->{percent_major_tickmark} )
        && !defined( $vrf->{major_tickmark_score} ) )
    {
        $stplen =
          10**( int( &logdec( $vrf->{_pbox_len} * 1000 ) ) - 2 ) /
          $vrf->{percent_major_tickmark};
    }
    else {
        $stplen = 10**( int( &logdec( $vrf->{_pbox_len} * 1000 ) - 3 ) );
    }
    defined( $vrf->{major_tickmark_score} ) || do {
        $vrf->{major_tickmark_score} = $stplen;
    };    # defined($vrf->{major_tickmark_score})
    defined( $vrf->{minor_tickmark_score} ) || do {
        $vrf->{minor_tickmark_score} =
          $vrf->{major_tickmark_score} / $vrf->{percent_minor_tickmark};
    };    # defined($vrf->{minor_tickmark_score})
    &report( 'TICKDONE',
        "Major step $vrf->{major_tickmark_nucleotide} - "
        . "Minor step $vrf->{minor_tickmark_nucleotide}",
        "Major step $vrf->{major_tickmark_score} - "
        . "Minor step $vrf->{minor_tickmark_score}" );
    return;
}    # set_tickmark_vars

sub find_prop_ratio() {
    my ($vrf) = @_;
    RATIO: {
        defined( $vrf->{aplot_xy_scale} ) && do {
            ( $vrf->{aplot_xy_scale} > 0 && $vrf->{aplot_xy_scale} < 10 )
              && do {
                $vrf->{aplot_xy_same_length} = $F;
                last RATIO;
            };
        };    # {aplot_xy_scale}
        $vrf->{aplot_xy_same_length} && do {
            $vrf->{aplot_xy_scale} = -1;
            last RATIO;
        };    # {aplot_xy_same_length}
        $vrf->{aplot_xy_scale} =
          sprintf( "%5.4f", ( $vrf->{_x_len} / $vrf->{_y_len} ) );
    };    # RATIO
    &report( 'SCALEDONE',
        ( ( $vrf->{aplot_xy_scale} ne -1 )
        ? "Current plot has an XY-RATIO of $vrf->{aplot_xy_scale} ..."
        : "XY-RATIO is not taken into account for this plot ..." ),
        ( $vrf->{aplot_xy_same_length} ? 'YES' : 'NO' ) );
    return;
}    # find_prop_ratio

sub set_misc_vars() {
    my ($vrf) = @_;
    my %lbllen;
    %lbllen = (
        xfl => $vrf->{feature_x_label_length},
        xgl => $vrf->{group_x_label_length},
        yfl => $vrf->{feature_y_label_length},
        ygl => $vrf->{group_y_label_length},
    );
    $vrf->{_plot_x_sequence} && do {
        defined( $lbllen{xfl} ) || do {
            $lbllen{xfl} =
              $GFF_DATA{ $vrf->{_x_sequence_name} }[$_counter][$_flw];
        };    # feature_x_label_length
        defined( $lbllen{xgl} ) || do {
            $lbllen{xgl} =
              $GFF_DATA{ $vrf->{_x_sequence_name} }[$_counter][$_glw];
        };    # group_x_label_length
    };    # _plot_x_sequence
    $vrf->{_plot_y_sequence} && do {
        defined( $lbllen{yfl} ) || do {
            $lbllen{yfl} =
              $GFF_DATA{ $vrf->{_y_sequence_name} }[$_counter][$_flw];
        };    # feature_y_label_length
        defined( $lbllen{ygl} ) || do {
            $lbllen{ygl} =
              $GFF_DATA{ $vrf->{_y_sequence_name} }[$_counter][$_glw];
        };    # feature_y_label_length
    };    # _plot_y_sequence
    foreach my $u ( keys %lbllen ) {
        $lbllen{$u} =
          ( defined( $lbllen{$u} ) && $lbllen{$u} >= 0 ) ? $lbllen{$u} : 0;
    };    # foreach $u
    $vrf->{feature_x_label_length} = $lbllen{xfl};
    $vrf->{group_x_label_length}   = $lbllen{xgl};
    $vrf->{feature_y_label_length} = $lbllen{yfl};
    $vrf->{group_y_label_length}   = $lbllen{ygl};
}    # set_misc_vars

sub set_page_labels() {
    my ($vrf) = @_;
    defined( $vrf->{title} ) || do {
        ( $vrf->{title} = $vrf->{_alignment_name} ) =~ s/:/ x /o;
    };    # defined($vrf->{title})
    defined( $vrf->{subtitle} ) || ( $vrf->{subtitle} = "" );
    defined( $vrf->{x_label} ) || do {
        $vrf->{x_label} = $vrf->{_x_sequence_name};
    };    # defined($vrf->{x_label})
    defined( $vrf->{y_label} ) || do {
        $vrf->{y_label} = $vrf->{_y_sequence_name};
    };    # defined($vrf->{y_label})
    defined( $vrf->{percent_box_label} ) || ( $vrf->{percent_box_label} = '' );
    defined( $vrf->{percent_box_sublabel} )
      || ( $vrf->{percent_box_sublabel} = '' );
    defined( $vrf->{extra_box_label} ) || ( $vrf->{extra_box_label} = '' );
    defined( $vrf->{extra_box_sublabel} )
      || ( $vrf->{extra_box_sublabel} = '' );
    $vrf->{show_ps_warnings} && do {
        $vrf->{_plot_x_sequence}
          || ( $vrf->{x_label} .= " ... X-Sequence Annotation NOT FOUND ... " );
        $vrf->{_plot_y_sequence}
          || ( $vrf->{y_label} .= " ... Y-Sequence Annotation NOT FOUND ... " );
        $vrf->{_plot_aln} || do {

            # PS --> angle string valn haln {fz fn} tcol x y  xtxt
            my ( $epo, $epe ) =
              ( ( $vrf->{_x_end} - $vrf->{_x_start} ) / 2,
              ( $vrf->{_y_end} - $vrf->{_y_start} ) / 2, );
            push @{ $SpecialCache{POST} },
              [ $epo, $epo, $epe, $epe,
                '[ 45 (GFF-record Alignment Features NOT FOUND) '
                . '(bv) (ch) 18 pt /Helvetica-Bold _f '
                . ( $epo - 10 ) . ' '
                . ( $epe + 10 ) . ' xtxt' ];
            push @{ $SpecialCache{POST} },
              [ $epo, $epo, $epe, $epe,
                '[ 45 (...Please, check your command-line and/or GFF '
                . 'input files...) (tv) (ch) 15 pt /Helvetica _f '
                . ( $epo + 10 ) . ' '
                . ( $epe - 10 ) . ' xtxt' ];
        };    # _plot_aln
    };    # $vrf->{show_ps_warnings}
    return;
}    # set_page_labels

sub make_plot() {
    &header("WRITING POSTSCRIPT TO STDOUT");

    &ps_header;
    &ps_colors;
    &ps_page_formats;
    &ps_variables;
    &ps_main;

    &ps_plot;

    &ps_trailer;

    &footer("WRITING POSTSCRIPT FINISHED");
}    # make_plot

sub prt_help() {
    open( HELP, "| more" );
    print HELP <<"+++EndOfHelp+++";
PROGRAM:
                        $PROGRAM - $VERSION

    Converting GFF files for pairwise alignments to PostScript.

USAGE:        $PROGRAM [options] <GFF_files|STDIN>

DESCRIPTION:

    This program draws color-filled alignment plots from GFF
    files for that alignment and two sequences annotations.

REQUIRES:

    $PROGRAM needs the following Perl modules installed in 
    your system, we used those available from the standard 
    Perl distribution. Those that are not in the standard 
    distribution are marked with an '(*)', in such cases 
    make sure that you already have downloaded them from 
    CPAN (http://www.perl.com/CPAN) and installed.

    "Getopt::Long" - processing command-line options.
    "Data::Dumper" - pretty printing data structures for debugging (*).
    "Benchmark" - checking and comparing running times of code.

ENVIRONMENT VARIABLES:

        There are two environmental variables that can be set by 
    users to their preferences:
     + You can specify the path where $PROGRAM can find the default
      files with the shell variable \"GFF2APLOT_CUSTOMPATH\". Default
      value is the path where you are running $PROGRAM.
     + You can also define the default custom filename you will like
      with the variable \"GFF2APLOT_CUSTOMFILE\", program default
      filename for custom file is \".gff2aplotrc\".
     + Now $PROGRAM does not need to write any temporary file, 
      so that previous versions default temporary directory path
      variable (\"GFF2APLOT_TMP\") is no longer used.
     + Setting those vars in Bourne-shell and C-shell:
       o Using a Bourne-Shell (e.g. bash):
            export GFF2APLOT_CUSTOMPATH=\"path\"
            export GFF2APLOT_CUSTOMFILE=\"file_name\"
       o Using a C-Shell:
            setenv GFF2APLOT_CUSTOMPATH \"path\"
            setenv GFF2APLOT_CUSTOMFILE \"file_name\"

COMMAND-LINE OPTIONS:

    A double dash on itself "--" signals end of the options
    and start of file names (if present). You can use a single
    dash "-" as STDIN placeholder. Available options and a
    short description are listed here:

    -h, --help   Shows this help.
    --version    Shows current version and exits.
    -v, --verbose
          Verbose mode, a full report is sent to standard error 
          (default is set to showing only WARNINGS).
    -V, --logs-filename  <filename>
          Report is written to a log file.
    -q, --quiet
          Quiet mode, do not show any message/warning
          to standard error (reporting only ERRORS).
    -P, --page-bbox  <width,height>
          Setting a user-defined page size, <width> and
          <height> are set to points if no unit is given,
          you can use pt, mm, cm or in. This option overrides ANY
          'page-size' definition (from command-line or custom file).
    -p, --page-size  <format_name>
          Setting a page size among pre-defined ones
          (see below for a list of available page formats).
    --margin-left   <length>
    --margin-right  <length>
    --margin-top    <length>
    --margin-bottom <length>
          Setting each page margin to <length>. If no units are provided,
          points are assumed, you can use points, milimeters, centimeters
          or inches (pt, mm, cm or in, respectively).
    -B, --background-color  <color>    Background color.
    -F, --foreground-color  <color>    Foreground color.
    -T, --title <string>   Definning plot title.
    -t, --subtitle <string>   Definning plot subtitle.
    -X, --x-label <string>   Defining X-axis label.
    -Y, --y-label <Y-Label>   Defining Y-axis label.
    -L, --percent-box-label <string>   Defining percent-box label.
    -l, --extra-box-label <XBox-Label>   Definning Extra-Box Label.
    -x, --x-sequence-coords <pos..pos> 
    -S, --start-x-sequence <pos>  Sets X-sequence first nucleotide.
    -E, --end-x-sequence <pos>    Sets X-sequence last nucleotide.
    -y, --y-sequence-coords <pos..pos> 
    -s, --start-y-sequence <pos>  Sets Y-sequence first nucleotide.
    -e, --end-y-sequence <pos>    Sets Y-sequence last nucleotide.
    --x-sequence-zoom <pos..pos> 
    --y-sequence-zoom <pos..pos> 
    -Z, --zoom [ [-S <pos>] [-E <pos>] [-s <pos>] [-e <pos>] ]
                   This option zooms an area you have selected
                   with -S,-E,-s,-e (all 4 are optional).
    -z, --zoom-area [ [-S <pos>] [-E <pos>] [-s <pos>] [-e <pos>] ]
                   This option marks a zoom area on your plot,
                   but does not make a zoom.
    -A, --alignment-name <SeqXName:SeqYName>
         Defining which alignment is going to be plotted 
         if you have more than one alignment in your gff input.
    -N, --x-sequence-name <SeqXName>
         Defining which sequence is going to be plotted at X-axis.
    -n, --y-sequence-name <SeqYName>
         Defining which sequence is going to be plotted at Y-axis.
    -r, --aplot-xy-noteq
          By default X and Y axes have same length, this option 
          disables such behaviour, so X and Y sequence will have
          axes-lengths proportional to their nucleotide lengths.
    -R, --xy-axes-scale  <X/Y ratio>
          This option allows to set a different scale between X
          and Y axes lengths (by default is '1'). Below 1 values
          make Y larger than X, and larger than 1 result in getting
          X larger than Y. # Must be explained better.
    -W, --aln-scale-width   Scaling score on width for Aplot lines.
    -w, --aln-scale-color   Scaling score on color for Aplot lines.
    -K, --show-ribbons <ribbon_type>
          Force Ribbons for all features on axes to be:
             (N)one, (L)ines, (R)ibbons, (B)oth.
    -G, --show-grid
          Switches 'on' grid (default is 'off').
    -g, --hide-grid
          Switches 'off' grid (if switched on from customization files).
    -I, --show-percent-box
          Switches 'on' Percent box (default is 'off').
    -i, --hide-percent-box
          Switches 'off' Percent box (if set to 'on' on custom files).
    -O, --show-extra-box
          Switches 'on' Extra box (default is 'off').
    -o, --hide-extra-box
          Switches 'off' Extra box (if set to 'on' on custom files).
    -D, --aplot-box-color <color>   Aplot main box background color.
    -d, --percent-box-color <color>   Percent box background color.
    -b, --extra-box-color <color>   Extra box background color.
    --nopswarnings
          Switch off warnings that may appear on the final PostScript figure
          when X sequence, Y sequence and/or alignment input data is missing.
    -a, --hide-credits
          Switch off $PROGRAM credits line on plot.
    --debug
        Reporting variable contents when testing the program.
        Requires that log report file option was also activated.
      --layout-var '<variable=value>'
    --sequence-var '<sequence::variable=value>'
      --source-var '<source::variable=value>'
      --strand-var '<strand::variable=value>' 
       --group-var '<group::variable=value>'
     --feature-var '<feature::variable=value>'
             Loading a feature/group/strand/source/sequence/layout
             customization variable from command-line. You can set
             several variables by repeating any of these options, 
             i.e.:
               ... --feature-var 'cds::feature_shape=box' \
                   --feature-var 'cds::feature_color=blue' ... 
    -C, --custom-filename <filename>
          Loading customization parameters from a given file (if 
          default \".gff2aplotrc\" exists is loaded before it).

    Those are the colors defined in $PROGRAM:
    + Default Colors: black white.
    + Color Palette: 
          grey pink magenta violet blue skyblue cyan seagreen
             green limegreen yellow orange red brown
      You can get up to nine color shades from Variable Colors with
        \"verydeep\", \"deep\", \"verydark\", \"dark\",
        \"light\", \"verylight\", \"pale\" and \"verypale\" prefixes,
      as example: 
        verydeepblue, deepblue, verydarkblue, darkblue,
        blue, lightblue, verylightblue, paleblue and verypaleblue.
    + Dynamic colors: gradient rainbow.

    The following page sizes are available: from A0 to A10, 
    from B0 to B10, 10x14, executive, folio, ledger, legal, 
    letter, quarto, statement and tabloid.

BUGS:    Report any problem to 'jabril\@imim.es'.

AUTHORS:

            $PROGRAM is under GNU-GPL
         Josep Francesc ABRIL FERRANDO  
                 Thomas WIEHE                   
                Roderic GUIGO SERRA       

            Copyright (C) 1999-2003

+++EndOfHelp+++
    close(HELP);
    exit(1);
}    # prt_help
     #
     # GENERAL FUNCTIONS
     #
sub close_logfile() { close(LOGFILE) if $LogFile };
sub get_units() { return $UNITS{ $_[0]->[1] }->( $_[0]->[0] ); }
sub tobool()    { return ( $_[0] ? '_t' : '_f' ); }

sub tostring() {
    ( defined( $_[0] ) && $_[0] ne '' ) || return '()';
    $_[0] =~ s{[\(]}{\\050}og;
    $_[0] =~ s{[\)]}{\\051}og;
    return '(' . $_[0] . ')';
};    # tostring
sub logdec() { return ( log( $_[0] ) / log(10) ); }

sub trap_signals() {
    &prt_to_logfile( $Messages{'USER_HALT'} );
    &close_logfile();
    die ( $Messages{'USER_HALT'} );
}    # trap_signals

sub trap_signals_prog() {
    &prt_to_logfile( $Messages{'PROCESS_HALT'} );
    &close_logfile();
    die ( $Messages{'PROCESS_HALT'} );
}    # trap_signals_prog

sub warn() {
    my $type       = shift @_;
    my $screen_flg = shift @_;
    my $comment    = sprintf( $Messages{$type}, @_ );

    # ALWAYS to STDERR if $screen_flg==$T unless $Quiet==$T
    $screen_flg && ( $Quiet || print STDERR $comment );
    &prt_to_logfile($comment);
}    # warn
sub prt_to_logfile() { $LogFile && ( print LOGFILE $_[0] ) }
sub prt_to_stderr()  { $Verbose && ( $Quiet || print STDERR $_[0] ) }

sub report() {
    my $type = shift @_;
    my $comment = sprintf( $Messages{$type}, @_ );
    &prt_to_stderr($comment);
    &prt_to_logfile($comment);
}    # report

sub header() {
    my $comment = $line;
    foreach my $ln (@_) {
        $comment .= "### " . &fill_mid( "$ln", 72, " " ) . " ###\n";
    }
    $comment .= $line;
    &prt_to_stderr($comment);
    &prt_to_logfile($comment);
}    # header

sub footer() {
    $total_time = &timing($F);
    &header( @_, $total_time );
    &prt_to_stderr("###\n");
    &prt_to_logfile("###\n");
}    # footer

sub timing() {
    push @Timer, ( new Benchmark );

    # partial time 
    $_[0]
      || (
        return timestr( timediff( $Timer[$#Timer], $Timer[ ( $#Timer - 1 ) ] ) )
    );

    # total time
    return timestr( timediff( $Timer[$#Timer], $Timer[0] ) );
}    # timing
     #

sub max() {
    my $z = shift @_;
    foreach my $l (@_) { $z = $l if $l > $z }
    return $z;
}    # max

sub min() {
    my $z = shift @_;
    foreach my $l (@_) { $z = $l if $l < $z }
    return $z;
}    # min
     #
sub fill_right() { $_[0] . ( $_[2] x ( $_[1] - length( $_[0] ) ) ) }
sub fill_left() { ( $_[2] x ( $_[1] - length( $_[0] ) ) ) . $_[0] }

sub fill_mid() {
    my $l = length( $_[0] );
    my $k = int( ( $_[1] - $l ) / 2 );
    ( $_[2] x $k ) . $_[0] . ( $_[2] x ( $_[1] - ( $l + $k ) ) );
}    # fill_mid
     #

sub counter {    # $_[0]~current_pos++ $_[1]~char
    my $str;
    $str = "$_[1]";
    ( ( $_[0] % 50 ) == 0 )
      && ( $str .= " [" . &fill_left( $_[0], 6, "0" ) . "]\n" );
    &prt_to_stderr($str);
    &prt_to_logfile($str);
}    # counter
     #

sub counter_end {    # $_[0]~current_pos   $_[1]~char
    my $str;
    ( ( $_[0] % 50 ) != 0 ) && do {
        $str = " [" . &fill_left( $_[0], 6, "0" ) . "]\n";
        &prt_to_stderr($str);
        &prt_to_logfile($str);
      }
}    # counter_end
     #
     # POSTSCRIPT CODE
     #

sub ps_plot() {
    &report('PSPLOT');

    #
    &ps_open_page;

    #
    &ps_block_aplot;

    #
    $Vars{LAYOUT}{show_percent_box} && &ps_block_percent;

    #
    $Vars{LAYOUT}{show_extra_box} && &ps_block_extra;

    #
    &ps_close_page;

    #
    &report('PSPLOTDONE');
}    # ps_plot

sub ps_block_aplot() {
    my $vlr = \%{ $Vars{LAYOUT} };

    #
    &report( 'BLOCKTODO', 'APLOT-BOX' );

    #
    $StringCache{ $vlr->{_x_sequence_name} } = {} if $vlr->{_plot_x_sequence};
    $StringCache{ $vlr->{_y_sequence_name} } = {} if $vlr->{_plot_y_sequence};
    $StringCache{ $vlr->{_alignment_name} } = '' if $vlr->{_plot_aln};

    # zoom area
    $vlr->{zoom_area} && do {
        my %kt;
        $kt{color} =
          defined( $vlr->{zoom_area_fill_color} )
          ? "$vlr->{zoom_area_fill_color} " . &tobool($T)
          : &tobool($F);
        $kt{xo} = &max( $vlr->{_x_zoomarea_start}, $vlr->{_x_start} );
        $kt{xe} = &min( $vlr->{_x_zoomarea_end},   $vlr->{_x_end} );
        $kt{yo} = &max( $vlr->{_y_zoomarea_start}, $vlr->{_y_start} );
        $kt{ye} = &min( $vlr->{_y_zoomarea_end},   $vlr->{_y_end} );
        unshift @{ $SpecialCache{PRE} },
          [
            $kt{xo},
            $kt{xe},
            $kt{yo},
            $kt{ye},
            join (
                " ",
                $vlr->{zoom_area_mark_color},
                &get_units( $vlr->{zoom_area_mark_width} ),
                $vlr->{zoom_area_mark_style},
                $kt{color},
                $kt{xo},
                $kt{xe},
                $kt{yo},
                $kt{ye}
            ) . " xbox"
        ];
    };    # zoom_area
          #
    print STDOUT << '+++MAINProcs+++';
%
%%%%%%%% START NEW PLOT
%
DoInit
DoHeader
%
%%%%%%%% ALIGNMENT PLOT - BOX
%
begindata
+++MAINProcs+++

    # tag definitions: Box,Join,Arrow,Banner,Ribbons;
    &ps_plot_features( $vlr->{_x_sequence_name}, $T )
      if $vlr->{_plot_x_sequence};
    &ps_plot_features( $vlr->{_y_sequence_name}, $F )
      if $vlr->{_plot_y_sequence};

    # nice ribbon lines finishing
    &ps_plot_outlines( $vlr->{_x_sequence_name}, $T )
      if $vlr->{_plot_x_sequence};
    &ps_plot_outlines( $vlr->{_y_sequence_name}, $F )
      if $vlr->{_plot_y_sequence};

    # extras: PRE (under aln features)
    &ps_plot_addons('PRE') if scalar( @{ $SpecialCache{PRE} } ) > 0;

    # alignment
    &ps_plot_aln( $vlr->{_alignment_name} ) if $vlr->{_plot_aln};

    # extras: POST (over aln features)
    &ps_plot_addons('POST') if scalar( @{ $SpecialCache{POST} } ) > 0;

    #
    print STDOUT "\%\nenddata\n\%\n";
    &report( 'BLOCKDONE', 'APLOT-BOX' );
}    # ps_block_aplot

sub ps_plot_features() {
    my ( $seq, $flg ) = @_;
    my (
        $src,  $str,   $grp,  $qlst, $slst, $nlst, $glst,
        $elst, @rflst, %this, $vaf,  $go,   $orf,  $what
    );
    $vaf = \%{ $Vars{LAYOUT} };
    ( $go, $orf ) = &get_ordered_ary( 'GFF', $seq );
    $what = $flg ? 'X-sequence' : 'Y-sequence';
    &report( 'BLOCKXY', $what, $seq );
    $go && do {
        print STDOUT "\%\n\%\%\%\% DATA ---> $what Annotations\n\%\n"
          . ( $flg ? 'beginXseq' : 'beginYseq' ) . "\n";
        $qlst = \@{ $GFF_DATA{$seq} };
        for ( my $a = 0 ; $a <= $#{$orf} ; $a++ ) {
            $src  = $orf->[$a][0];
            $slst = \@{ $qlst->[$_element]{$src} };
            &show_element($slst) || next;    # source:{hide}
            for ( my $b = 0 ; $b <= $#{ $orf->[$a][1] } ; $b++ ) {
                $str  = $orf->[$a][1][$b][0];
                $nlst = \@{ $slst->[$_element]{$str} };
                &show_element($nlst) || next;    # strand:{hide}
                for ( my $c = 0 ; $c <= $#{ $orf->[$a][1][$b][1] } ; $c++ ) {
                    $grp  = $orf->[$a][1][$b][1][$c];
                    $glst = \@{ $nlst->[$_element]{$grp} };
                    &show_element($glst) || next;    # group:{hide}
                    $this{go}    = $glst->[$_counter][$_ori];
                    $this{ge}    = $glst->[$_counter][$_end];
                    $this{sflag} = $flg
                      ? &check_seqlim(
                        $this{go},        $this{ge},
                        $vaf->{_x_start}, $vaf->{_x_end}
                      )
                      : &check_seqlim(
                        $this{go},        $this{ge},
                        $vaf->{_y_start}, $vaf->{_y_end}
                    );
                    $this{sflag} && do {
                        for ( my $d = 0 ; $d < $glst->[$_counter][$_elemNum] ;
                            $d++ )
                        {
                            $elst = \@{ $glst->[$_element][$d] };
                            &show_element($elst) || next;    # feature:{hide}
                            $this{so}    = $elst->[$_ftori];
                            $this{se}    = $elst->[$_ftend];
                            $this{sflag} = $flg
                              ? &check_seqlim(
                                $this{so},        $this{se},
                                $vaf->{_x_start}, $vaf->{_x_end}
                              )
                              : &check_seqlim(
                                $this{so},        $this{se},
                                $vaf->{_y_start}, $vaf->{_y_end}
                            );
                            $this{sflag} && do {
                                @rflst = ( $elst, $glst, $nlst, $slst, $qlst );

                                #
                                ( $this{strnd}, undef ) = split //o, $str, 2;

                                # setting score
                                $this{sco} = &norm_score(
                                    $elst->[$_ftsco],
                                    $qlst->[$_counter][$_mnsco],
                                    $qlst->[$_counter][$_mxsco]
                                );

# {feature_color} {feature_scale_color} {feature_scale_heigth}
                                $this{col} =
                                  &get_var_value( 'feature_color', \@rflst );
                                $this{colscl} =
                                  &tobool(
                                    &get_var_value( 'feature_scale_color',
                                        \@rflst ) );
                                $this{col} =~ /$gradregex/o && do {
                                    $this{col}    = "$this{sco} $this{col}";
                                    $this{colscl} = &tobool($F);
                                };
                                $this{col} .= " "
                                  . $this{colscl} . " "
                                  . &tobool(
                                    &get_var_value( 'feature_scale_heigth',
                                        \@rflst ) );

                                # {feature_shape}
                                $this{shape} =
                                  $Shapes{ &get_var_value( 'feature_shape',
                                      \@rflst ) };

                                # {feature_label} {show_feature_label}
                                $this{lbl_tmp} =
                                  &get_var_value( 'feature_label', \@rflst );
                                $this{lbl_flg} =
                                  &get_var_value( 'show_feature_label',
                                    \@rflst );
                                $this{lbl} =
                                  $this{lbl_flg}
                                  ? ( ( defined( $elst->[$_ftid] )
                                    && $elst->[$_ftid] ne '' )
                                    ? "($elst->[$_ftid]) _t"
                                    : ( ( defined( $this{lbl_tmp} )
                                        && $this{lbl_tmp} ne '' )
                                        ? "($this{lbl_tmp}) _t"
                                        : "(" . ( $d + 1 ) . ") _t" ) )
                                  : "_f";    # defined($this{lbl_tmp})
                                             # {show_ribbons} {ribbon_style}
                                @{ $this{rbnsty} } =
                                  @{ $Ribbons{ &get_var_value( 'ribbon_style',
                                      \@rflst ) } };
                                ( !$this{rbnsty}[0] && $this{rbnsty}[1] )
                                  && ( $this{rbnflg} = $T );
                                $this{rbnflg} =
                                  &get_var_value( 'show_ribbons', \@rflst );
                                ( !$this{rbnsty}[0] && !$this{rbnsty}[1] )
                                  && ( $this{rbnflg} = $F );

                                #
                                print STDOUT
                                  "\/$this{shape} $this{col} $this{lbl} "
                                  . "$this{sco} $this{so} $this{se} ($this{strnd}) shp\n";

                                # {ribbon_color}
                                $this{rbnflg} && do {
                                    $this{col} =
                                      &get_var_value( 'ribbon_color', \@rflst );
                                    $this{rbnsty}[0] && do {
                                        defined( $StringCache{$seq}{rbn} )
                                          || ( $StringCache{$seq}{rbn} = '' );
                                        $StringCache{$seq}{rbn} .=
"$this{so} $this{se} $this{col} rbx\n";
                                    };    # style:show_ribbons
                                    $this{rbnsty}[1] && do {
                                        defined( $StringCache{$seq}{lns} )
                                          || ( @{ $StringCache{$seq}{lns} } =
                                            () );
                                        push @{ $StringCache{$seq}{lns} },
                                          $this{so}, $this{se};
                                    };    # style:show_ribbons
                                };    # show_ribbons
                            };    # check_nuclim for features
                        };    # for $d
                        @rflst = ( $glst, $nlst, $slst, $qlst );

                        #
                        ( $this{strnd}, undef ) = split //o, $str, 2;

                        # {group_color} {group_shape} {show_group_limits}
                        $this{col} = &get_var_value( 'group_color', \@rflst );
                        $this{shape} =
                          $GroupShapes{ &get_var_value( 'group_shape', \@rflst )
                          };
                        $this{grid_flg} =
                          &tobool(
                            &get_var_value( 'show_group_limits', \@rflst ) );

                        # {group_label} {show_group_label}
                        $this{lbl_tmp} =
                          &get_var_value( 'group_label', \@rflst );
                        $this{lbl_flg} =
                          &get_var_value( 'show_group_label', \@rflst );
                        $this{lbl} =
                          $this{lbl_flg}
                          ? (
                            defined( $this{lbl_tmp} ) 
                            ? "($this{lbl_tmp}) _t"
                            : "($grp) _t" )
                          : "_f";

                        #
                        print STDOUT
"\/$this{shape} $this{grid_flg} $this{col} $this{lbl} "
                          . "$this{go} $this{ge} ($this{strnd}) grp\n";

                        #
                    };    # check_nuclim for groups
                };    # for $c
            };    # for $b
        };    # for $a
        print STDOUT "endseq\n";
        &ps_plot_ribbons( $seq, $flg );
    };    # $go
    $go || &report( 'NOGFFDATA', $what );
}    # ps_plot_features

sub ps_plot_addons() {
    my $vky = shift;
    my $vsf = \%{ $Vars{LAYOUT} };
    my %they;
    scalar( @{ $SpecialCache{$vky} } ) > 0 || return;
    &report('BLOCKFX');
    print STDOUT "\%\n\%\%\%\% Extra features from custom files ($vky)...\n"
      . "\%\nbeginaln\n";
    foreach my $ilm ( @{ $SpecialCache{$vky} } ) {
        ( $they{xo}, $they{xe}, $they{yo}, $they{ye}, $they{string} ) = @{$ilm};
        ( &check_seqlim( $they{xo}, $they{xe}, $vsf->{_x_start},
              $vsf->{_x_end} )
          && &check_seqlim( $they{yo}, $they{ye}, $vsf->{_y_start},
              $vsf->{_y_end} ) ) && do {
            print STDOUT "$they{string}\n";
        };    # check_nuclim
    };    # foreach $ilm
    print STDOUT "endseq\n";
}    # ps_plot_addons

sub ps_plot_aln() {
    my ($aln) = @_;
    my (
        $src,  $str,   $grp,  $qlst, $slst, $nlst, $glst,
        $elst, @rflst, %this, $vaf,  $go,   $orf
    );
    $vaf = \%{ $Vars{LAYOUT} };
    ( $go, $orf ) = &get_ordered_ary( 'ALN', $aln );
    $StringCache{$aln} = '';
    &report( 'BLOCKA', $aln );
    $go && do {
        print STDOUT "\%\n\%\%\%\% DATA ---> Alignment\n\%\nbeginaln\n";
        $qlst = \@{ $ALN_DATA{$aln} };
        for ( my $a = 0 ; $a <= $#{$orf} ; $a++ ) {
            $src  = $orf->[$a][0];
            $slst = \@{ $qlst->[$_element]{$src} };
            &show_element($slst) || next;    # source:{hide}
            for ( my $b = 0 ; $b <= $#{ $orf->[$a][1] } ; $b++ ) {
                $str  = $orf->[$a][1][$b][0];
                $nlst = \@{ $slst->[$_element]{$str} };
                &show_element($nlst) || next;    # strand:{hide}
                for ( my $c = 0 ; $c <= $#{ $orf->[$a][1][$b][1] } ; $c++ ) {
                    $grp  = $orf->[$a][1][$b][1][$c];
                    $glst = \@{ $nlst->[$_element]{$grp} };
                    &show_element($glst) || next;    # group:{hide}
                    for ( my $d = 0 ; $d < $glst->[$_counter][$_elemNum] ;
                        $d++ )
                    {
                        $elst = \@{ $glst->[$_element][$d] };
                        &show_element($elst) || next;    # feature:{hide}
                        if ( $vaf->{_swap_aln} ) {
                            $this{xo} = $elst->[$_ftnori];
                            $this{xe} = $elst->[$_ftnend];
                            $this{yo} = $elst->[$_ftori];
                            $this{ye} = $elst->[$_ftend];
                        }
                        else {
                            $this{xo} = $elst->[$_ftori];
                            $this{xe} = $elst->[$_ftend];
                            $this{yo} = $elst->[$_ftnori];
                            $this{ye} = $elst->[$_ftnend];
                        }
                        $this{aflag} = (
                          &check_seqlim(
                              $this{xo},        $this{xe},
                              $vaf->{_x_start}, $vaf->{_x_end}
                          ) && &check_seqlim(
                              $this{yo},        $this{ye},
                              $vaf->{_y_start}, $vaf->{_y_end}
                          )
                        );
                        $this{aflag} && do {
                            ( $this{xn}, $this{yn} ) = ( split //o, $str );
                            defined( $this{yn} ) || ( $this{yn} = $this{xn} );
                            $this{xn} eq '-' && do {
                                ( $this{xo}, $this{xe} ) =
                                  ( $this{xe}, $this{xo} );
                            };
                            $this{yn} eq '-' && do {
                                ( $this{yo}, $this{ye} ) =
                                  ( $this{ye}, $this{yo} );
                            };
                            @rflst = ( $elst, $glst, $nlst, $slst, $qlst );

# {alignment_color} {alignment_scale_color} {alignment_scale_width}
                            $this{col} =
                              &get_var_value( 'alignment_color', \@rflst ) . " "
                              . &tobool(
                                &get_var_value( 'alignment_scale_color',
                                    \@rflst ) ) . " "
                              . &tobool(
                                &get_var_value( 'alignment_scale_width',
                                    \@rflst ) );
                            $this{sco} = &norm_score(
                                $elst->[$_ftsco],
                                $vaf->{_aln_min_score},
                                $vaf->{_aln_max_score}
                            );
                            print STDOUT
                              "$this{xo} $this{yo} $this{xe} $this{ye}"
                              . " $this{col} $this{sco} aln\n";
                            $vaf->{show_percent_box} && do {
                                $this{ysc} = $elst->[$_ftsco];
                                &check_inlim(
                                    $this{ysc},
                                    $vaf->{_pbox_min_score},
                                    $vaf->{_pbox_max_score}
                                  ) && do {
                                    $StringCache{$aln} .=
"$this{xo} $this{ysc} $this{xe} $this{ysc}"
                                      . " $this{col} $this{sco} aln\n";
                                };    # if check_inlim(ysc)
                            };    # show_percent_box
                        };    # check_nuclim aln
                    };    # for $d
                };    # for $c
            };    # for $b
        };    # for $a
        print STDOUT "endseq\n";
    };    # $go
    $go || &report('NOAPLOTDATA');
}    # ps_plot_aln

sub ps_block_percent() {
    my $vlr = \%{ $Vars{LAYOUT} };
    &report( 'BLOCKTODO', 'PERCENT-BOX' );
    print STDOUT << "+++MAINProcs+++";
%
%%%%%%%% MATCHES PERCENT - BOX
%
beginmatches
%
+++MAINProcs+++
    $vlr->{_plot_x_sequence} && do {
        &ps_plot_ribbons( $vlr->{_x_sequence_name},  $T );
        &ps_plot_outlines( $vlr->{_x_sequence_name}, $T );
    };    # show_ribbons && _plot_x_sequence
    ( $vlr->{show_percent_box_grid} && !$vlr->{show_grid} ) && do {
        print STDOUT "\%\npvtick\n";
    };    # show_percent_box_grid
    $vlr->{_plot_aln} && do {
        &ps_plot_pctaln( $vlr->{_alignment_name} );
    };    # _plot_aln
    print STDOUT "\%\nendmatches\n\%\n";
    &report( 'BLOCKDONE', 'PERCENT-BOX' );
}    # ps_block_percent

sub ps_plot_pctaln() {
    my ($aln) = @_;
    &report( 'BLOCKA', $aln );
    print STDOUT "\%\n\%\%\%\% DATA ---> Alignment\n\%\nbeginaln\n";
    print STDOUT $StringCache{$aln};
    print STDOUT "endseq\n";
}    # ps_plot_pctaln

sub ps_block_extra() {
    my $vlr = \%{ $Vars{LAYOUT} };
    &report( 'BLOCKTODO', 'EXTRA-BOX' );
    print STDOUT << "+++MAINProcs+++";
%
%%%%%%%% EXTRA DATA - BOX
%
% numlines <- vlr->{_x_sequence_sources} % NGROUPS 
beginextra
%
+++MAINProcs+++
    $vlr->{_plot_x_sequence} && do {
        &ps_plot_ribbons( $vlr->{_x_sequence_name},  $T );
        &ps_plot_outlines( $vlr->{_x_sequence_name}, $T );
    };    # show_ribbons && _plot_x_sequence

    print STDOUT "\%\nendextra\n\%\n";
    &report( 'BLOCKDONE', 'EXTRA-BOX' );
}    # ps_block_extra

sub ps_plot_ribbons() {
    my ( $seq, $flg ) = @_;
    defined( $StringCache{$seq}{rbn} ) && do {
        my $what = $flg ? 'X-sequence' : 'Y-sequence';
        &report( 'BLOCKRB', $what );
        print STDOUT "\%\n\%\%\%\% $what Ribbons\n\%\n"
          . ( $flg ? 'beginXseq' : 'beginYseq' ) . "\n";
        print STDOUT $StringCache{$seq}{rbn};
        print STDOUT "endseq\n";
    };    # defined ribbons
}    # ps_plot_ribbons

sub ps_plot_outlines() {
    my ( $seq, $flg ) = @_;
    ( defined( $StringCache{$seq}{lns} )
      && scalar( @{ $StringCache{$seq}{lns} } ) > 0 )
      && do {
        my $what = $flg ? 'X-sequence' : 'Y-sequence';
        &report( 'BLOCKRBL', $what );
        print STDOUT "\%\n\%\%\%\% Finishing $what Ribbons\n\%\n"
          . ( $flg ? 'beginXseq' : 'beginYseq' ) . "\n";
        $n = 0;
        foreach my $pos ( @{ $StringCache{$seq}{lns} } ) {
            ( $n % 10 == 0 ) && print STDOUT '[ ';
            print STDOUT "$pos ";
            $n++;
            ( $n % 10 == 0 ) && print STDOUT "Rln\n";
        };    # foreach $pos
        ( $n % 10 != 0 ) && print STDOUT "Rln\n";
        print STDOUT "endseq\n";
    };    # defined ribbons lines
}    # ps_plot_outlines

sub get_ordered_ary() {
    my ( $main, $name ) = @_;
    my ( $k, $h );
    for ( $k = 0 ; $k <= $#{ $Order{$main} } ; $k++ ) {
        $h = $Order{$main}[$k][0];
        $h eq $name && do {
            return ( $T, \@{ $Order{$main}[$k][1] } );
        };
    }
    return ( $F, undef );
}    # get_ordered_ary

sub show_element() {
    my $rf = shift;
    ( defined( $rf->[$_prop]{hide} ) && $rf->[$_prop]{hide} ) && return $F;
    return $T;
}    # show_element

sub check_seqlim() {
    my ( $o, $e, $O, $E ) = @_;
    ( $o < $E && $e > $O ) && return $T;
    return $F;
}    # check_seqlim

sub check_inlim() {
    my ( $p, $l, $u ) = @_;
    ( $p >= $l && $p <= $u ) && return $T;
    return $F;
}    # check_inlim

sub get_var_value() {
    my ( $name, $refs ) = @_;
    my ( $out,  $v );

    # print STDERR ("#"x40)."\n";
    foreach $v (qw/ FEATURE GROUP STRAND SOURCE SEQUENCE /) {
        defined( $Defaults{$v}{$name} ) && ( $out = $Defaults{$v}{$name} );

        # print STDERR "### GET VAR VALUE: $v ($name) : ".
#              (defined($Defaults{$v}{$name})?$Defaults{$v}{$name}:"undef")."\n";
    };    # foreach
    foreach $v (@$refs) {
        defined( $v->[$_prop]{$name} ) && ( $out = $v->[$_prop]{$name} );

        # print STDERR "### GET VAR VALUE: $v ($name) : ".
#              (defined($v->[$_prop]{$name})?$v->[$_prop]{$name}:"undef")."\n";
    };    # foreach
    ( defined( $Vars{LAYOUT}{$name} )
      && ( !defined( $Defaults{LAYOUT}{$name} )
          || $Vars{LAYOUT}{$name} ne $Defaults{LAYOUT}{$name} ) )
      && ( $out = $Vars{LAYOUT}{$name} );

    # print STDERR "### GET VAR VALUE: Vars{LAYOUT} ($name) : ".
#              (defined($Vars{LAYOUT}{$name})?$Vars{LAYOUT}{$name}:"undef")."\n".
    #              "### GET VAR VALUE: Defaults{LAYOUT} ($name) : ".
#              (defined($Defaults{LAYOUT}{$name})?$Defaults{LAYOUT}{$name}:"undef")."\n";
# print STDERR "### GET VAR VALUE: $name = ".(defined($out)?$out:"undef")."\n";
    return $out;
}    # get_var_value

sub norm_score() {
    my ( $score, $minsco, $maxsco ) = @_;
    $minsco == $maxsco && return 1;
    $score < $minsco && ( $score = $minsco );
    $score > $maxsco && ( $score = $maxsco );
    return sprintf( "%5.3f", ( ( $score - $minsco ) / ( $maxsco - $minsco ) ) );
}    # norm_score

sub ps_header() {
    my $vr = \%{ $Vars{LAYOUT} };
    print STDOUT << "+++HEADER+++";
%!PS-Adobe-3.0
%%Title: $vr->{title}
%%Creator: $PROGRAM
%%Version: $VERSION
%%CreationDate: $DATE
%%For: $USER
%%Pages: 1
%%Orientation: Portrait
%%BoundingBox: 0 0 $vr->{_page_width} $vr->{_page_height}
%%EndComments
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                          GFF2APLOT                               %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%    Converting alignments in GFF format to PostScript dotplots.
% 
%     Copyright (C) 1999 - Josep Francesc ABRIL FERRANDO  
%                                  Thomas WIEHE                   
%                                 Roderic GUIGO SERRA       
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $LAST_UPDATE
%
% Report BUGS to: jabril\@imim.es 
%
%%BeginProlog
%
%%BeginProcSet: Short_names 1.0 0
%
/B { bind def } bind def /X { exch def } B /D { def } B /_t { true } B /_f { false } B /S { gsave } B /R { grestore } B /F { scale } B /T { translate } B /m { moveto } B /rm { rmoveto } B /l { lineto } B /rl { rlineto } B /K { stroke } B /ie { ifelse } B /cmyk { setcmykcolor } B /slw { setlinewidth } B /solid { [ ] 0 setdash } D /dotted { [ 1 ] 0 setdash } D /ldotted { [ 1 2 ] 0 setdash } D /dashed { [ 3 3 ] 0 setdash } D /ldashed { [ 6 4 ] 0 setdash } D /ddashed { [ 4 3 1 3 ] 0 setdash } D
%
%%EndProcSet:   Short_names 1.0 0
%
%%BeginProcSet: Constants 1.0 0
%
/pt { } B /cm { 28.35 mul } B /icm { 28.35 div } B /in { 72 mul } B
/OST 0.25 cm D /OSB 0.25 cm D /OSL 0.25 cm D /OSR 0.25 cm D
%
%%EndProcSet:   Constants 1.0 0
%
%%BeginProcSet: Global_functions 1.0 0
%
/vflg _f D /tflg @{[ &tobool($Debug) ]} D tflg { mark } if
/msg { print (\\n) print flush } B /mst { print counttomark dup { dup index 20 string cvs print ( ) print 1 sub } repeat (\\n) print pop flush } B /msc { mst counttomark 1 add { pop } repeat } B
/bbox { 4 copy 3 1 roll exch 6 2 roll 8 -2 roll m l l l closepath } B
/obox { 2 div dup dup dup 7 1 roll 5 1 roll 3 1 roll add 7 1 roll add 6 1 roll sub 5 1 roll sub 4 1 roll bbox } B /ibox { 2 div dup dup dup 7 1 roll 5 1 roll 3 1 roll sub 7 1 roll sub 6 1 roll add 5 1 roll add 4 1 roll bbox } B
/min { 2 copy gt { exch } if pop } B /max { 2 copy lt { exch } if pop } B
/clim { min max } B /zlim { dup 0 eq { 3 { pop } repeat 0 } { clim } ie } B
/hip { dup mul exch dup mul add sqrt } B
%
%%EndProcSet:   Global_functions 1.0 0
%
%%BeginProcSet: Text_functions 1.0 0
%
/sfnt { findfont exch scalefont setfont } B
/glblw { sfnt (Z) stringwidth pop mul } B
/chrh { S newpath 0 0 m _f charpath flattenpath pathbbox exch pop 3 -1 roll pop R } B
/strh { 2 dict begin /lly 0.0 D /ury 0.0 D { ( ) dup 0 4 -1 roll put chrh dup ury gt { /ury X } { pop } ie dup lly lt { /lly X } { pop } ie } forall ury end } B
/ttshw { vflg { 9 copy 8 { pop } repeat mark exch dup length 20 gt { 0 20 getinterval (...) } if (ttshw -> ) msc } if S cmyk sfnt 8 dict begin /h X /v X /lbl X /angle X /y X /x X /hs lbl stringwidth pop D /vs lbl strh D x y T angle rotate h (rh) eq { hs } { h (ch) eq { hs 2 div } { 0 } ie } ie neg v (tv) eq { vs } { v (cv) eq { vs 2 div } { 0 } ie } ie neg m lbl show end R } B
%
%%EndProcSet:   Text_functions 1.0 0
%
%%BeginProcSet: Pseudohash_dicts 1.0 0
%
+++HEADER+++
    &report('PSHEADER');
}    # ps_header

sub ps_colors() {
    my ( $colornum, $j, @tmp, @colary, @gradary );
    @tmp     = ();
    @colary  = ();
    @gradary = ();
    @colary  = ( keys %{ $USED{'COLORS'} } );
    ($Debug) && ( @colary = ( keys %COLORS ) );
    @gradary = ( keys %{ $USED{'GRADIENTS'} } );
    ($Debug) && ( @gradary = ( keys %GRADIENTS ) );
    $colornum = scalar(@colary) + scalar(@gradary);

    # &min(scalar(@colary), $colors);
    print STDOUT "%% Fixed Color Variables (CMYK)\n";
    print STDOUT "/colordict "
      . ( $colornum + 28 )
      . " dict D colordict begin %% using "
      . $colornum
      . " colors + 28 definitions\n";
    @tmp =
      map { $_->[1] }
      sort { $a->[0] <=> $b->[0] } map { [ $COLORS{$_}->[0], $_ ] } @colary;

    # print STDERR "### colary: @colary\n### colary-tmp: @tmp\n";
    for ( $j = 0 ; $j <= $#tmp ; $j++ ) {
        my $name = $tmp[$j];
        my $ref  = \$COLORS{$name};
        my $cmyk = "$$ref->[1] $$ref->[2] $$ref->[3] $$ref->[4]";
        print STDOUT "/" . $name . " { $cmyk } D"; # (&fill_right($name,20," "))
        print STDOUT ( ( $j % 5 == 0 && $j > 0 ) ? "\n" : " " );
    };    # for $j
    print STDOUT "\n" if ( $j % 5 > 0 );
    &ps_color_gradients( \@gradary ) if scalar(@gradary) > 0;
    print STDOUT "end %% colordict\n";
    &ps_color_pattern;
    &report('PSCOLORS');
}    # ps_colors

sub ps_color_gradients() {
    my $gary = shift;
    foreach my $gr (@$gary) {
        print STDOUT "$GRADIENTS{$gr}\n";
    }
}    # ps_color_gradients

sub ps_color_pattern() {
    print STDOUT "%% Ribbon color mode (cmyk/pattern)\n";
    $Vars{LAYOUT}{ribbon_color_merge} || do {
        print STDOUT '/rbxcol { cmyk } B' . "\n";
        return;
    };
    $Vars{LAYOUT}{color_merge_factor} =
      &max( 0.07, &min( 2, $Vars{LAYOUT}{color_merge_factor} ) );
    print STDOUT '/PG '
      . $Vars{LAYOUT}{color_merge_factor}
      . ' D % "alpha" color blending factor ' . "\n";
    print STDOUT << '+++PatternProcs+++';
/PS PG 2 mul D
<< /PaintType 2 /PatternType 1 /TilingType 1 /BBox [ PG neg PG neg PG PG ] /XStep PS /YStep PS /PaintProc { pop PG slw 0 setlinecap PG 2 div 0 m 0 PG rl K PG 2 div neg 0 m 0 PG neg rl K } >> matrix makepattern /HP X
<< /PaintType 2 /PatternType 1 /TilingType 1 /BBox [ PG neg PG neg PG PG ] /XStep PS /YStep PS /PaintProc { pop PG slw 0 setlinecap PG 2 div neg 0 m 0 PG rl K PG 2 div 0 m 0 PG neg rl K } >> matrix makepattern /VP X
[ /Pattern /DeviceCMYK ] setcolorspace
/rbxcol { H { HP } { VP } ie setcolor } B
+++PatternProcs+++
}    # ps_color_pattern

sub ps_page_formats() {
    my ( $formatnum, $j, @tmp, @fmtary );
    @tmp    = ();
    @fmtary = ();
    @fmtary = ( keys %{ $USED{FORMATS} } );
    ($Debug) && ( @fmtary = ( keys %FORMATS ) );
    $formatnum = scalar(@fmtary);    # &min(scalar(@fmtary), $formats);
    print STDOUT "%% Paper Sizes (in points)\n";
    print STDOUT "/pagedict "
      . ( $formatnum + 2 )
      . " dict D pagedict begin %% "
      . $formatnum
      . " formats + 2 definitions\n";
    @tmp =
      map { $_->[1] }
      sort { $a->[0] <=> $b->[0] } map { [ $FORMATS{$_}->[0], $_ ] } @fmtary;

    for ( $j = 0 ; $j <= $#tmp ; $j++ ) {
        my $name = $tmp[$j];
        my $ref  = \$FORMATS{$name};
        my $pgsz = "$$ref->[1] $$ref->[2]";

        # &fill_left($$ref->[1],4," ").&fill_left($$ref->[2],5," ");
        print STDOUT "/pg" . $name
          . "{ $pgsz } D";    # (&fill_right($name,10," "))
        print STDOUT ( ( $j % 5 == 0 && $j > 0 ) ? "\n" : " " );
    }
    print STDOUT "\n" if ( $j % 5 > 0 );
    print STDOUT "end %% pagedict\n";
    &report('PSFORMATS');
}    # ps_page_formats

sub ps_variables() {
    my $vr = \%{ $Vars{LAYOUT} };
    print STDOUT << "+++PSVARS+++";
%
%%EndProcSet:   Pseudohash_dicts 1.0 0
%
tflg { (%%% Basic settings and functions were defined %%%) msg } if
%
%%BeginProcSet: Setting_Vars 1.0 0
%
tflg { (%%% Setting variables %%%) msg } if
/H _t D /HT { /H _t D } B /VT { /H _f D } B /flgcrd @{[ &tobool(! $vr->{hide_credits}) ]} D /flglscape @{[ &tobool($F) ]} D /Pb @{[ &tobool($vr->{show_percent_box}) ]} D /Qb @{[ &tobool($vr->{show_extra_box}) ]} D /axesb @{[ &tobool($vr->{aplot_xy_same_length}) ]} D /XYR $vr->{aplot_xy_scale} D /ZM @{[ &tobool($vr->{zoom_marks}) ]} D
/Dpage { pagedict begin pg$vr->{page_size} flglscape { exch } if end } D /PPC pagedict begin pga4 end hip Dpage hip exch div D /Ch 4 pt 15 pt 4 pt PPC mul clim D
/MT @{ $vr->{margin_top} } D /MB @{ $vr->{margin_bottom} } D /ML @{ $vr->{margin_left} } D /MR @{ $vr->{margin_right} } D
/bg { colordict begin @{[ ($vr->{background_color} ne 'bg') ? $vr->{background_color} : 'white' ]} end } D /fg { colordict begin @{[ ($vr->{foreground_color} ne 'fg') ? $vr->{foreground_color} : 'black' ]} end } D /ABc { colordict begin $vr->{aplot_box_bgcolor} end } D /PBc { colordict begin $vr->{percent_box_bgcolor} end } D /QBc { colordict begin $vr->{extra_box_bgcolor} end } D
/XO $vr->{_x_start} D /XE $vr->{_x_end} D /XD XE XO sub 1 add D /YO $vr->{_y_start} D /YE $vr->{_y_end} D /YD YE YO sub 1 add D /PO $vr->{_pbox_min_score} D /PE $vr->{_pbox_max_score} D /PD PE PO sub D
/QO 0 D /QE 1 D /QD PE PO sub D
/Ph @{ $vr->{percent_box_height} } D /Qh @{ $vr->{extra_box_height} } D
/KGb @{[ &tobool($vr->{show_grid}) ]} D /KLb @{[ &tobool($vr->{show_tickmark_label}) ]} D /KAXb @{[ &tobool($vr->{show_aplot_x_ticks}) ]} D /KAYb @{[ &tobool($vr->{show_aplot_y_ticks}) ]} D /KPXb @{[ &tobool($vr->{show_percent_x_ticks}) ]} D /KPYb @{[ &tobool($vr->{show_percent_y_ticks}) ]} D /KQXb @{[ &tobool($vr->{show_extrabox_x_ticks}) ]} D /KQYb @{[ &tobool($vr->{show_extrabox_y_ticks}) ]} D
/KASm $vr->{major_tickmark_nucleotide} D /KASn $vr->{minor_tickmark_nucleotide} D /KPSm $vr->{major_tickmark_score} D /KPSn $vr->{minor_tickmark_score} D /KQSm 10 D /KQSn 10 D
/TTLb @{[ &tobool($vr->{show_title}) ]} D /TTLl @{[ &tostring($vr->{title}) ]} D /TTLz @{ $vr->{title_fontsize} } D /TTLf { TTLz $Fonts{$vr->{title_font}} } D /TTLp 10 pt 25 pt TTLz 2 div zlim D /TTlb @{[ &tobool($vr->{show_subtitle}) ]} D /TTll @{[ &tostring($vr->{subtitle}) ]} D /TTlz @{ $vr->{subtitle_fontsize} } D /TTlf { TTlz $Fonts{$vr->{subtitle_font}} } D /TTlp 10 pt 20 pt TTlz zlim D
/XSLb @{[ &tobool($vr->{show_x_label}) ]} D /XSLl @{[ &tostring($vr->{x_label}) ]} D /XSLz @{ $vr->{x_label_fontsize} } D /XSLf { XSLz $Fonts{$vr->{x_label_font}} } D /XSLp 2.5 pt 10 pt XSLz 0.25 mul zlim D /YSLb @{[ &tobool($vr->{show_y_label}) ]} D /YSLl @{[ &tostring($vr->{y_label}) ]} D /YSLz @{ $vr->{y_label_fontsize} } D /YSLf { YSLz $Fonts{$vr->{y_label_font}} } D /YSLp 2.5 pt 10 pt YSLz 0.25 mul zlim D
/PBLb @{[ &tobool($vr->{show_percent_box_sublabel}) ]} D /PBLl @{[ &tostring($vr->{percent_box_label}) ]} D /PBLz @{ $vr->{percent_box_label_fontsize} } D /PBLf { PBLz $Fonts{$vr->{percent_box_label_font}} } D /PBLp 10 pt 20 pt PBLz zlim D /PBlb @{[ &tobool($vr->{show_percent_box_sublabel}) ]} D /PBll @{[ &tostring($vr->{percent_box_sublabel}) ]} D /PBlz @{ $vr->{percent_box_sublabel_fontsize} } D /PBlf { PBlz $Fonts{$vr->{percent_box_sublabel_font}} } D /PBlp 10 pt 20 pt PBlz zlim D
/QBLb @{[ &tobool($vr->{show_extra_box_label}) ]} D /QBLl @{[ &tostring($vr->{extra_box_label}) ]} D /QBLz @{ $vr->{extra_box_label_fontsize} } D /QBLf { QBLz $Fonts{$vr->{extra_box_label_font}} } D /QBLp 10 pt 20 pt QBLz zlim D /QBlb @{[ &tobool($vr->{show_extra_box_sublabel}) ]} D /QBll @{[ &tostring($vr->{extra_box_sublabel}) ]} D /QBlz @{ $vr->{extra_box_sublabel_fontsize} } D /QBlf { QBlz $Fonts{$vr->{extra_box_sublabel_font}} } D /QBlp 10 pt 20 pt QBlz zlim D
/FTz @{ $vr->{feature_label_fontsize} } D /FTf { FTz $Fonts{$vr->{feature_label_font}} } D /FTp 2 pt 10 pt FTz 0.5 mul zlim D /FTXw FTz 0 gt { $vr->{feature_x_label_length} FTf glblw } { 0 } ie D /FTXa $vr->{feature_x_label_angle} D /FTXrb @{[ &tobool($vr->{feature_x_label_rotate}) ]} D /FTYw FTz 0 gt { $vr->{feature_y_label_length} FTf glblw } { 0 } ie D /FTYa $vr->{feature_y_label_angle} D /FTYrb @{[ &tobool($vr->{feature_y_label_rotate}) ]} D
/GPz @{ $vr->{group_label_fontsize} } D /GPf { GPz $Fonts{$vr->{group_label_font}} } D /GPp 4 pt 16 pt GPz 0.75 mul zlim D /GPpt GPp 0.25 mul D /GPXw GPz 0 gt { $vr->{group_x_label_length} GPf glblw } { 0 } ie D /GPXa $vr->{group_x_label_angle} D /GPXrb @{[ &tobool($vr->{group_x_label_rotate}) ]} D /GPYw GPz 0 gt { $vr->{group_y_label_length} GPf glblw } { 0 } ie D /GPYa $vr->{group_y_label_angle} D /GPYrb @{[ &tobool($vr->{group_y_label_rotate}) ]} D
/DRh 10 pt D /DRhb DRh 0.1 mul D /DRht DRh DRhb sub D /DRFTlh DRh FTp add D /DRGPlh DRFTlh FTz add GPp add D /DRGPlnh DRGPlh GPpt 2 mul sub D
/Blw 2 pt D /Alw 1 pt D /Clw 0.1 pt D /Clw2 Clw 2 mul D /BClw Blw Clw sub D
/KLz 8 pt D /KLf { KLz /Helvetica } D /KWx 1.5 pt D /KWn KWx 0.1 mul D /KHx KLz 0.75 mul D /KHn KLz 0.50 mul D /KHt KLz 1.0 mul D /KLp KLz 0.5 mul D /KLh KHt KLz add KLp add D /KLw @{[ &max(length($vr->{_x_end}), length($vr->{_y_end}), length($vr->{_pbox_max_score})) ]} KLf glblw D
%
%%EndProcSet:   Setting_Vars 1.0 0
%
+++PSVARS+++
    &report('PSVARS');
}    # ps_variables

sub ps_main() {
    print STDOUT << '+++MAINProcs+++';
%%BeginProcSet: Page_layout 1.0 0
%
tflg { (%%% Computing page layout %%%) msg } if
/MT MT OST max D /MB MB OSB max D /ML ML OSL max D /MR MR OSR max D
/FL { Dpage MT MB add flgcrd { Ch add } if sub exch ML MR add sub exch } D /FXO ML D /FYO Dpage exch pop MT sub D
/getXt { dup sin 0 eq { pop (ch) } { 90 sub cos 0 ge { (lh) } { (rh) } ie } ie } B /getYt { dup cos 0 eq { pop (cv) } { 90 sub sin 0 ge { (tv) } { (bv) } ie } ie } B
/FTXh FTz FTXa cos mul FTXw FTXa sin abs mul add FTp add D /FTXt { FTXa getYt FTXa getXt } D /FTYh FTz FTYa cos mul FTYw FTYa sin abs mul add FTp add D /FTYt { FTYa getYt FTYa getXt } D
/GPXh GPz GPXa cos mul GPXw GPXa sin abs mul add GPp add D /GPXt { GPXa getYt GPXa getXt } D /GPYh GPz GPYa cos mul GPYw GPYa sin abs mul add GPp add D /GPYt { GPYa getYt GPYa getXt } D
/ATX GPXh GPp add FTXh FTp add add DRh add Blw add D /ATY GPYh GPp add FTYh FTp add add DRh add Blw add D
/AXL YSLz YSLp add ATY add D /AXR FL pop KAYb KPYb or KQYb or { KLw KHt add sub } if D /AX AXR AXL sub D 
/AYT TTLz TTLp add TTlz TTlp add add XSLz XSLp add add ATX add D /AYB AYT AX add D /AY AX D
XYR 0 le { /XYR 1 D } if axesb not { XYR 1 ge { /AY AY XYR div D } { /AX AX XYR mul D } ie } if
/PYT AYB KAXb { KLb { KLh } { KHt } ie } { KLp } ie add D
/rspc FL exch pop PYT sub Pb { KPXb { KLb { KLh } { KHt } ie } { KLp } ie sub } if Qb { KQXb { KLb { KLh KLp sub } { KHt } ie } if sub } if D
/PY Pb { Ph 0 gt { 1.5 cm rspc Qb { 2 div } if Ph clim } { 2.0 cm } ie } { 0 } ie D /PYB PYT PY add D
/QYT PYB Pb { KPXb { KLb { KLh } { KHt } ie } { KLp } ie add } if D /QY Qb { rspc Pb { PY sub } if Qh 0 gt { 2 cm exch Qh clim } if } { 0 } ie D /QYB QYT QY add D
%
%%EndProcSet:   Page_Layout 1.0 0
%
%
%%BeginProcSet: Aplot_dict 1.0 0
%
tflg { (%%% Setting aplot dictionary %%%) msg } if
/aplotdict 120 dict D aplotdict begin
/Xscm { Xscale mul } B /Xscme { Xscm exch } B
/Yscm { Yscale mul } B /Yscme { Yscm exch } B
/Hscm { H { Xscm } { Yscm } ie } B
/rF { 1 exch dup 0 eq { pop 10 -9 exp } if div } B 
/Fm { Yscme Xscme m } B /Fl { Yscme Xscme l } B
/corner { S 4 -2 roll m 0 0 4 1 roll 2 { S rl fg cmyk Blw slw [Blw 3] 0 setdash K R } repeat R } B 
/zoomtk { 4 copy 3 1 roll exch ATX neg ATY corner ATX ATY neg corner ATX ATY corner ATX neg ATY neg corner } D
/gradcol { dup dup dup 7 1 roll 5 1 roll 3 1 roll mul 7 1 roll mul 6 1 roll mul 5 1 roll mul 4 1 roll } B
/shp { S 7 dict begin /s_s X Hscm /s_e X Hscm /s_o X /s_p X /f_x { s_e s_o sub } D /s_m { f_x 2 div s_o add } D 0 BClw T { H { FTXa } { FTYa } ie s_m DRFTlh 4 -2 roll exch H { FTXt } { FTYt } ie FTf fg ttshw } if S s_o 0 T { DRht s_p mul DRhb add } { DRh } ie /f_y X f_x f_y F s_s (-) eq { 1 0 T -1 1 F } if { s_p gradcol } if cmyk cvx exec S fill R f_x rF f_y rF F fg cmyk Clw slw 2 setmiterlimit K R end R } B
/rbx { S rbxcol H { YH neg } { AX neg } ie 3 1 roll Hscm 3 1 roll Hscm 0 Clw2 ibox fill R } B /rln { S fg cmyk Hscm H { YH neg } { AX neg } ie 0 3 -1 roll dup 4 1 roll exch m l Clw slw K R } B /Rln { counttomark { rln } repeat pop } B
/grp { S 5 dict begin /s_s X Hscm /s_e X Hscm /s_o X /f_x { s_e s_o sub } D /s_m { f_x 2 div s_o add } D 0 BClw T { H { GPXa } { GPYa } ie s_m DRGPlh 4 -2 roll exch H { GPXt } { GPYt } ie GPf fg ttshw } if S cmyk { S s_o DRGPlnh s_o 0 s_e DRGPlnh s_e 0 2 { m l dotted Clw slw K } repeat R } if s_s (-) eq { s_e DRGPlnh T -1 1 F } { s_o DRGPlnh T } ie f_x exch Clw2 slw 2 setmiterlimit cvx exec R end R } B
/aln { S dup 4 1 roll exch { Alw mul Alw add slw } { pop Alw slw } ie { gradcol } { pop } ie cmyk Fm Fl K R } B
/xbox { S Yscme Yscme 4 -2 roll Xscme Xscm 4 1 roll exch bbox { S cmyk fill R } if cvx exec slw cmyk K R } B
/xcir { S Yscme Xscme T rotate Yscme Xscme 2 copy 5 2 roll F 0 0 1 0 360 arc closepath { S 6 2 roll cmyk fill R } if 2 { 1 exch div exch } repeat F cvx exec slw cmyk K R } B
/xtxt { S Yscme Xscme counttomark 2 roll not { fg } if ttshw pop R } B
/xarw { S Yscme Xscme T dup neg counttomark 1 roll rotate 2 dict begin /la X /lb X cmyk 0 0 m lb 6 mul dup la gt { pop la 0.75 mul } if dup dup dup 3 div dup neg 3 1 roll l l closepath fill 0 m la 0 l lb slw K la lb add end 0 T 0 0 xtxt R } B
/chktk { 3 dict begin /pk X /ck 0 D /ek 1 D { ck 5 gt { exit } if /ek 10 ck exp D pk ek mul cvi ek div pk eq { exit } if /ck ck 1 add D } loop ck end } B /mymod { 2 copy chktk exch chktk max dup 0 eq { pop 2 { cvi exch } repeat mod } { 1 dict begin /ff 10 3 -1 roll exp D 2 { ff mul round cvi exch } repeat mod ff div end } ie } B /rndtk { dup abs 10 -6 exp le { pop 0 } if } B
/mkgrid { S 0 0 m rl fg cmyk { dashed KWn 2 mul slw } { ldotted KWn slw } ie K R } B /mjrtick { S 0 0 m { 0 KHx neg } { KHx 0 } ie rl KWx slw fg cmyk K R } B /mnrtick { S 0 0 m { 0 KHn neg } { KHn 0 } ie rl KWn slw fg cmyk K R } B
/hticklbl { S 10 string cvs 0 KHt neg 0 4 -1 roll (tv) (ch) KLf fg ttshw R } B /htickslbl { S 10 string cvs (tv) 3 -1 roll ZM { { (lh) KHn } { (rh) KHn neg } ie } { pop (ch) 0 } ie KHt neg 0 6 -3 roll KLf fg ttshw R } B
/htick { S 10 dict begin /tmn X /tmx X /te X /to X /yt X /lxo to to tmx mymod sub D /lxe tmx te tmx mymod sub te add D /lno to to tmn mymod sub D /lne tmn te tmn mymod sub te add D S /cnt lxo D to neg Xscm 0 T { cnt te ge { exit } if cnt to gt cnt te lt and { cnt Xscm S 0 T _t mjrtick KGb { _t 0 yt mkgrid } if cnt rndtk hticklbl R } if /cnt cnt tmx add D } loop R S /cnt lno D to neg Xscm 0 T { cnt te ge { exit } if cnt to gt cnt te lt and cnt tmx mymod 0 eq not and { cnt Xscm S 0 T _t mnrtick KGb { _f 0 yt mkgrid } if R } if /cnt cnt tmn add D } loop S to Xscm 0 T _t mjrtick _f to htickslbl R S te Xscm 0 T _t mjrtick _t te htickslbl R R end R } B
/vticklbl { S 10 string cvs KHt 0 0 4 -1 roll (cv) (lh) KLf fg ttshw R } B /vtickslbl { S 10 string cvs KHt ZM { KHn } { 0 } ie 0 5 -1 roll ZM { { (bv) } { (bv) } ie } { pop (cv) } ie 5 -1 roll exch (lh) KLf fg ttshw R } B
/vtick { S 10 dict begin /tmn X /tmx X /te X /to X /xt X /lxo to to tmx mymod sub D /lxe tmx te tmx mymod sub te add D /lno to to tmn mymod sub D /lne tmn te tmn mymod sub te add D xt 0 T S /cnt lxo D 0 to neg Yscm T { cnt te ge { exit } if cnt to gt cnt te lt and { cnt Yscm S 0 exch T _f mjrtick KGb { _t xt neg 0 mkgrid } if cnt rndtk vticklbl R } if /cnt cnt tmx add D } loop /cnt lno D { cnt te ge { exit } if cnt to gt cnt te lt and cnt tmx mymod 0 eq not and { cnt Yscm S 0 exch T _f mnrtick KGb { _f xt neg 0 mkgrid } if R } if /cnt cnt tmn add D } loop R S 0 to neg Yscm T S 0 to Yscm T _f mjrtick _f to vtickslbl R S 0 te Yscm T _f mjrtick _t te vtickslbl R R end R } B
/pvtick { S 6 dict begin /tmx KPSm D /te PE D /to PO D /xt AX D /lxo to to tmx mymod sub D xt 0 T S /cnt lxo D 0 to neg Yscm T { cnt te ge { exit } if cnt to gt cnt te lt and { cnt Yscm S 0 exch T _t xt neg 0 mkgrid R } if /cnt cnt tmx add D } loop R end R } B
/nucltick { XO XE KASm KASn htick } B
/mkshp { m { l } repeat closepath } B
/fbox { 0.0 1.0 1.0 1.0 1.0 0.0 3 0.0 0.0 mkshp } B /hbox { 0.0 0.5 1.0 0.5 1.0 0.0 3 0.0 0.0 mkshp } B
/fraw { 0.75 0.0 0.75 -0.1 1.00 0.5 0.75 1.1 0.75 1.0 0.00 1.0 6 0.0 0.0 mkshp } B /hraw { 1.00 0.0 0.75 1.1 0.75 1.0 0.00 1.0 4 0.0 0.0 mkshp } B
/flaw { 1 0 T -1 1 F fraw } B /hlaw { 1 0 T -1 1 F hraw } B
/frae { 0.25 0.5 0.00 1.0 1.00 1.0 1.00 0.0 4 0.0 0.0 mkshp } B /hrae { 0.00 1.0 1.00 1.0 1.00 0.0 3 0.25 0.0 mkshp } B
/flae { 1 0 T -1 1 F frae } B /hlae { 1 0 T -1 1 F hrae } B
/frsg { 0.75 0.0 0.75 -0.1 1.00 0.5 0.75 1.1 0.75 1.0 0.00 1.0 0.25 0.5 7 0.00 0.0 mkshp } B
/hrsg { 1.00 0.0 0.75 1.1 0.75 1.0 0.00 1.0 4 0.25 0.0 mkshp } B
/flsg { 1 0 T -1 1 F frsg } B /hlsg { 1 0 T -1 1 F hrsg } B
/frtn { 0.0 1.0 1.0 0.5 2 0.0 0.0 mkshp } B /hrtn { 0.0 1.0 1.0 0.0 2 0.0 0.0 mkshp } B
/fltn { 1 0 T -1 1 F frtn } B /hltn { 1 0 T -1 1 F hrtn } B
/futn { 0.5 1.0 1.0 0.0 2 0.0 0.0 mkshp } B /fdtn { 0.0 1.0 1.0 1.0 2 0.5 0.0 mkshp } B
/fdmd { 0.0 0.5 0.5 1.0 1.0 0.5 3 0.5 0.0 mkshp } B /fmdm { 1.0 0.0 0.0 1.0 1.0 1.0 3 0.0 0.0 mkshp } B
/fcir { 0.5 0.5 0.5 0 360 arc closepath } B /farc { 1 2 F 0.5 0 0.5 0 180 arc closepath } B /harc { 0.5 0 0.5 0 180 arc closepath } B
/gmid { 2 div dup 0 m GPpt l } B /gbas { S { dotted } if 0 0 m 0 l K R } B
/glns { dup _f gbas gmid K } B /glnd { dup _t gbas gmid K } B /gbrk { dup dup _f gbas gmid dup 0 m GPpt neg l 0 0 m 0 GPpt neg l K } B /gbrc { 4 dict begin /_l X /_lh _l 2 div D /_a GPpt 2 mul D /_b GPpt 4 mul D 0 GPpt neg m 0 _a sub _b _lh _b neg _lh GPpt curveto _l GPpt neg m _l _a add _b _lh _b neg _lh GPpt curveto K end } B
/gseg { exch dup 3 -1 roll gbas dup gmid dup GPpt m GPpt neg l 0 GPpt m 0 GPpt neg l K } B /gsgm { _f gseg } B /gsgd { _t gseg } B
/garw { 3 -1 roll dup 3 -1 roll gbas dup gmid 0 GPpt neg m 0 GPpt l K 0 T GPpt GPpt F { 0 0 m -2.5 1 l -2.5 -1 l closepath S fill R } { 0 0 m -2.5 1 l 0 0 m -2.5 -1 l } ie K } B /garh { 3 -1 roll dup 3 -1 roll gbas dup gmid 0 0 m 0 GPpt l K 0 T GPpt GPpt F { 0 0 m -2.5 1 l -2.5 0 l closepath S fill R } { 0 0 m -2.5 1 l } ie K } B
/gsaw { _f _f garw } B /gsad { _f _t garw } B /gfaw { _t _f garw } B /gfad { _t _t garw } B
/ghaw { _t _f garh } B /ghad { _t _t garh } B /ghsw { _f _f garh } B /ghsd { _f _t garh } B
%
/DoHeader { TTLb { S AXL TTLz neg 0 TTLl (bv) (lh) TTLf fg ttshw R } if TTlb { S AXL TTLz TTLp add TTlz add neg 0 TTll (bv) (lh) TTlf fg ttshw R } if } B
/shwcrd { S 2 dict begin /Cz Ch 1 pt PPC 2 div mul sub D /Cf { Cz /Courier } D FL Ch add neg 2 copy 0 (This plot has been obtained using GFF2APLOT. The most recent version is freely available at \042http:\057\057www1.imim.es\057software\057gfftools\057GFF2APLOT.html\042. Copyright\040\0401999-2003 by Josep F. ABRIL, Thomas WIEHE & Roderic GUIGO) (bv) (rh) Cf fg ttshw exch S Cf sfnt (\0401999-2003 by Josep F. ABRIL, Thomas WIEHE & Roderic GUIGO) stringwidth pop sub R exch 0 (\343) (bv) (ch) Cz /Symbol fg ttshw end R } B
/begindata { vflg { (### BEGINdata) msg } if S AXL AYB neg T S 0 0 AX AY 4 copy bbox ABc cmyk fill Blw obox fg cmyk Blw slw K R /Xscale AX XD div D /Yscale AY YD div D /YH AY D /YY YO D S KAXb { AY nucltick } if R S KAYb { AX YO YE KASm KASn vtick } if R S XSLb { AX 2 div AY ATX add XSLp add 0 XSLl (bv) (ch) XSLf fg ttshw } if YSLb { ATY YSLp add neg AY 2 div 90 YSLl (bv) (ch) YSLf fg ttshw } if R } B
/enddata { ZM { S 0 0 AX AY zoomtk R /ZM _f D } if R vflg { (### ENDdata) msg } if } B
/beginXseq { S 0 0 AX YH ATX add bbox clip newpath XO neg Xscm YH T colordict begin HT } B /beginYseq { S ATY neg 0 AX YH bbox clip newpath 0 YO neg Yscm T 90 rotate colordict begin VT } B /beginaln { S 0 Alw 2 div neg AX YH Alw add bbox clip newpath XO neg Xscm YY neg Yscm T colordict begin } B /endseq { end R } B
/beginmatches { vflg { (### BEGINmatches) msg } if S AXL PYB neg T S 0 0 AX PY 4 copy bbox PBc cmyk fill Blw obox fg cmyk Blw slw K R /Yscale PY PD div D /YH PY D /YY PO D S KPXb { PY nucltick } if R S KPYb { AX PO PE KPSm KPSn vtick } if R S PBLb { DRFTlh PBlz add PBLp add neg PY 2 div 90 PBLl (bv) (ch) PBLf fg ttshw } if PBLb { DRFTlh neg PY 2 div 90 PBll (bv) (ch) PBlf fg ttshw } if R S } D
/endmatches { R R vflg { (### ENDmatches) msg } if } D
/beginextra { vflg { (### BEGINextra) msg } if S AXL QYB neg T S 0 0 AX QY 4 copy bbox QBc cmyk fill Blw obox fg cmyk Blw slw K R /Yscale QY D /YH QY D /YY QO D S KQXb { QY nucltick } if R S QBLb { DRFTlh QBlz add QBLp add neg PY 2 div 90 QBLl (bv) (ch) QBLf fg ttshw } if QBLb { DRFTlh neg PY 2 div 90 QBll (bv) (ch) QBlf fg ttshw } if R S } B
/endextra { R R vflg { (### ENDextra) msg } if } D
end % aplotdict
%
%%EndProcSet:   Aplot_dict 1.0 0
%
%%BeginProcSet: Openings 1.0 0
%
/DoInit { tflg { (%%% START NEW PLOT %%%) msg /vflg _t D } if aplotdict begin Dpage 0 0 bbox S bg cmyk fill R clip newpath FXO FYO T tflg { S 0 0 FL neg bbox fg cmyk K R } if flgcrd { shwcrd } if } B
/DoEnd { end tflg { (%%% PLOT DONE %%%) msg } if } B
%
%%EndProcSet:   Openings 1.0 0
%
%
%%EndProlog
%
%%BeginSetup
%
% initgraphics % must not be used in EPS documents
% _t setpacking
_t setstrokeadjust
0.125 setlinewidth
0 setlinejoin
0 setlinecap
%
%%EndSetup
%
+++MAINProcs+++
    &report('PSCODE');
}    # ps_main

sub ps_open_page() {
    print STDOUT << '+++OPEN+++';
%%Page: 1 1
%%BeginPageSetup
%
% Saving current page settings
/pgsave save D
%
%%EndPageSetup
%
+++OPEN+++
}    # ps_open_page

sub ps_close_page() {
    print STDOUT << '+++CLOSE+++';
%
DoEnd
%
%%%%%%%% PLOT DONE
%
% grestoreall % must not be used in EPS documents
pgsave restore
showpage
%
% Page: 1 1
%%PageTrailer
%
+++CLOSE+++
}    # ps_close_page

sub ps_trailer() {
    my $vr = \%{ $Vars{LAYOUT} };
    print STDOUT << "+++EOF+++";
%%Trailer
%
%%Pages: 1
%%Orientation: Portrait
%%BoundingBox: 0 0 $vr->{_page_width} $vr->{_page_height}
%%EOF
+++EOF+++
}    # ps_trailer
