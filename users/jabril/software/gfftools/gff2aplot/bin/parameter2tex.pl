#!/usr/local/bin/perl -w
# This is perl, v5.6.1 built for i686-linux
# /usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
#line 14127 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
# 
# parameter2tex.pl
#
#line 15952 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
# $Id: parameter2tex.pl,v 1.2 2003-03-04 18:05:40 jabril Exp $
#line 14131 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
#
use strict;
use Getopt::Std;
use File::Basename;
#
#line 14151 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
my $USAGE = << "+++EOU+++";
################################################################################
###
#### $0
###
####\t\t@{[ (localtime)." : ".(defined($ENV{USER})?$ENV{USER}:"nouser") ]} 
###
#### USAGE:
###
###    parameter2tex [options] -o <output_file> -- <input_files>
###
###    -h              prints this help.
###    -v              process custom-file vars input.
###    -c              process command-line options input.
###    -o <path+file>  output base name for the whole path/file (required).
###
################################################################################
+++EOU+++
#
my %VARS = ();
my ($base_dir,$base_file,$base_ext) = ('./','out','.tex');
my $parseinput = \&parseVars;
#line 14214 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
my @files;
#line 14443 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
local (*LAYa,*LAYb,*LAYc,*LAYd,*LAYe,*LAYf,*LAYg,*LAYh,
       *SEQa,*SEQb,*SRCa,*SRCb,*STRa,*STRb,*GRPa,*GRPb,*FEAa,*FEAb,
       *XTRa,*XTRb);
my %varO = (
    #                    [ HANDLE, EXTENSION, LABEL                     ]
    LAYOUT            => [ \*LAYa, 'Layout'                             ],
    LAYOUTPageLayout  => [ \*LAYb, 'playout', 'Page Layout'             ],
    LAYOUTZoom        => [ \*LAYc, 'zoom',    'Zoom Options'            ],
    LAYOUTLabelsA     => [ \*LAYd, 'labela',  'Labels'                  ],
    LAYOUTLabelsB     => [ \*LAYe, 'labelb',  'Labels (Continued)'      ],
    LAYOUTTickmarks   => [ \*LAYf, 'ticks',   'Tick Marks'              ],
    LAYOUTAplot       => [ \*LAYg, 'aplot',   'Aplot Layout'            ],
    LAYOUTGeneral     => [ \*LAYh, 'general', 'General Definitions'     ],
    #
    SEQUENCE          => [ \*SEQa, 'Sequence'                           ],
    SEQUENCESequences => [ \*SEQb, 'seq',     'GFF-Sequence Attributes' ],
    #
    SOURCE            => [ \*SRCa, 'Source'                             ],
    SOURCESources     => [ \*SRCb, 'src',     'GFF-Source Attributes'   ],
    #
    STRAND            => [ \*STRa, 'Strand'                             ],
    STRANDStrands     => [ \*STRb, 'str',     'GFF-Strand Attributes'   ],
    #
    GROUP             => [ \*GRPa, 'Group'                              ],
    GROUPGroups       => [ \*GRPb, 'grp',     'GFF-Group Attributes'    ],
    #
    FEATURE           => [ \*FEAa, 'Feature'                            ],
    FEATUREFeatures   => [ \*FEAb, 'feat',    'GFF-Feature Attributes'  ],
    #
    EXTRA             => [ \*XTRa, 'Extra'                              ],
    EXTRAExtra        => [ \*XTRb, 'extra',   'Special Customization Features'  ],
    );
my @varMain = qw / LAYOUT SEQUENCE SOURCE STRAND GROUP FEATURE EXTRA /;
my @varSub  = qw /
                   LAYOUTPageLayout  LAYOUTZoom
                   LAYOUTLabelsA     LAYOUTLabelsB
                   LAYOUTTickmarks   LAYOUTAplot     LAYOUTGeneral
                   SEQUENCESequences SOURCESources   STRANDStrands
                   GROUPGroups       FEATUREFeatures EXTRAExtra
                 /;
#line 14137 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
#
# MAIN
&getcmdlineopts;
&$parseinput;
exit(0);
#
# SUBS
#line 14178 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
sub getcmdlineopts() {
    our($opt_v,$opt_c,$opt_o,$opt_h);
    getopts('vco:h');
    $opt_h && do {
        print STDERR $USAGE;
        exit(1);
    }; # $opt_h
    $opt_v && ($parseinput = \&parseVars);
    $opt_c && ($parseinput = \&parseOpts);
    ($base_file,$base_dir,$base_ext) = fileparse($opt_o, '\..*');
    # print STDERR "::$base_dir\::$base_file\::$base_ext\::\n";
    (-d $base_dir) ||
        die("##\n## ERROR ## OUTPUT DIRECTORY DOES NOT EXIST... $!\n##\n"); 
    $base_file eq '' && ($base_file = 'out');
	$base_ext eq '' && ($base_ext = '.tex');
    #
    @files = ();
    scalar(@ARGV) > 0 && do {
        foreach my $fl (@ARGV) {
            $fl eq '-' && do {
                push @files, "-";
                next;
            };
            -e $fl || do {
                print STDERR "### FILE $fl DOES NOT EXIST, SKIPPING... $!\n";
                next;
            };
            push @files, $fl;
		}; # foreach $fl
        return;
    }; # scalar(@ARGV) > 0
    push @files, "-";
} # getcmdlineopts
#line 14222 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
sub parseinput() {
    print STDERR "### Parsing input records...\n";
    my $char;
    my $c = 0;
    my $extra = 1000;
    foreach my $fl (@files) {
        my ($a,$b,$lstelm,$srtelm,$order);
        open(IFILE,"< $fl");
        $lstelm = $srtelm = 0;
        $order = undef;
        while (<IFILE>) {
            next if /^\s*$/o;
            chomp;
            $_ =~ s/\s*$//o;
            $char = '.';
            
#line 14249 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
$_ =~ /^\#\#\#EOR\#\#\#/o && do {
    $lstelm = $srtelm = 0;
    $order = undef;
    print STDERR "\#\n";
    $char = ''; # $char = '#';
    next;
}; # ###EOR###
$_ =~ /^LDE/o && do {
    $_ =~ /JOIN-NEXT\s*$/o && ($VARS{$order}{JOIN} = 1);
    $_ =~ s/^LDE:(\s+JOIN-NEXT\s+)?//o;
    $lstelm = 1;
    $srtelm = 0;
}; # LDE
$lstelm && do {
    $_ =~ s/^\s*//o;
    $char = ':'; 
    $VARS{$order}{LDE} .= " $_" unless $_ eq "";
    next;
};
$_ =~ /^SDE/o && do {
    $_ =~ /JOIN-NEXT\s*$/o && ($VARS{$order}{JOIN_SHORT} = 1);
    $_ =~ s/^SDE:(\s+JOIN-NEXT\s+)?//o;
    $srtelm = 1;
}; # SDE
$srtelm && do {
    $_ =~ s/^\s*//o;
    $char = '+'; 
    $VARS{$order}{SDE} .= " $_" unless $_ eq "";
    next;
};
#line 14238 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
            
#line 14284 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
($a,$b) = split /:\s*/o, $_, 2;
defined($order) || do {
    CHECK: {
      $a =~ /^ORD/o && do {
          $order = (defined($b) && $b ne '') ? $b : ++$extra;
          $char = 'O';
          last CHECK;
      }; # ORD
      $order = ++$extra;
      $char = 'X';
    }; # CHECK
    $order = &fill_left($order,4,'0');
    %{ $VARS{$order} } = ();        
    $VARS{$order}{JOIN} = 0;
    $VARS{$order}{JOIN_SHORT} = 0;
    print STDERR &fill_left(++$c,3,'0')." ";
    next; 
}; # $order
(defined($b) && $b ne '') && do {
    $VARS{$order}{$a} = $b;
}; # $b
#line 14239 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
		} continue {
            print STDERR $char; # &counter(++$c,$char);
        }; # while read $fl
        close(IFILE);
    }; # foreach $fl
    print STDERR "\#\n";        # &counter_end($c,$char);
} # parseinput
#line 14310 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
sub rep_varop()  { $_[0] =~ s/\[(.+?)\]/\\op\{$1\}/og;     return $_[0]; } 
# sub repq_varop() { $_[0] =~ s/\[(.+?)\]/\'\\op\{$1\}\'/og; return $_[0]; } 
sub repq_varop() { $_[0] =~ s/\[(.+?)\]/\\op\{$1\}/og;     return $_[0]; } 
#line 14316 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
sub rep_param()  { $_[0] =~ s/<(.+?)>/\\pp\{$1\}/og;     return $_[0]; }
# sub repq_param() { $_[0] =~ s/<(.+?)>/\'\\pp\{$1\}\'/og; return $_[0]; }
sub repq_param() { $_[0] =~ s/<(.+?)>/\\pp\{$1\}/og;     return $_[0]; }
#line 14322 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
sub rep_value()  { $_[0] =~ s/\|(.+?)\|/\\vp\{$1\}/og;     return $_[0]; }
# sub repq_value() { $_[0] =~ s/\|(.+?)\|/\'\\vp\{$1\}\'/og; return $_[0]; }
sub repq_value() { $_[0] =~ s/\|(.+?)\|/\\vp\{$1\}/og;     return $_[0]; }
#line 14328 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
sub rep_chars()  { 
    $_[0] =~ s/\%\{/\[/og;  # to recover LaTeX []
    $_[0] =~ s/\%\}/\]/og;  # 
    $_[0] =~ s/\\n\\/\n/og; # to recover LaTeX newline
    return $_[0];
} # rep_chars
#line 14337 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
sub tex_header() {
    my ($a,$b) = @_;
    $a =~ m{/([^/]*)$}o && ($a = $1);
    return "\%\n\% $a\n\%\n\% $b\n\%\n\%".
           ' $Id: parameter2tex.pl,v 1.2 2003-03-04 18:05:40 jabril Exp $ '.
           "\n\%\n"; 
} # tex_header
#line 14365 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
sub parseVars() {
    my ($k,$c);
    print STDERR "###\n### PARSING CUSTOM-FILE VARS DEFINITIONS\n###\n";
    &parseinput;
    print STDERR "### Opening Files ${base_dir}${base_file}\_\*${base_ext}\n";
    &open_varfiles;
    $c = 0;
    foreach $k (sort keys %VARS) {
        print STDERR "$k".((++$c % 10) ? ".." : "\n");
        ((defined($VARS{$k}{SEC}) && $VARS{$k}{SEC} ne '')
                                  &&
         (defined($VARS{$k}{SUB}) && $VARS{$k}{SUB} ne '')) || do {
             print STDERR "\n##\n## ERROR ## SECTION/SUBSECTION ".
                          "NOT DEFINED ON INPUT for $k ... $!\n##\n";
             next;
         }; # SEC/SUB NOT DEFINED ON INPUT
        my ($sec,$subsec, 
            $tmpdsc,$strdsc,$tmptbl,$strtbl,$par,$vareq,$vareqj);
        ($sec,$subsec) = ($VARS{$k}{SEC}, "$VARS{$k}{SEC}$VARS{$k}{SUB}");
        (defined($varO{$sec}) && defined($varO{$subsec})) || do {
             print STDERR "\n##\n## ERROR ## SECTION/SUBSECTION ".
                          "NOT DEFINED ON PROGRAM for $k ... $!\n##\n";
             next;
         }; # SEC/SUB NOT DEFINED ON PROGRAM
        if (defined($VARS{$k}{PAR}) && $VARS{$k}{PAR} ne '') {
            $par = (&rep_chars(&rep_varop(&rep_param($VARS{$k}{PAR}))));
        } else {
            $par = '';
        }; # PAR
        (defined($VARS{$k}{OPT}) && $VARS{$k}{OPT} ne '') || do {
             print STDERR "\n##\n## ERROR ## VARIABLE NAME ".
                          "NOT DEFINED for $k ... $!\n##\n";
             next;         
        }; # DEFINED OPT
        #
        $tmpdsc = '\op{'.$VARS{$k}{OPT}.'} $\,=\,$ '.$par;
        $tmptbl = '\op{'.$VARS{$k}{OPT}.'}';
        #
        $vareqj = ' \hfill ';
        if (defined($VARS{$k}{DEF})) {
            $vareq  = (&rep_chars(&rep_varop(
                           &rep_param(&rep_value($VARS{$k}{DEF})))));
        } else {
            $vareq  = '\bydef';
        }; # DEF
        defined($VARS{$k}{SDE}) || ($VARS{$k}{SDE} = '{\tbdef}');
        defined($VARS{$k}{LDE}) || ($VARS{$k}{LDE} = ' {\tbdef}');
        if ($VARS{$k}{JOIN_SHORT}) {
            $strtbl = '\rvjoin{'.$tmptbl.'}{'.$vareq.'}'."\n\%\n";
        } else {
            $strtbl = '\rvdesc{'.$tmptbl.'}{'.$vareq.'}'."\n".'   { '.
                (&rep_chars(&repq_varop(
                     &repq_param(&repq_value($VARS{$k}{SDE}))))).' }'.
                "\n\%\n";
        }; # JOIN_SHORT
        if ($VARS{$k}{JOIN}) {
            $strdsc = '\ijoin{'.$tmpdsc.$vareqj.'[ '.$vareq.' ]}'."\n\%\n";
        } else {
            $strdsc = '\idesc{'.$tmpdsc.$vareqj.'[ '.$vareq.' ]}'."\n".'   {'.
                (&rep_chars(&repq_varop(
                     &repq_param(&repq_value($VARS{$k}{LDE}))))).' }'.
                "\n\%\n";
        }; # JOIN
        #
        $strdsc =~ s/(_|\%\-)/\\_/og;
        $strtbl =~ s/(_|\%\-)/\\_/og;
        print { $varO{$sec}[0] }    $strdsc;
        print { $varO{$subsec}[0] } $strtbl;
    }; # foreach keys %VARS
    ($c % 10) && print STDERR "\#\#\#\#\n";
    print STDERR "\#\#\#\n";
    &close_varfiles;
} # parseVars
#line 14495 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
sub open_varfiles() {
    my ($flh,$vlh,$mf);
    foreach $flh (@varMain) {
        $vlh = "_dsc".$varO{$flh}[1];
        $mf = "${base_dir}${base_file}${vlh}${base_ext}";
        open($varO{$flh}[0],"> $mf");
        print { $varO{$flh}[0] } 
              &tex_header($mf,"Description items for $flh");
	}; # foreach $flh
    foreach $flh (@varSub) {
        $vlh = "_tbl".$varO{$flh}[1];
        $mf = "${base_dir}${base_file}${vlh}${base_ext}";
        open($varO{$flh}[0],"> $mf");
        print { $varO{$flh}[0] } 
              &tex_header($mf,"Summary table for $flh: $varO{$flh}[2]");
        print { $varO{$flh}[0] } "\%\n".'\begin{tabular}{p{5cm}p{3cm}p{15cm}}';
        print { $varO{$flh}[0] } "\n\%\n".'\rvdef{'.$varO{$flh}[2].'}'."\n\%\n";
	}; # foreach $flh
} # open_varfiles
sub close_varfiles() {
    my $flh;
    foreach $flh (@varMain) {
        close($varO{$flh}[0]);
	}; # foreach $flh
    foreach $flh (@varSub) {
        print { $varO{$flh}[0] } '\end{tabular}'."\n";
        close($varO{$flh}[0]);
	}; # foreach $flh
} # close_varfiles
#line 14544 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
sub parseOpts() {
    my ($k,$c);
    print STDERR "###\n### PARSING CMD-LINE OPTIONS DEFINITIONS\n###\n";
    &parseinput;
    print STDERR "### Writing to ${base_dir}${base_file}\_\*${base_ext}\n";
    open(FDSC,"> ${base_dir}${base_file}\_dsc${base_ext}");
    open(FTBL,"> ${base_dir}${base_file}\_tbl${base_ext}");
    $c = 0;
    print FTBL '\begin{tabular}{rl}'."\n\%\n";
    foreach $k (sort keys %VARS) {
        print STDERR "$k".((++$c % 10) ? ".." : "\n");
        my ($tmpdsc,$strdsc,$tmptbl,$strtbl,$par,$varlng,$varsht,$vareq);
        if (defined($VARS{$k}{PAR}) && $VARS{$k}{PAR} ne '') {
            $par = &rep_param($VARS{$k}{PAR});
        } else {
            $par = '';
        }; # PAR
        (defined($VARS{$k}{LNG}) && $VARS{$k}{LNG} ne '') || do { 
             print STDERR "\n##\n## ERROR ## LONG OPTION NAME ".
                          "NOT DEFINED for $k ... $!\n##\n";
             next;
        }; # LONG OPTION NOT DEFINED
        $varlng = '\op{-\/-'.$VARS{$k}{LNG}.'} '.$par;
        if (defined($VARS{$k}{OPT}) && $VARS{$k}{OPT} ne '') { 
            $varsht  = '\op{-'.$VARS{$k}{OPT}.'}';
            $tmptbl  = $varsht.'{\x}'.$varlng;
            $varsht .= " $par";
            $tmpdsc  = '\shortstack[l]{\ '.$varsht.' \\\\ '.$varlng.'}';
        } else {
            $tmpdsc = $varlng;
            $tmptbl = $varlng;
        }; # OPT
        $vareq = ' \hfill ';
        (defined($VARS{$k}{EQV})) && 
            ( $vareq .= (&rep_varop(&rep_param(&rep_value($VARS{$k}{EQV})))) );
        defined($VARS{$k}{SDE}) ||
            ($VARS{$k}{SDE} = '{\tbdef}');
        defined($VARS{$k}{LDE}) ||
            ($VARS{$k}{LDE} = ' {\tbdef}');
        if ($VARS{$k}{JOIN_SHORT}) {
            $strtbl = '\rjoin{'.$tmptbl.'}'."\n\%\n";
        } else {
            $strtbl = '\rdesc{'.$tmptbl.'}'."\n".'   { '.
                (&rep_chars(&repq_varop(
                     &repq_param(&repq_value($VARS{$k}{SDE}))))).' }'.
                "\n\%\n";
        }; # JOIN_SHORT
        if ($VARS{$k}{JOIN}) {
            $strdsc = '\ijoin{'.$tmpdsc.$vareq.'}'."\n\%\n";
        } else {
            $strdsc = '\idesc{'.$tmpdsc.$vareq.'}'."\n".'   {'.
                (&rep_chars(&repq_varop(
                     &repq_param(&repq_value($VARS{$k}{LDE}))))).' }'.
                "\n\%\n";
        }; # JOIN
        #
        $strdsc =~ s/(_|\%\-)/\\_/og;
        $strtbl =~ s/(_|\%\-)/\\_/og;
        print FDSC $strdsc;
        print FTBL $strtbl;
    }; # foreach keys %VARS
    ($c % 10) && print STDERR "\#\#\#\#\n";
    print STDERR "\#\#\#\n";
    print FTBL '\end{tabular}'."\n";
    close(FDSC);
    close(FTBL);
} # parseOpts
#line 15912 "/home/ug/jabril/development/softjabril/gfftools/gff2aplot/gff2aplot.nw"
#
sub fill_right() { $_[0].($_[2] x ($_[1] - length($_[0]))) }
sub fill_left()  { ($_[2] x ($_[1] - length($_[0]))).$_[0] }
sub fill_mid()   { 
    my $l = length($_[0]);
    my $k = int(($_[1] - $l)/2);
    ($_[2] x $k).$_[0].($_[2] x ($_[1] - ($l+$k)));
} # fill_mid
