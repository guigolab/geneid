#!/usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
#line 1960 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
# $Id: coords2gff.pl,v 1.1 2001-10-01 14:09:27 jabril Exp $
#line 1786 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
#
use strict;
#line 1363 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
#
# coords2gff.pl seqname < coords_file > GFF_file
#
#     Converting raw coords files to GFF
#
# use Data::Dumper;
# local $Data::Dumper::Purity = 0;
# local $Data::Dumper::Deepcopy = 1;

my ($seqname,%gff,$string);
$seqname = shift @ARGV;
%gff = ();
$string = ("\%s\t" x 8)."\%s\n";

my $c = 0;
while (<STDIN>) {
    next if /^\s*$/o;
    my ($l,$a,$b,$s,$r);
    chomp;
    $l = $_;
	# print STDERR "$. (X): $l\n";
    $l =~ /^[><]/o && do {
	    # print STDERR "$. (A): $l\n";
        $l =~ s/^([><])\s*//o && ($s = $1);
        $r = \%{ $gff{++$c} };
        ($r->{gn_start},$r->{gn_end},$r->{gn_name}) = split /\s+/og, $l, 3;
        $r->{gn_name} =~ s/\s+$//og;
        $r->{gn_name} =~ s/\s+/_/og;
        $r->{strand} = ($s eq '>') ? '+' : (($s eq '<') ? '-' : '.');
        # $r->{exons} = ();
        next;
	}; # $newgene ?
	# print STDERR "$. (B): $l\n";
    ($a,$b,undef) = split /\s+/og, $l, 3;
    push @{ $gff{$c}{exons} }, [ $a, $b ];
}; # while
# print STDERR Dumper(\%gff);

foreach my $i (1..$c) {
    my $r = \%{ $gff{$i} };
    my $exnum = scalar(@{ $r->{exons} });
    print STDOUT "\# $seqname - $r->{gn_name}: ".
        ($exnum)." exons\n";
    printf STDOUT $string, 
        $seqname,'annotation','Gene',$r->{gn_start},$r->{gn_end},
        '.',$r->{strand},'.',$r->{gn_name};
    $r->{strand} eq '-' && (@{ $r->{exons} } = reverse @{ $r->{exons} });
    my ($lastori,$lastend,$lastfrm) = (0,0,0);
    foreach my $j (0..$#{ $r->{exons} }) {
        my ($feat,$ori,$end,$frm);
        ($ori,$end) = @{ ${ $r->{exons} }[$j] };
	  THIS: {
        ($j == 0) && do {
            $feat = 'First';
            $exnum == 1 && ($feat = 'Single');
            $frm = 0;
            last THIS;
        }; # $d == 0
        $frm = (($lastend - $lastori + 1) + $lastfrm) % 3;
        ($j == $#{ $r->{exons} }) && do {
            $feat = 'Terminal';
            last THIS;
        }; # $d == $#{ $r->{exons} }
        $feat = 'Internal';
      }; # THIS
        # print STDERR "$r->{gn_name} ($j): $ori $end $frm : ".
        #              "$lastori $lastend $lastfrm\n";
        ($lastori,$lastend,$lastfrm) = ($ori,$end,$frm);
        printf STDOUT $string, 
            $seqname,'annotation',$feat,$ori,$end,
            '.',$r->{strand},$frm,$r->{gn_name};
    }; # foreach @j
}; # foreach $i

exit(0);
