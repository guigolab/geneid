#!/usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
#line 1960 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
# $Id: getfastamasked.pl,v 1.1 2001-10-01 14:09:27 jabril Exp $
#line 1786 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
#
use strict;
#line 1446 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
#
# getfastamasked.pl < fastafile > GFFfile
#
#     Retrieving masked regions coords from fasta files
#
use Bio::Seq;
use Bio::SeqIO;
#line 1807 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
use Benchmark;
  
#line 1815 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
my @Timer = (new Benchmark);
#line 1454 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
my $PROGRAM = 'getfastamasked.pl';
my ($T,$F) = (1,0);
my $DATE = localtime;
my $USER = defined($ENV{USER}) ? $ENV{USER} : '??????';
my $host = `hostname`;
chomp($host);
my $line = ('#' x 80)."\n";
my $s = '### ';
#
my ($id,$seq) = ('','');
my ($total_time);

print STDERR << "+++EOR+++";
$line$s\n$s Running $PROGRAM\n$s
$s HOST: $host
$s USER: $USER
$s DATE: $DATE\n$s\n$line$s
+++EOR+++

&main();

$total_time = &timing($T);
print STDERR << "+++EOR+++";
$s\n$line$s\n$s $PROGRAM FINISHED\n$s
$s TOTAL TIME: $total_time\n$line
+++EOR+++

exit(0);

sub main() {
    my $seqin = Bio::SeqIO->new(-format => 'FASTA', -fh => \*STDIN);
    while (my $sequence = $seqin->next_seq()) {
        my ($sid,$len,$seq,@nuc,@coords,$masked_flg,$match,$msk_num);
        @coords = ();
        
#line 1497 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
print STDERR "### READING FASTA............\n";
$sid  = $sequence->display_id();
$len  = $sequence->length();
$seq  = $sequence->seq();
#line 1489 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
        
#line 1504 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
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
#line 1490 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
        
#line 1527 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
$msk_num = scalar(@coords) / 2;
print STDERR "###         WRITING GFF COORDS: $msk_num masked regions found.\n";
for (my $n = 0; $n <= $#coords; $n+=2) {
    my $GFFstring = ("\%s\t" x 5).(".\t" x 3).".\n";
    printf STDOUT $GFFstring, $sid, "masked", "masked", @coords[$n..($n + 1)];
}; # for coords in @coords
#line 1491 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
    }; # while 
} # main
#line 1819 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
sub timing() {
    push @Timer, (new Benchmark);
    # partial time 
    $_[0] || 
        (return timestr(timediff($Timer[$#Timer],$Timer[($#Timer - 1)])));
    # total time
    return timestr(timediff($Timer[$#Timer],$Timer[0]));
} # timing
