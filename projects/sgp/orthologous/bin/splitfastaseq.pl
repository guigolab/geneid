#!/usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
#line 1960 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
# $Id: splitfastaseq.pl,v 1.1 2001-10-01 14:09:27 jabril Exp $
#line 1786 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
#
use strict;
#line 1541 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
#
# splitfastaseq.pl \
#     seqlength overlap < fastafile > output
#
#     Breaking large fasta sequences to build 
#     databases for running tblastx faster
# 
use lib qw( /usr/lib/perl5/site_perl/5.005/ );
use Bio::Seq;
use Bio::SeqIO;
use Benchmark;
my ($T,$F) = (1,0);
my @Timer = (new Benchmark);
my $PROGRAM = 'splitfastaseq.pl';
my $DATE = localtime;
my $USER = defined($ENV{USER}) ? $ENV{USER} : '??????';
my $host = `hostname`;
chomp($host);
my $line = ('#' x 80)."\n";
my $s = '### ';
#
my ($id,$ln,$sq) = ('',0,'');
my ($total_time,$seq);
my ($maxlen,$overlap) = @ARGV;

print STDERR << "+++EOR+++";
$line$s\n$s Running $PROGRAM\n$s
$s HOST: $host
$s USER: $USER
$s DATE: $DATE\n$s\n$line$s
+++EOR+++

&getseq();
&splitseq();

$total_time = &timing($T);
print STDERR << "+++EOR+++";
$s\n$line$s\n$s $PROGRAM FINISHED\n$s
$s TOTAL TIME: $total_time\n$line
+++EOR+++

exit(0);

sub getseq() { # assuming here single sequence input fasta files
    print STDERR "$s Processing fasta file.\n";
    my $seqin = Bio::SeqIO->new(-format => 'FASTA', -fh => \*STDIN);
    while (my $iseq = $seqin->next_seq()) {
        $id = $iseq->display_id();
        $ln = $iseq->length();
        $sq = $iseq->seq();
        last; # to make sure that we only catch a single fasta sequence
    }; # while next_seq
    $seq = Bio::Seq->new( -seq => $sq , -id => $id );
    print STDERR "$s Processing DONE: ".(&timing($F))."\n$s\n";
} # getseq
#
sub splitseq() {
    my ($e,$sid,$ssq,$nseq,$wseq);
    my ($t,$sqlen) = (1,($maxlen + $overlap - 1));
    print STDERR "$s Creating splitted-sequence fasta file ($ln nt).\n";
    my $seqout = Bio::SeqIO->new(-format => 'FASTA', -fh => \*STDOUT);
    while ($t < $ln) {
         $e = $t + $sqlen;
         ($e > $ln) && ($e = $ln);
         $sid = "$id\_$t\_$e";
         print STDERR "$s --> $id : from $t to $e (".($e - $t + 1)." nt)\n";
         $ssq = $seq->subseq($t,$e);
         $t += $maxlen;
         #
         $wseq = Bio::Seq->new( -seq => $ssq , -id => $sid );
         $seqout->write_seq($wseq);
    }; # while  
    print STDERR "$s Splitting DONE: ".(&timing($F))."\n$s\n";
} # splitseq
#
sub timing() {
    push @Timer, (new Benchmark);
    # partial time
    $_[0] ||
        (return timestr(timediff($Timer[$#Timer],$Timer[($#Timer - 1)])));
    # total time
    return timestr(timediff($Timer[$#Timer],$Timer[0]));
} # timing
