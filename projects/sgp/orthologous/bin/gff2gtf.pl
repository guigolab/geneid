#!/usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
#line 1960 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
# $Id: gff2gtf.pl,v 1.1 2001-10-01 14:09:27 jabril Exp $
#line 1786 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
#
use strict;
#line 1632 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
#
# gff2gtf.pl < infile > outfile
#
#     Converting GFF formated records into GTF
#
# use lib qw( /usr/lib/perl5/site_perl/5.005/ /usr/lib/perl5/5.00503/ ) ;
use Benchmark;
# use Data::Dumper;
# local $Data::Dumper::Purity = 0;
# local $Data::Dumper::Deepcopy = 1;
my ($T,$F) = (1,0);
my @Timer = (new Benchmark);
my $PROGRAM = 'gff2gtf.pl';
my $DATE = localtime;
my $USER = defined($ENV{USER}) ? $ENV{USER} : '??????';
my $host = `hostname`;
chomp($host);
my $line = ('#' x 80)."\n";
my $s = '### ';
#
print STDERR << "+++EOR+++";
$line$s\n$s Running $PROGRAM\n$s
$s HOST: $host
$s USER: $USER
$s DATE: $DATE\n$s\n$line$s
+++EOR+++

my %genes;
my @f = ();
my $gcnt = 0;
while (<>) {
    /^\# Gene/o && do {
        $gcnt++;
        $genes{$gcnt}{COMMENT} = $_;
        next;
    };
    next if /^\#/o;
    next if /^\s*$/o;
    chomp;
    @f = split /\s+/og, $_;
    push @{ $genes{$gcnt}{FEATURES} }, [ @f[0..8] ];
}; # while

# print STDERR Dumper(\%genes);

my $GTFrec = ("\%s\t" x 8)."\%s\n";
foreach my $i (1..$gcnt) {
    print STDOUT $genes{$i}{COMMENT};
    my ($main,$tail_O,$tail_E,$groupflg,$gtfgroup) = ('','','',$T,'');
    foreach my $j (0..$#{ $genes{$i}{FEATURES} }) {
        my ($seq,$source,$feat,$start,$end,
            $score,$strand,$frame,$group,
            $cdsfeat,$cdsori,$cdsend,$t);
        ($seq,$source,$feat,$start,$end,
            $score,$strand,$frame,$group) =
                @{ $genes{$i}{FEATURES}[$j] };
        $score = '.';
        $groupflg && do {
            $group =~ s/gene_//o;
            $t = &fill_left($group,3,"0");
            $gtfgroup = "gene_id $t; transcript_id $t.1";
            $groupflg = $F;
        }; # $groupflg
      FEATS: {
          $feat eq 'Single' && do {
              ($cdsfeat,$cdsori,$cdsend) = &get_start($start,$end,$strand);
              $tail_O = sprintf($GTFrec,$seq,$source,$cdsfeat,$cdsori,$cdsend,
                                $score,$strand,'0',$gtfgroup);
              ($cdsfeat,$cdsori,$cdsend) = &get_final($start,$end,$strand);
              $tail_E = sprintf($GTFrec,$seq,$source,$cdsfeat,$cdsori,$cdsend,
                                $score,$strand,'0',$gtfgroup);
              last FEATS;
          }; # Single
          $feat eq 'First' && do {
              ($cdsfeat,$cdsori,$cdsend) = &get_start($start,$end,$strand);
              $tail_O = sprintf($GTFrec,$seq,$source,$cdsfeat,$cdsori,$cdsend,
                                $score,$strand,'0',$gtfgroup);              
              last FEATS;
          }; # First
          $feat eq 'Terminal' && do {
              ($cdsfeat,$cdsori,$cdsend) = &get_final($start,$end,$strand);
              $tail_E = sprintf($GTFrec,$seq,$source,$cdsfeat,$cdsori,$cdsend,
                                $score,$strand,'0',$gtfgroup);
                  # only if nucleotides of stop codon not included
                  #     $strand eq '+' && ($end = $end + 3);
                  #     $strand eq '-' && ($start = $start - 3);
              # last FEATS;
          }; # Terminal
        }; # FEATS
        $feat = 'CDS';
        $main .= sprintf($GTFrec,$seq,$source,$feat,$start,$end,
                         $score,$strand,$frame,$gtfgroup); 
    }; # foreach $j
    print STDOUT "$main$tail_O$tail_E";
}; # foreach $i

my $total_time = &timing($T);
print STDERR << "+++EOR+++";
$s\n$line$s\n$s $PROGRAM FINISHED\n$s
$s TOTAL TIME: $total_time\n$line
+++EOR+++

exit(0);

sub get_start() {
    my ($o,$e,$s) = @_;
    my $str = "start_codon";
    $s eq '+' && do {
       return $str, $o, ($o + 2);
    }; # forward
    $s eq '-' && do {
       return $str, ($e - 2), $e;
    }; # reverse
    die("### ERROR: Strand not defined... ($o $e : $s)... $!");
} # 
#
sub get_final() {
    my ($o,$e,$s) = @_;
    my $str = "stop_codon";
    $s eq '+' && do {
        return $str, ($e - 2), $e;
           # only if nucleotides of stop codon not included
           #     return $str, ($e + 1), ($e + 3);
    }; # forward
    $s eq '-' && do {
        return $str, $o, ($o + 2);
           # only if nucleotides of stop codon not included
           #     return $str, ($o - 3), ($o -1);
    }; # reverse
    die("### ERROR: Strand not defined... ($o $e : $s)... $!");
} # 
#
sub fill_left() { ($_[2] x ($_[1] - length($_[0]))).$_[0] } 
#
sub timing() {
    push @Timer, (new Benchmark);
    # partial time
    $_[0] ||
        (return timestr(timediff($Timer[$#Timer],$Timer[($#Timer - 1)])));
    # total time
    return timestr(timediff($Timer[$#Timer],$Timer[0]));
} # timing
