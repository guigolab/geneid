#!/usr/local/bin/perl

use warnings;
use Tk;
use Tk::Tree;
use Data::Dumper;

my ($mv, $tree, $dirfile, @dirs, $ary, $last_lvl, @stack);

@dirs = ();
$dirfile = shift @ARGV || undef;
exit(1) unless defined($dirfile);

$last_lvl = -1;
open(DIRS,"< $dirfile");
while (<DIRS>) {
    my ($pos,$lvl,$data,$full,@F);
    next if /^\s*$/o;
    next if /files\s*$/o;
    chomp;
    $pos = index($_, "-- ") + 3;
    $pos = $pos > 2 ? $pos : 0;
    $lvl = $pos / 4;
    $data = substr($_, $pos);
    @F = split /\s+/o, $data;
    if ($F[0] eq '.') { $full = "SEQ";
    } else {
        $full = $F[0];
    };
    $full =~ s/\./_/og;
    $lvl > $last_lvl && do { push @stack, $full; };
    $lvl == $last_lvl && do { $stack[$lvl] = $full; };
    $lvl < $last_lvl && do { @stack = (@stack[0..($lvl - 1)], $full); };
    # $data =~ s/^([^\s]+\s)/$1.(" " x $lvl)/oe;
    push @dirs, [ join(".", @stack), $data ]; # @F ];
    $last_lvl = $lvl;
    print "$lvl ($last_lvl) -> @F : ".join(".", @stack)."\n";
};
close(DIRS);
# print STDERR Data::Dumper->Dump([ \@dirs ],[ qw/ *dirs / ]);

$mw = MainWindow->new(-title => 'SEQ USAGE');

$tree = $mw->Tree(-font => [ -family => 'courier', -size => 8 ]
                  )->pack(-fill => 'both', -expand => 1);

foreach $ary (@dirs) {
    # my ($branch, $name, $size, $ndirs, $nfiles, $txt);
    # ($branch, $name, $size, $ndirs, $nfiles) = @$ary;
    # $txt = " $name\t$size\t$ndirs\t$nfiles";
    my ($branch, $txt);

    ($branch, $txt) = @$ary;

    $tree->add($branch, -text => "$txt");
};

$tree->autosetmode();
MainLoop;

# exit(0);
