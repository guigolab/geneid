#!/usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
#line 796 "/home/jabril/development/webprecarios/D-recerca_web.nw"
# $Id: HtmlSiteBuilder.pl,v 1.1 2001-07-27 13:44:02 jabril Exp $
#line 636 "/home/jabril/development/webprecarios/D-recerca_web.nw"
#
use strict;
#line 255 "/home/jabril/development/webprecarios/D-recerca_web.nw"
#
# HtmlSiteBuilder.pl
#
#     Building web pages sharing the same page layout,
#     global paths to pages can be defined separately from HTML code. 
#
# Usage:
# 
#     HtmlSiteBuilder.pl <OutputPath> <HierarchyFile> <ContentsFile> \
#                                               [ ... <ContentsFile> ]
#
#<Use Modules - HSB>>
#<Use Modules - Dumper>>
#
#line 277 "/home/jabril/development/webprecarios/D-recerca_web.nw"
my $PROGRAM = 'HtmlSiteBuilder.pl';
my $USAGE = '<OutputPath> <HierarchyFile> <ContentsFile> [ ... <ContentsFile> ]';
#line 293 "/home/jabril/development/webprecarios/D-recerca_web.nw"
my ($output_dir,$hierarchy_file,@contents_file);
#line 323 "/home/jabril/development/webprecarios/D-recerca_web.nw"
my (%HTMLvar,%HTMLtree);
#line 270 "/home/jabril/development/webprecarios/D-recerca_web.nw"
#
#line 282 "/home/jabril/development/webprecarios/D-recerca_web.nw"
print STDERR "###\n### RUNNING $PROGRAM\n###\n".
             "### $ENV{USER} - ".(`date`)."###\n#\n";
&parse_args();
&read_hierarchy($hierarchy_file);
foreach my $file (@contents_file) {
    &read_contents($file);
}; # foreach
exit(0);
#line 272 "/home/jabril/development/webprecarios/D-recerca_web.nw"
#
#line 297 "/home/jabril/development/webprecarios/D-recerca_web.nw"
sub parse_args() {
    scalar(@ARGV) < 3 && do {
        print "\nUSAGE: \n    $PROGRAM $USAGE\n\n";
        die("!!! ERROR - NOT ENOUGH COMMAND-LINE PARAMETERS. $!\n");
    };
    $output_dir = shift @ARGV;
    $output_dir =~ s%/$%%o;
    ( -e $output_dir && -d _ ) ||
        die("!!! $output_dir DOES NOT EXIST...\n");
    $hierarchy_file = shift @ARGV;
    ( -e $hierarchy_file ) ||
        die("!!! $hierarchy_file DOES NOT EXIST...\n");
    foreach my $file (@ARGV) {
        ( -e $file ) ||
            print STDERR "!!! $file DOES NOT EXIST...\n";
        push @contents_file, $file;
    }; # foreach
    scalar(@contents_file) == 0 && do {
        print "\nUSAGE: \n    $PROGRAM $USAGE\n\n";
        die("!!! No content files were given...\n\n");
    };
    @ARGV = ();
} # parse_args
#line 327 "/home/jabril/development/webprecarios/D-recerca_web.nw"
sub read_hierarchy() {
    my $infile = $_[0];
    print STDERR "### READING WEB HIERARCHY FROM: $infile\n";
    open(HFILE,"< $infile") ||
        die("!!! CANNOT OPEN FILE: $infile $!");
    while (<HFILE>) {
        my @f;
        
#line 695 "/home/jabril/development/webprecarios/D-recerca_web.nw"
next if /^\#/o;
next if /^\s*$/o;
chomp;
#line 335 "/home/jabril/development/webprecarios/D-recerca_web.nw"
        @f = split /\s+/og;
        $f[0] =~ m/^>#</ && do {  # HTML main variables definition
            scalar(@f) < 3 && (next);
            $HTMLvar{$f[1]} = $f[2];
            next;
        }; #             # HTML main variables definition
        scalar(@f) < 3 && do {
            
        };
        $HTMLtree{$f[0]}{PARENT} = $f[1];
        $HTMLtree{$f[0]}{BASEHREF} = $f[2];
        push @{ $HTMLtree{$f[1]}{CHILD} }, $f[0];
    }; # while
    close(HFILE);
    print STDERR "### READING WEB HIERARCHY... DONE!!!\n";
    &show_hierarchy('ROOT',1);
} # read_hierarchy
#line 355 "/home/jabril/development/webprecarios/D-recerca_web.nw"
sub show_hierarchy() {
    my ($node,$tab) = @_;
    print STDERR ('   |' x $tab)."-- $node\n";
    return unless defined($HTMLtree{$node});
    foreach my $child (@{ $HTMLtree{$node}{CHILD} }) {
        &show_hierarchy($child,($tab + 1));
    }; # foreach
} # show_hierarchy
#line 366 "/home/jabril/development/webprecarios/D-recerca_web.nw"
sub read_contents() {
    my $infile = $_[0];
    print STDERR "### PROCESSING WEB CONTENTS FROM: $infile\n";
    open(CFILE,"< $infile") ||
        die("!!! CANNOT OPEN FILE: $infile $!");
    while (<CFILE>) {
        my @f;
        
#line 695 "/home/jabril/development/webprecarios/D-recerca_web.nw"
next if /^\#/o;
next if /^\s*$/o;
chomp;
#line 374 "/home/jabril/development/webprecarios/D-recerca_web.nw"
    }; # while
    close(CFILE);
    print STDERR "### PROCESSING WEB CONTENTS... DONE!!!\n";
} # read_contents
