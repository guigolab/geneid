#!/usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
#
#line 1475 ".//deploy.nw"
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %                             DEPLOY                               %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
#    .
# 
#     Copyright (C) 2001 - Josep Francesc ABRIL FERRANDO  
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
#line 1307 ".//deploy.nw"
#
#line 1469 ".//deploy.nw"
# $Id: deploy.pl,v 1.1 2001-08-21 13:15:04 jabril Exp $
#line 1309 ".//deploy.nw"
#
use strict;
#line 196 ".//deploy.nw"
#
# MODULES
#

#
# VARIABLES
#
#line 218 ".//deploy.nw"
my $USAGE = "\nUSAGE:\n\tdeploy.pl <projectname>\n".
            "(It asumes that you are in the right directory)\n\n";
my @working_dirs = qw(
                       RCS
                       bin  bin/param
                       data
                       docs docs/psfigures docs/tables docs/html
                       tests
                       );
my $PROJECT;
# my $HOME = $ENV{HOME};
my $CWD  = `pwd`;
chomp($CWD);
my $PATH = $CWD;
# $PATH =~ s%^$HOME/%%o;
#line 204 ".//deploy.nw"
#
# MAIN LOOP
#
#line 236 ".//deploy.nw"
&parse_argvs();

print STDERR "###\n### RUNNING $PROGRAM..........\n###\n".
             "### User: $USER\n".
             "### Date: $DATE\n###\n".
             "### Current Working Directory: $CWD\n".
             "### Setting PATH to: $PATH\n".
             "### Project NAME: $PROJECT\n###\n";

&make_dirs();
&new_noweb_doc();
&extract_files();

print STDERR "###\n### RUNNING deploy.pl............ DONE\n###\n";

exit(0);
#line 208 ".//deploy.nw"
#
# FUNCTIONS
#
#line 255 ".//deploy.nw"
sub parse_argvs() {
    @ARGV > 0 || do {
        print STDERR $USAGE;
	    exit(1);
    };
    $PROJECT = shift @ARGV;
} # 
#line 265 ".//deploy.nw"
sub make_dirs() {
    print STDERR "###\n### Creating Project Subdirectories...\n###\n";
    foreach my $d (@working_dirs) {
        print STDERR "### ... $d\n";
		system("mkdir $d") unless (-e $d && -d _);
	};
    print STDERR "###\n### Project Subdirectories............ DONE\n###\n";
} # make_dirs
#line 276 ".//deploy.nw"
sub new_noweb_doc() {
    my $file = "$PROJECT.nw";
    (-e $file && -f _) && do {
         print STDERR "###\n### Project file \"$file\" does exist...\n".
                      "### EXITING PROGRAM !!!\n";
	     exit(1);   
	};
    print STDERR "###\n### Writing Project NOWEB file: $file\n###\n";
	open(NOWEB,"> $file");
	while (<DATA>) {
        my ($FINDPATH,$FINDPROJECT) = ('@@@PATH@@@','@@@PROJECT@@@');
        my $l = $_;
        $l =~ /$FINDPATH/o && do {
            $l =~ s/$FINDPATH/$PATH/o;
		}; 
		$l =~ /$FINDPROJECT/o && do {
            $l =~ s/$FINDPROJECT/$PROJECT/o;
		};
        print NOWEB $l;
    }; 
	close(NOWEB);
    print STDERR "###\n### NOWEB file........................ DONE\n###\n";
} # new_noweb_doc
#line 302 ".//deploy.nw"
sub extract_files() {
    print STDERR "###\n### Extracting Files from NOWEB file...\n###\n";
    # my $WORK = '$HOME/'.$PATH;
    my $WORK = $PATH;
    my $nwfile = "$PROJECT.nw";
    system << "+++EOS+++" ;
notangle -R\'BASH Environment Variables\' $WORK/$nwfile > $WORK/.bash_VARS ; 
notangle -R\'CSH Environment Variables\'  $WORK/$nwfile > $WORK/.csh_VARS ; 
notangle -Rweaving  $WORK/$nwfile > $WORK/nw2tex ;
notangle -RLaTeXing $WORK/$nwfile > $WORK/ltx ;
chmod a+x $WORK/nw2tex ;
chmod a+x $WORK/ltx ;
ci -l -i0.1 -t-\'\t\t$nwfile: NOWEB file for $PROJECT\' \\
   -m'BASIC TEMPLATE for THIS PROJECT' $nwfile ;
emacs $nwfile \&
$WORK/nw2tex ;
$WORK/ltx ;
/usr/X11R6/bin/ghostview -color -title -magstep -1 \\
                         -portrait -a4 $WORK/docs/$PROJECT.ps \&
+++EOS+++
    print STDERR "###\n### File Extraction................. DONE\n###\n";
} # extract_files
