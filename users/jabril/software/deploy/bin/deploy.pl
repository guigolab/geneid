#!/usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
#
#line 277 "/home/ug/jabril/development/softjabril/deploy/deploy.nw"
# deploy.pl [options] project_name template_file
#
#     Preparing all files we need to work with noweb on a project.
#line 1509 "/home/ug/jabril/development/softjabril/deploy/deploy.nw"
#
#line 1679 "/home/ug/jabril/development/softjabril/deploy/deploy.nw"
# #----------------------------------------------------------------#
# #                             DEPLOY                             #
# #----------------------------------------------------------------#
# 
# Creates basic file set to work with noweb literate programming tool.
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
# #----------------------------------------------------------------#
#line 1511 "/home/ug/jabril/development/softjabril/deploy/deploy.nw"
#
#line 1673 "/home/ug/jabril/development/softjabril/deploy/deploy.nw"
# $Id: deploy.pl,v 1.5 2001-10-03 11:21:02 jabril Exp $
#line 1513 "/home/ug/jabril/development/softjabril/deploy/deploy.nw"
#
use strict;
#
#line 272 "/home/ug/jabril/development/softjabril/deploy/deploy.nw"
my $PROGRAM = 'deploy.pl';
my $VERSION = '1.0_alpha';
#line 1517 "/home/ug/jabril/development/softjabril/deploy/deploy.nw"
my $DATE = localtime;
my $USER = defined($ENV{USER}) ? $ENV{USER} : 'Child Process';
my $host = `hostname`;
chomp($host);
#
#line 250 "/home/ug/jabril/development/softjabril/deploy/deploy.nw"
#
# MODULES
#

#
# VARIABLES
#
#line 283 "/home/ug/jabril/development/softjabril/deploy/deploy.nw"
my $USAGE = "\nUSAGE:\n\tdeploy.pl <projectname> <template>\n".
            "(It asumes that you are in the right directory)\n\n";
my @working_dirs = qw(
                       bin  bin/param
                       data
                       docs docs/psfigures docs/tables docs/html
                       tests
                       );
my ($PROJECT, $TEMPLATE);
# my $HOME = $ENV{HOME};
my $CWD  = `pwd`;
chomp($CWD);
my $PATH = $CWD;
# $PATH =~ s%^$HOME/%%o;
#line 258 "/home/ug/jabril/development/softjabril/deploy/deploy.nw"
#
# MAIN LOOP
#
#line 300 "/home/ug/jabril/development/softjabril/deploy/deploy.nw"
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
#line 262 "/home/ug/jabril/development/softjabril/deploy/deploy.nw"
#
# FUNCTIONS
#
#line 319 "/home/ug/jabril/development/softjabril/deploy/deploy.nw"
sub parse_argvs() {
    scalar(@ARGV) == 2 || do {
        print STDERR $USAGE;
	    exit(1);
    };
    $PROJECT = shift @ARGV;
    $TEMPLATE = shift @ARGV;
} # 
#line 330 "/home/ug/jabril/development/softjabril/deploy/deploy.nw"
sub make_dirs() {
    print STDERR "###\n### Creating Project Subdirectories...\n###\n";
    foreach my $d (@working_dirs) {
        print STDERR "### ... $d\n";
		system("mkdir $d") unless (-e $d && -d _);
	};
    print STDERR "###\n### Project Subdirectories............ DONE\n###\n";
} # make_dirs
#line 341 "/home/ug/jabril/development/softjabril/deploy/deploy.nw"
sub new_noweb_doc() {
    my $file = "$PROJECT.nw";
    (-e $file && -f _) && do {
         print STDERR "###\n### Project file \"$file\" does exist...\n".
                      "### EXITING PROGRAM !!!\n";
	     exit(1);   
	};
    print STDERR "###\n### Writing Project NOWEB file: $file\n###\n";
	open(NOWEB,"> $file");
	open(DATA,"< $TEMPLATE") ||
        die ("#### ERROR #### Template File does not exists: $TEMPLATE . $!\n");
    (my $FIXPROJ = $PROJECT) =~ s{([\W_])}{\\$1}og;
	while (<DATA>) {
        my ($FINDPATH,$FINDPROJECT,$FINDQUOTED) =
           ('@@@PATH@@@','@@@PROJECT@@@','@@@QUOTED@@@');
        my $l = $_;
      CURRENT: {
          $l =~ /$FINDPATH/o && do {
              $l =~ s/$FINDPATH/$PATH/o;
              last CURRENT;
          }; 
		  $l =~ /$FINDPROJECT/o && do {
              $l =~ s/$FINDPROJECT/$PROJECT/o;
              last CURRENT;
          };
          $l =~ /$FINDQUOTED/o && do {
              $l =~ s/$FINDQUOTED/$FIXPROJ/o;
		  };
	  }; # CURRENT
        print NOWEB $l;
    }; 
	close(DATA);
	close(NOWEB);
    print STDERR "###\n### NOWEB file........................ DONE\n###\n";
} # new_noweb_doc
#line 379 "/home/ug/jabril/development/softjabril/deploy/deploy.nw"
sub extract_files() {
    print STDERR "###\n### Extracting Files from NOWEB file...\n###\n";
    # my $WORK = '$HOME/'.$PATH;
    my $WORK = $PATH;
    my $nwfile = "$PROJECT.nw";
    system << "+++EOS+++" ;
notangle -R\'BASH Environment Variables\' $WORK/$nwfile > $WORK/.bash_VARS ; 
notangle -Rweaving  $WORK/$nwfile > $WORK/nw2tex ;
notangle -RLaTeXing $WORK/$nwfile > $WORK/ltx ;
chmod a+x $WORK/nw2tex ;
chmod a+x $WORK/ltx ;
# ci -l -i0.1 -t-\'\t\t$nwfile: NOWEB file for $PROJECT\' \\
#    -m'BASIC TEMPLATE for THIS PROJECT' $nwfile ;
emacs $nwfile \&
$WORK/nw2tex ;
$WORK/ltx ;
/usr/X11R6/bin/ghostview -color -title -magstep -1 \\
                         -portrait -a4 $WORK/docs/$PROJECT.ps \&
+++EOS+++
    print STDERR "###\n### File Extraction................. DONE\n###\n";
} # extract_files
