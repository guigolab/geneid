#!/usr/bin/perl -w 


#=============================================================
# This program convert a file to the FASTA format.
# Usage: FileToFasta file.in > file.out 
#=============================================================

use strict;

use Text::Wrap qw($columns &wrap);

my $IDENTITATE="ID";
my $DEFINITION="DE";
my $SEQUENCE="SQ";
my $ACCESION="AC";

my @line;
my ($var,$varSEQ,$temp)=("","","");

    while (<STDIN>) {
	chomp;
	if ($_=~ /^$IDENTITATE/) {
	    @line = ();
	    @line = split / +/,$_;
	    $var .= "$line[1] ";    
	}#if    
	elsif ($_=~ /^$ACCESION/) {
	    @line = ();
	    @line = split / +/,$_;
	    $line[1]=~ s{;}{}o;
	    $var .= "$line[1] ";
	}#elsif
	elsif ($_=~ /^$DEFINITION/) {
	    ($temp=$_)=~ s{$DEFINITION( )+}{}o;
	    $var .= $temp;
	}#elsif
	elsif ($_=~ /^$SEQUENCE/) {
	    SEQ: while (<STDIN>) {
		last SEQ if m{^//$}o;
		chomp;
		$_ =~ s{[0-9]+|\s+}{}go;
		$varSEQ .= $_ ;
	    }#while
            print ">$var\n";		
	    $columns = 60;
	    print wrap("","",$varSEQ),"\n";
	    ($var,$varSEQ,$temp)=("","","");
	}#elsif    	
    }#while

