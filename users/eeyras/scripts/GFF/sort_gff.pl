#!/usr/local/bin/perl -w

use strict;
use Getopt::Long;


############################################################
# script to conver split a gff into different files according to 
# the target id

my %targets;
while (<STDIN>){
    chomp;
    my @e = split;
    next unless @e;
    push( @{$targets{$e[0]}}, \@e );
}

my $t_id;
foreach my $target ( sort { $a cmp $b } keys %targets ){
    foreach my $entry ( sort { $a->[3] <=> $b->[3] } @{$targets{$target}}  ){
	
	############################################################
	# evaluation separator
	if ( defined ($t_id) && !($t_id eq $entry->[0])){
	    print "#\$\n";
	}

	my $line = join "\t",@$entry;
	print $line."\n";
	$t_id = $entry->[0];
    }
}
