#!/usr/local/bin/perl

use strict;

while(<>){

    chomp;
    my @list = split '\[';
    
    shift @list;
    
    my $clone = shift @list;
    $clone =~ s/\]//;
    my ($c,$clone_lib) = split ":", $clone;

    my $id = shift @list;
    $id =~ /dbEST\s+ID\:(\d+)/;
    my $dbEST_id = $1;

    foreach my $e ( @list ){
	
	$e =~ /\:([a-zA-Z0-9]+)\]/;
	my $est_id = $1;
	
	print $clone_lib."\t";
	print $dbEST_id."\t";
	print "$est_id\n";
	#print "*******************************\n";
    }


}
