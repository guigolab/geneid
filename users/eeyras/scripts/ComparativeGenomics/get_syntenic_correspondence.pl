#!/usr/local/bin/perl -w

use strict;

my $dir = $ARGV[0];

open (IN, "ls -al $dir | ");

my %synteny;
while(<IN>){
    chomp;
    next unless /EN/;
    
    my @i = split;
    
    $i[-1] =~/(EN\S+)_(\d+)/;

    push (@{$synteny{$1}}, $i[-1] );
}

close(IN);

foreach my $key ( keys %synteny ){
    foreach my $s ( @{$synteny{$key}} ){
	print $key."\t".$s."\n";
    }
}
