#!/usr/local/bin/perl -w

use strict;

# this program reads in one shot all the annotations in gff format contained in the
# file given as a first argument, and intersects them with the annotations given in another GFF format file
# as a second parameter. the intersection is given by the standard output STDOUT

if (scalar(@ARGV) < 2) {
    print "intersectannotations.pl hsps.gff annotations.gff\n";
    exit(1);
}

my $hsp_file = $ARGV[0];
my $ann_file = $ARGV[1]; 

my %counter;
my %annotations;

open (ANN,"<$ann_file") or die("cannot open $ann_file");
while (<ANN>){
    chomp;
    my @e = split;
    next unless @e;

    my $t_id = $e[0];
    my $s    = $e[3];
    my $e    = $e[4];
    my $id   = $e[8];
    
    unless ( $counter{$t_id} ){
	$counter{$t_id} = 0;
    }
    unless ( $annotations{$t_id}{$id} ){
	$annotations{$t_id}{$id} = [$s,$e,$id];
    }
    if ( $annotations{$t_id}{$id} ){
	if ( $s < $annotations{$t_id}{$id}[0] ){
	    $annotations{$t_id}{$id}[0] = $s;
	}
	if ( $e > $annotations{$t_id}{$id}[1] ){
	    $annotations{$t_id}{$id}[1] = $e;
	}
    }
    
}
close(ANN);

my %sorted_ids;
# order the annotations
foreach my $t_id ( keys %annotations ){
    my @sorted = sort { $annotations{$t_id}{$a}[0] <=> $annotations{$t_id}{$b}[0] } keys %{$annotations{$t_id}};
    $sorted_ids{$t_id} = \@sorted;
    #foreach my $id ( @sorted ){
#	print "$t_id @{$annotations{$t_id}{$id}}\n";
#    }
}

open (HSP, "<$hsp_file") or die("cannot open $hsp_file");

# read the other hsps in GFF format
while (<HSP>) {
    chomp;
    my @a = split;
    
    my $chr    = $a[0];
    my $source = $a[1];
    my $feature= $a[2];
    my $s      = $a[3];
    my $e      = $a[4];
    my $score  = $a[5];
    my $strand = $a[6];
    my $id     = $a[8];
    
    my $i = search_annotation($chr,$s,$e);
    
    if ($i >= 0){
	my $string = join "\t",@a;
	print $string."\n";
	#print "\tOVERLAPS @{$annotations{$chr}{$sorted_ids{$chr}[$i]}}\n";
    }
    
    #foreach my $id (@{$sorted_ids{$chr}}){
    #	if ( !( $s > $annotations{$chr}{$id}[1] || $e < $annotations{$chr}{$id}[0] ) ){
    #	    my $string = join "\t",@a;
    #	    print "DOUBLECHECKL ".$string."\tOVERLAPS @{$annotations{$chr}{$id}}\n";
    #	}
    #}
}
    

close(HSP);


############################################################
# function: search_annotation
# purpose:  searches which annotation intersects with a given begin and end positions by
#           performing binary search (log cost)

sub search_annotation {
    my ($chr, $s, $e ) = @_;

    unless ( $sorted_ids{$chr} ){
	return -1;
    }
    my @sorted_ids = @{$sorted_ids{$chr}};
    my $min    = 0;                       # left bound of the window
    my $max    = scalar(@sorted_ids)-1;   # right bound of the window
    my $found  = 0;                       # flag to store whether we've found an intersecting annotation
    my $iannot = -1;
    
    my $flanking = 1000;

    while ($min <= $max && !$found) {
	my $middle = int( ($min+$max)/2 ); # go the the middle of the window
	
	# decide to switch to the left half or the right half of the window
	
	if ($s > ($annotations{$chr}{$sorted_ids[$middle]}[1] + $flanking) ) {
	    $min = $middle + 1;
	} 
	else {
	    if ($e < ($annotations{$chr}{$sorted_ids[$middle]}[0] - $flanking) ) {
		$max = $middle - 1;
	    } 
	    else {
		$found = 1;
		$iannot = $middle;
	    }
	}
    }
    return $iannot;
}
