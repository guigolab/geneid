#!/usr/local/bin/perl -w

use strict;
use ClusterMerge::Exon;
use ClusterMerge::ExonUtils;
use ClusterMerge::GFFTools;

# this program reads in one shot all the annotations in gff format contained in the
# file given as a first argument, and intersects them with the annotations given in another GFF format file
# as a second parameter. the intersection is given by the standard output STDOUT

if (scalar(@ARGV) < 2) {
    print "$0 hsps1.gff hsps2.gff\n";
    print "Description: script to obtain the intersections of exons/HSPs purely on coordinate and strand level\n";
    print "             It reports the HSPs from hsps1.gff that intersects with HSPs in hsps2.gff\n";
    exit(1);
}

my $pred_file = $ARGV[0];
my $ann_file = $ARGV[1]; 

############################################################

my %forward_hsps;
my %reverse_hsps;

my %tag;
my $ann_type;
open (ANN,"<$ann_file") or die("cannot open $ann_file");
while (<ANN>){
    #chr1	SGP_v1.0	Internal	32133	32262	31.70	-	0	chr1_1
    
    chomp;
    my @e = split;      
    next unless ($e[3] && $e[4]);
    
    my $hsp = ClusterMerge::GFFTools->exon_from_gff ($_);

    if( $e[6] eq '+') {
	push ( @{$forward_hsps{$e[0]}}, $hsp );
    } 
    elsif ( $e[6] eq '-' ){
	push ( @{$reverse_hsps{$e[0]}}, $hsp );
    }
    else{
	push ( @{$forward_hsps{$e[0]}}, $hsp );
	push ( @{$reverse_hsps{$e[0]}}, $hsp );
    }
    unless ($ann_type){
	$ann_type = $e[1];
    }
}
close(ANN);

############################################################

my $pred_type;
open (PRED,"<$pred_file") or die("cannot open $pred_file");
while (<PRED>){
    chomp;
    my @e = split;      
    next unless ($e[3] && $e[4]);
    my $hsp = ClusterMerge::GFFTools->exon_from_gff($_);
    if( $e[6] eq '+') {
	push ( @{$reverse_hsps{$e[0]}}, $hsp );
    } 
    elsif ( $e[6] eq '-' ){
	push ( @{$forward_hsps{$e[0]}}, $hsp );
    }
    else{
	push ( @{$reverse_hsps{$e[0]}}, $hsp );
    	push ( @{$forward_hsps{$e[0]}}, $hsp );
    }
    unless ( $pred_type ){
	$pred_type = $e[1];
    }
}
close(PRED);

############################################################

#print "anntype = $ann_type predtype = $pred_type\n";

foreach my $tid ( keys %forward_hsps ){
    
    my ($clusters,$exon2cluster)= ClusterMerge::ExonUtils->_cluster_Exons(@{$forward_hsps{$tid}} );
    
    #print "FWD: ".scalar(@$clusters)." clusters found for $tid\n";
    foreach my $c ( @$clusters ){
	my @ann  = grep { $_->source_tag eq $ann_type }  @{$c->get_Exons};
	my @pred = grep { $_->source_tag eq $pred_type } @{$c->get_Exons};
	
	if ( @ann && @pred ){
	    foreach my $hsp (@pred){
		print $hsp->gff_string."\n";
	    }
	}
    }
    delete $forward_hsps{$tid};
}

foreach my $tid ( keys %reverse_hsps ){
    
    my ($clusters,$exon2cluster)= ClusterMerge::ExonUtils->_cluster_Exons(@{$reverse_hsps{$tid}} );

    #print "REV: ".scalar(@$clusters)." clusters found for $tid\n";
    foreach my $c ( @$clusters ){
	my @ann  = grep { $_->source_tag eq $ann_type }  @{$c->get_Exons};
	my @pred = grep { $_->source_tag eq $pred_type } @{$c->get_Exons};
	
	if ( @ann && @pred ){
	    foreach my $hsp (@pred){
		print $hsp->gff_string."\n";
	    }
	}
    }
    delete $reverse_hsps{$tid};
}


############################################################
