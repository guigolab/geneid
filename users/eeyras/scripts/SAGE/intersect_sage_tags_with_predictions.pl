#!/usr/local/bin/perl -w

use strict;
use ClusterMerge::GFFTools;
use ClusterMerge::Transcript;


# this program reads the sage tags and a set of predictions in GFF format
# and calculates 
# the sensitivity of the predictions with respect to the
# sage set, i.e., how many predictions have a sage tag overlapping/dowstream
# at a given distance,
# and the specificity, i.e. which predictions have at least one sage tag.
#
# The sage tags found are divided into two clases, whether they are downstream
# or whether they are intragenic. Likewise, the intragenic tags are
# divided into exonic (if they overlap an exon) and intronic.
# the exonic ones give also the amount of overlap (<= 21bp of the tag)

if (scalar(@ARGV) < 2) {
    print STDERR "$0 sage_tags.gff predictions.gff\n";
    exit(1);
}

my $sage_file       = $ARGV[0]; 
my $prediction_file = $ARGV[1];

my %counter;
my %annotations;

# read first the sage tags:
#chrE26C13	sage	tag	7338	7358	.	+	.	44845_CATGAGGCAGGGAGGCTGCAG_1
#chrE26C13	sage	tag	37531	37551	.	+	.	66726_CATGGAAACCGCACTGGCCTA_1
#chrE26C13	sage	tag	40569	40589	.	-	.	43387_CATGGGTGTAGGCTGTGTGAT_6

open (ANN,"<$sage_file") or die("cannot open $sage_file");
while (<ANN>){
    chomp;
    my @e = split;      ;
    next unless ($e[3] && $e[4]);
    
    my $t_id    = $e[0]; # chromosome id
    my $s       = $e[3]; # start
    my $e       = $e[4]; # end
    my $id      = $e[8]; # tag identifier
    my $strand  = $e[6]; # strand is important here

    my $uniq_id = $t_id.".".$id.".".$s."-".$e;  # unique id
    
    unless ( $counter{$t_id} ){
	$counter{$t_id} = 0;
    }
    $strand = 1 if $strand eq "+";
    $strand = -1 if $strand eq "-";
    unless ( $annotations{$t_id}{$strand}{$uniq_id} ){
	$annotations{$t_id}{$strand}{$uniq_id} = [$s,$e,$id];
    }
}
close(ANN);

############################################################
# sort the sage tags in each chromosome according to the start coordinate.
my %sorted_ids;
# order the annotations
my $total = 0;
foreach my $t_id ( keys %annotations ){
    #my @sorted = sort { $annotations{$t_id}{$a}[0] <=> $annotations{$t_id}{$b}[0] } keys %{$annotations{$t_id}};
 
    $total += scalar( keys %{$annotations{$t_id}{"1"}} ) + scalar( keys %{$annotations{$t_id}{"-1"}} );
    print scalar( keys %{$annotations{$t_id}{"1"}} )." tags in chromosome $t_id positive strand\n";
    print scalar( keys %{$annotations{$t_id}{"-1"}} )." tags in chromosome $t_id negative strand\n";
    #$sorted_ids{$t_id} = \@sorted;
    
    #foreach my $id ( @sorted ){
    #	print "$t_id @{$annotations{$t_id}{$id}}\n";
    #}
}
print "total number of tags = $total\n";

############################################################
# read the predictions:
open (HSP, "<$prediction_file") or die("cannot open $prediction_file");

# read the other hsps in GFF format
my %tag2transcript;
while (<HSP>) {
    chomp;
    my $exon = ClusterMerge::GFFTools->exon_from_gff ($_);    
    unless ( $tag2transcript{$exon->group_tag} ){
	$tag2transcript{$exon->group_tag} =  ClusterMerge::Transcript->new();
	$tag2transcript{$exon->group_tag}->dbID( $exon->group_tag );
    }
    $tag2transcript{$exon->group_tag}->add_Exon($exon);
}
close(HSP);

############################################################
# scan the predictions one by one
my $total_pred = 0;
foreach my $id ( keys %tag2transcript ){
    my $trans = $tag2transcript{$id}; 
    search_sage_tags( $trans, \%annotations );
    $total_pred++;
}
print "total predictions tested: $total_pred\n";
    
############################################################
# function: search_annotation
# purpose:  searches which annotation intersects with a given begin and end positions by
#           performing binary search (log cost)

sub search_sage_tags{
    my ($trans, $sage) = @_;
    my %sage = %$sage;

    ############################################################
    # get info from this prediction
    my @exons = sort { $a->start <=> $b->start } @{$trans->get_all_Exons};
    my $chr   = $exons[0]->seqname;
    my $strand= $exons[0]->strand;
    my $start = $exons[0]->start;
    my $end   = $exons[-1]->end;

    # $annotations{$t_id}{$strand}{$uniq_id} = [$s,$e,$id];

    ############################################################
    # get the tags for this chromosome and this strand:
    my @tags = values %{$sage{$chr}{$strand}};
    #print "scanning ".scalar(@tags)." tags in chr $chr\n";
    my $closest;
    my $min_distance = 5001;

  TAG:
    foreach my $tag ( @tags ){

	# tags are identified by the chromosome name, the identifier from the mapping (number_sequence_multiplicity, where number is
	# an internal database identifier and multuplicity is the number of times the tag has been found in the genome according to the
	# database creators), and the start and end in the chromosome.

	my $uniq_id   = $chr.".".$tag->[2].".".$tag->[0]."-".$tag->[1];  # unique id
	my $tag_start = $tag->[0];
	my $tag_end   = $tag->[1];

	############################################################
	# forward strand
	if ($strand == +1 ){
	    
	    # only allow tags that are downstream of the first exon
	    next TAG unless ( $tag_start >= $start );

	    # is intragenic?
	    if ( $tag_start <= $end ){
		
		# does it overlap exons?
		my $max_overlap = 0;
		foreach my $e ( @exons ){
		    if ( !( $tag_start > $e->end || $tag_end < $e->start) ){
			my $overlap = ( min($tag_end,$e->end) - max( $tag_start,$e->start) + 1 );
			$max_overlap = $overlap if ( $overlap > $max_overlap );
		    }
		} 
		if ( $max_overlap ){
		    print "INTRAGENIC-EXONIC\t".$trans->dbID."\t".$uniq_id."\t"."overlap: $max_overlap\n";
		}
		else{
		    print "INTRAGENIC-INTRONIC\t".$trans->dbID."\t\t".$uniq_id."\n";
		}
	    }
	    else{
		# is it downstream?
		# calculate distance as the difference of bases between the end of the 3' most exon and the start of the tag
		my $distance = $tag_end - $end;		
		# a tag is considered near if the distance to the downstream exon is smaller than or equal to 5000bp
		print "DOWNSTREAM SENSITIVITY\t".$trans->dbID."\t".$uniq_id."\t"."distance: $distance"."\t".$tag->[2]."\n" if ($distance <= 5000);
		if ( $distance < $min_distance ){
		    $min_distance = $distance;
		    $closest = $tag;
		}
	    }
	}
	############################################################
	# reverse strand         start|                   |end
	#                             ###---###---###---###
	#          tag_start|    tag_end|
	# tag1              xxxxxxxxxxxxx           INTRAGENIC_EXONIC
	#
	#                |-----------|              
	# tag2 xxxxxxxxxxx    distance              DOWNSTREAM
	#
	else{
	    
	    # only allow tags that are downstream of the first exon
	    next TAG unless ( $tag_end <= $end );
	    
	    # is it intragenic?
	    if ( $tag_end >= $start ){
		
		# does it overlap exons?
		my $max_overlap = 0;
		foreach my $e ( @exons ){
		    if ( !( $tag_start > $e->end || $tag_end < $e->start) ){
			my $overlap = ( min($tag_end,$e->end) - max( $tag_start,$e->start) + 1 );
			$max_overlap = $overlap if ( $overlap > $max_overlap );
		    }
		}
		if ( $max_overlap ){
		    print "INTRAGENIC-EXONIC\t".$trans->dbID."\t".$uniq_id."\t"."overlap: $max_overlap\n";
		}
		else{
		    print "INTRAGENIC-INTRONIC\t".$trans->dbID."\t\t".$uniq_id."\n";
		}
	    }
	    else{
		# is it downstream?
		my $distance = $start - $tag_start;
		print "DOWNSTREAM SENSITIVITY\t".$trans->dbID."\t".$uniq_id."\t"."distance: $distance"."\t".$tag->[2]."\n" if ($distance <= 5000);
		if ( $distance < $min_distance ){
		    $min_distance = $distance;
		    $closest = $tag;
		}
	    }
	}	
    } # end of TAG
    
    if ( $closest ){
	my $uniq_id = $chr.".".$closest->[2].".".$closest->[0]."-".$closest->[1];  # unique id
	print "SPECIFICITY\t".$trans->dbID."\t".$uniq_id."\t"."closest_dist: $min_distance"."\t".$closest->[2]."\n";
    }
    
}



############################################################
sub min{
    my ($min, @values) = @_;
    foreach my $v (@values){
        if ($v < $min){
            $min = $v;
        }
    }
    return $min;
}
 
############################################################
sub max{
    my ($max, @values) = @_;
    foreach my $v (@values){
        if ($v > $max){
            $max = $v;
        }
    }
    return $max;
}
 
############################################################

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
