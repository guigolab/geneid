##########################################################################
#                                                                        #
#  This program is free software; you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation; either version 2 of the License, or     #
#  (at your option) any later version.                                   #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program; if not, write to the Free Software           #
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.             #
#                                                                        #
##########################################################################

package GeneComparison::GeneComparison;

use vars qw(@ISA);
use strict;

use GeneComparison::RangeCluster;
use GeneComparison::GeneCluster;
use ClusterMerge::TranscriptCluster;

############################################################

sub new{
    my ($class,@args) = @_;
  
    if (ref($class)){
	$class = ref($class);
    }
    my $self = {};
    bless($self,$class);
    
    return $self;
}

############################################################
# method to pair up clusters of type ClusterMerge::TranscriptCluster
# according to genomic overlap

sub pair_TranscriptClusters{
    my ($self, $pred_clusters, $ann_clusters ) = @_;
    
    my @sets;
    my @pred_clusters = sort { $a->start <=> $b->start } @$pred_clusters;
    my @ann_clusters  = sort { $a->start <=> $b->start } @$ann_clusters;
    
    # do a binary search:
    
    foreach my $ann_cluster ( @ann_clusters ){
	
	my @inter_pred_clusters = binary_search_intersection(@pred_clusters,$ann_cluster);
	if ( @inter_pred_clusters ){
	    
	}
	
    }
}

############################################################
# purpose:  searches which prediction/annotation intersects with a given begin and end positions by
#           performing binary search (log cost)

sub binary_search_intersection {
    my ($self,@clusters, $ann_cluster ) = @_;

    my $s = $ann_cluster->start;
    my $e = $ann_cluster->end;

    my $min    = 0;                       # left bound of the window
    my $max    = scalar(@clusters)-1;   # right bound of the window
    
    my @found;
    my $found  = 0;                       # flag to store whether we've found an intersecting annotation
    my $iannot = -1;
    
    while ($min <= $max && !$found) {
	my $middle = int( ($min+$max)/2 ); # go the the middle of the window
	
	# decide to switch to the left half or the right half of the window
	
	if ($s > ($clusters[$middle]->end ) ){
	    $min = $middle + 1;
	} 
	else {
	    if ($e < ($clusters[$middle]->start ) ){
		$max = $middle - 1;
	    } 
	    else {
		# $middle entry intersects with [$s,$e]
		push (@found, $clusters[$middle]);
		my $i = 1;
		while ( $clusters[$middle+$i] && $self->overlaps($clusters[$middle+$i],$ann_cluster) ){
		    push (@found, $clusters[$middle+$i]);
		    $i++;
		}
		my $j = 1;
		while ( $clusters[$middle-$j] && $self->overlaps($clusters[$middle-$j],$ann_cluster) ){
		    push (@found, $clusters[$middle-$j]);
		    $j++;
		}
		$found = 1;
	    }
	}
    }
    return @found;
}

############################################################

sub overlaps{
    my ($self, $cluster1, $cluster2 ) = @_;
    if ( $cluster1->strand == $cluster2->strand 
	 && 
	 !( $cluster1->end < $cluster2->start || $cluster1->start > $cluster2->end )
	 ){
	return 1;
    }
    else{
	return 0;
    }
}

############################################################
# cluster a set of transcripts according to genomic overlap
#

sub cluster_Transcripts {
    my ($self,$transcripts) = @_;
    
    my $verbose = 0;
    # first sort the transcripts by start and end coordinates
    my @transcripts = sort { my $result = ( $self->transcript_low($a) <=> $self->transcript_low($b) );
			     if ($result){
				 return $result;
			     }
			     else{
				 return ( $self->transcript_high($b) <=> $self->transcript_high($a) );
			     }
			 } @$transcripts;
    
    print "clustering ".scalar(@transcripts)." transcripts\n" if $verbose;
    
    # create a new cluster 
    my $cluster=ClusterMerge::TranscriptCluster->new();
    my $count = 0;
    my @cluster_starts;
    my @cluster_ends;
    my @clusters;
    
    # put the first transcript into these cluster
    $cluster->put_Transcripts( $transcripts[0] );
    #print STDERR "first cluster:\n";
    #ClusterMerge::TranscriptUtils->_print_SimpleTranscript( $transcripts[0] );
    
    $cluster_starts[$count] = $self->transcript_low( $transcripts[0]);
    $cluster_ends[$count]   = $self->transcript_high($transcripts[0]);
    
    # store the list of clusters
    push( @clusters, $cluster );
    
    # loop over the rest of the transcripts
  LOOP1:
    for (my $c=1; $c<=$#transcripts; $c++){
	
	#print STDERR "\nIn cluster ".($count+1)."\n";
	#print STDERR "start: $cluster_starts[$count] end: $cluster_ends[$count]\n";
	#print STDERR "comparing:\n";
	#ClusterMerge::TranscriptUtils->_print_SimpleTranscript( $transcripts[$c] );
	
	if ( !( $self->transcript_high($transcripts[$c]) < $cluster_starts[$count] ||
		$self->transcript_low($transcripts[$c])  > $cluster_ends[$count] ) ){
	    $cluster->put_Transcripts( $transcripts[$c] );
	    
	    # re-adjust size of cluster
	    if ($self->transcript_low($transcripts[$c]) < $cluster_starts[$count]) {
		$cluster_starts[$count] = $self->transcript_low($transcripts[$c]);
	    }
	    if ( $self->transcript_high($transcripts[$c]) > $cluster_ends[$count]) {
		$cluster_ends[$count] =   $self->transcript_high($transcripts[$c]);
	    }
	}
	else{
	    # else, create a new cluster with this feature
	    $count++;
	    $cluster = ClusterMerge::TranscriptCluster->new();
	    $cluster->put_Transcripts( $transcripts[$c] );
	    $cluster_starts[$count] = $self->transcript_low( $transcripts[$c]);
	    $cluster_ends[$count]   = $self->transcript_high($transcripts[$c]);
	    
	    # store it in the list of clusters
	    push(@clusters,$cluster);
	}
    }
    return @clusters;
}

############################################################
# Description: it returns the highest coordinate of a transcript

sub transcript_high{
    my ($self,$tran) = @_;
    my @exons = sort { $a->start <=> $b->start } @{$tran->get_all_Exons};
    return $exons[-1]->end;
}

############################################################
# Description: it returns the lowest coordinate of a transcript

sub transcript_low{
    my ($self,$tran) = @_;
    my @exons = sort { $a->start <=> $b->start } @{$tran->get_all_Exons};
    return $exons[0]->start;
}


############################################################
# cluster a set of transcripts according to genomic overlap
# but separatin by two types

# it returns a list of elements, where each element is a hash with two labels (suplied
# in the method). The hash with each label returns
# a list of clusters associated with each label

sub cluster_TranscriptsClusters_with_type {
    my ($self,$clusters1,$clusters2,$label1,$label2) = @_;

    my $verbose = 1;
    
    my %label;
    foreach my $c ( @$clusters1 ){
	$label{$c} = $label1;
    }
    foreach my $c ( @$clusters2 ){
	$label{$c} = $label2;
    }
    my @ranges = sort { my $result = ( $a->start <=> $b->start );
			if ($result){
			    return $result;
			}
			else{
			    return ( $b->end <=> $a->end );
			}
		    } 	(@$clusters1,@$clusters2);
    
    if ($verbose){
	print "sorted ranges\n";
	foreach my $range ( @ranges ){
	    print $range->start."-".$range->end."\n";
	}
    }


    my $count = 0;
    my @cluster_starts;
    my @cluster_ends;
    my @clusters;
    
    # put the first range into these cluster
    push( @{$clusters[0]{$label{$ranges[0]}}}, $ranges[0] );
    $cluster_starts[$count] = $ranges[0]->start;
    $cluster_ends[$count]   = $ranges[0]->end;
    
    # loop over the rest of the transcripts
  LOOP:
    for (my $c=1; $c<=$#ranges; $c++){
	
	print "current range: ".$ranges[$c]->start." - ".$ranges[$c]->end."\n" if $verbose;
	print "current cluster $cluster_starts[$count] - $cluster_ends[$count]\n" if $verbose;
	
	if ( !( $ranges[$c]->end  < $cluster_starts[$count] ||
		$ranges[$c]->start  > $cluster_ends[$count] ) ){
	    print "accept\n" if $verbose;
	    push( @{$clusters[0]{$label{$ranges[$c]}}}, $ranges[$c] );
	    
	    # re-adjust size of cluster
	    if ($ranges[$c]->start < $cluster_starts[$count]) {
		$cluster_starts[$count] = $ranges[$c]->start;
	    }
            if ( $ranges[$c]->end > $cluster_ends[$count]) {
		$cluster_ends[$count] =   $ranges[$c];
	    }
	}
	else{
	    print "make new cluster\n";
	    # else, create a new cluster with this feature
	    $count++;
	    push( @{$clusters[$count]{$label{$ranges[$c]}}}, $ranges[$c] );
	    $cluster_starts[$count] = $ranges[$c]->start;
	    $cluster_ends[$count]   = $ranges[$c]->end;
	}
    }
    return @clusters;
}


############################################################


### needs some fixing!!!!

sub cluster_Transcripts_with_type{
    my ($self,$trans1,$trans2,$label1,$label2) = @_;

    my $verbose = 1;
    
    my %label;
    foreach my $c ( @$trans1 ){
	$label{$c} = $label1;
    }
    foreach my $c ( @$trans2 ){
	$label{$c} = $label2;
    }
    
    my @transcripts = sort { my $result = ( $self->transcript_low($a) <=> $self->transcript_low($b) );
			if ($result){
			    return $result;
			}
			else{
			    return ( $self->transcript_high($b) <=> $self->transcript_high($a) );
			}
		    } 	(@$trans1,@$trans2);
    
    if ($verbose){
	print "sorted transcripts\n";
	foreach my $range ( @transcripts ){
	    print $range->start."-".$range->end."\n";
	}
    }


    my $count = 0;
    my @cluster_starts;
    my @cluster_ends;
    my @trans;
    
    # put the first range into these cluster
    push( @{$trans[0]{$label{$transcripts[0]}}}, $transcripts[0] );
    $cluster_starts[$count] = $transcripts[0]->start;
    $cluster_ends[$count]   = $transcripts[0]->end;
    
    # loop over the rest of the transcripts
  LOOP:
    for (my $c=1; $c<=$#transcripts; $c++){
	
	print "current range: ".$transcripts[$c]->start." - ".$transcripts[$c]->end."\n" if $verbose;
	print "current cluster $cluster_starts[$count] - $cluster_ends[$count]\n" if $verbose;
	
	if ( !( $transcripts[$c]->end  < $cluster_starts[$count] ||
		$transcripts[$c]->start  > $cluster_ends[$count] ) ){
	    print "accept\n" if $verbose;
	    push( @{$trans[0]{$label{$transcripts[$c]}}}, $transcripts[$c] );
	    
	    # re-adjust size of cluster
	    if ($transcripts[$c]->start < $cluster_starts[$count]) {
		$cluster_starts[$count] = $transcripts[$c]->start;
	    }
            if ( $transcripts[$c]->end > $cluster_ends[$count]) {
		$cluster_ends[$count] =   $transcripts[$c];
	    }
	}
	else{
	    print "make new cluster\n";
	    # else, create a new cluster with this feature
	    $count++;
	    push( @{$trans[$count]{$label{$transcripts[$c]}}}, $transcripts[$c] );
	    $cluster_starts[$count] = $transcripts[$c]->start;
	    $cluster_ends[$count]   = $transcripts[$c]->end;
	}
    }
    return @trans;
}

############################################################



############################################################
# cluster a set of TranscriptClusters according to genomic overlap
#

sub cluster_TranscriptClusters {
    my ($self,$trans_clusters) = @_;
    
    my $verbose = 0;
    # first sort the transcripts by start and end coordinates
    my @trans_clusters = sort { my $result = ( $a->start <=> $b->start );
				if ($result){
				    return $result;
				}
				else{
				    return ( $b->end <=> $a->end );
				}
			    } @$trans_clusters;
    
    print "GeneComparison:clustering ".scalar(@trans_clusters)." trans clusters\n" if $verbose;
    # create a new cluster 
    my $cluster = GeneComparison::RangeCluster->new();
    my $count = 0;
    my @cluster_starts;
    my @cluster_ends;
    my @clusters;
    
    # put the first TranscriptCluster into this cluster
    $cluster->put_Ranges( $trans_clusters[0] );
     
    $cluster_starts[$count] = $trans_clusters[0]->start;
    $cluster_ends[$count]   = $trans_clusters[0]->end;
    
    # store the list of clusters
    push( @clusters, $cluster );
    
    # loop over the rest of the trans_clusters
  LOOP1:
    for (my $c=1; $c<=$#trans_clusters; $c++){
	
	if ( !( $trans_clusters[$c]->end < $cluster_starts[$count] ||
		$trans_clusters[$c]->start  > $cluster_ends[$count] ) ){
	    $cluster->put_Ranges( $trans_clusters[$c] );
	    
	    # re-adjust size of cluster
	    if ( $trans_clusters[$c]->start < $cluster_starts[$count] ) {
		$cluster_starts[$count] = $trans_clusters[$c]->start;
	    }
	    if ( $trans_clusters[$c]->end > $cluster_ends[$count]) {
		$cluster_ends[$count] =   $trans_clusters[$c]->end;
	    }
	}
	else{
	    # else, create a new cluster with this feature
	    $count++;
	    $cluster = GeneComparison::RangeCluster->new();
	    $cluster->put_Ranges( $trans_clusters[$c] );
	    $cluster_starts[$count] = $trans_clusters[$c]->start;
	    $cluster_ends[$count]   = $trans_clusters[$c]->end;
	    
	    # store it in the list of clusters
	    push(@clusters,$cluster);
	}
    }
    return @clusters;
}

############################################################

1;
