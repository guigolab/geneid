############################################################
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

ClusterMerge::ExonUtils 

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONTACT

eae@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package ClusterMerge::ExonUtils;

use vars qw(@ISA);
use strict;

use ClusterMerge::Root;
use ClusterMerge::ExonCluster;

@ISA = qw(ClusterMerge::Root);

############################################################

=head2 _cluster_Exons
 
 Function: it cluster exons according to overlap,
=cut

sub _cluster_Exons{
  my ($self, @exons) = @_;

  #print STDERR "EST_GeneBuilder: clustering exons...\n";
  
  # no point if there are no exons!
  return unless ( scalar( @exons) > 0 );   

  # keep track about in which cluster is each exon
  my %exon2cluster;
  
  # main cluster feature - holds all clusters
  my $cluster_list;
  
  # sort exons by start coordinate
  @exons = sort { $a->start <=> $b->start } @exons;

  # Create the first exon_cluster
  my $exon_cluster = new ClusterMerge::ExonCluster;
  
  # Start off the cluster with the first exon
  $exon_cluster->add_Exon($exons[0],'EXPAND');
  $exon2cluster{$exons[0]} = $exon_cluster;
  
  $exon_cluster->strand($exons[0]->strand);    
  push(@$cluster_list, $exon_cluster);
  
  # Loop over the rest of the exons
  my $count = 0;
  
 EXON:
  foreach my $exon (@exons) {
    if ($count > 0) {
      
	# Add to cluster if overlap AND if strand matches
	if ( $exon_cluster->overlaps($exon) && ( $exon->strand == $exon_cluster->strand) ) { 
	    $exon_cluster->add_Exon($exon,'EXPAND');
	    $exon2cluster{$exon} = $exon_cluster;
	}  
	else {
	    # Start a new cluster
	    $exon_cluster = new ClusterMerge::ExonCluster;
	    $exon_cluster->add_Exon($exon,'EXPAND');
	    $exon_cluster->strand($exon->strand);
	    $exon2cluster{$exon} = $exon_cluster;
	    
	    # and add it to the main_cluster feature
	    push(@$cluster_list,$exon_cluster);
	}
    }
    $count++;
}
  return ($cluster_list,\%exon2cluster);
}

############################################################

sub _clone_Exon{
    my ($self,$exon) = @_;
    my $newexon = new ClusterMerge::Exon;
    
    $newexon->seqname    ($exon->seqname);
    $newexon->phase      ($exon->phase);
    $newexon->end_phase  ($exon->end_phase);
    $newexon->start      ($exon->start);
    $newexon->end        ($exon->end);
    $newexon->score      ($exon->score);
    $newexon->strand     ($exon->strand);

    return $newexon;
}

############################################################

sub print_Exons{
    my ($self,$exons) = @_;
    my @exons = sort {$a->start <=> $b->start } @$exons;
    foreach my $e (@exons){
	print $e->start."-".$e->end." ";
    }
    print "\n";
}

############################################################

=head2 eliminate_redundant_lists

  Arg: a list of lists of exons, where there are potentially redundant
    lists of exons that we want to get rif off.
    Note that this is based on an object-based comparison,
    hence you have to make sure that the exon objects are the same
    for exons that ARE the same.
    
    It returns the non redundant lists

=cut

sub eliminate_redundant_lists{
    my ($self,$lists) = @_;
    
    ############################################################
    # sort the lists by the number of exons in descending order
    my @lists = sort { scalar(@$b) <=> scalar(@$a) } @$lists;
    
    my @accepted_lists;
  LIST:
    while ( @lists ){
	my $list = shift @lists;
	my $found_embedding = 0;
	
      ACCEPTED_LIST:
	foreach my $accepted_list ( @accepted_lists ){
	    
	    # use a Boyer-Moore-like method
	    # to check the embedding of one into the other:
	    if ( $self->check_embedding( $list, $accepted_list ) ){
		next LIST;
	    }
	}
	
	# if we get to this point, it means that this list is genuine
	push( @accepted_lists, $list );
    }
    
    return @accepted_lists;
}

############################################################

=head2 check_embedding

     Arg[1]: a listref of ClusterMerge::Exon objects
     Arg[2]: a listref of ClusterMerge::Exon objects which potentially contains Arg[1]
description: this method checks whether Arg[1] is included in Arg[2], in the following sense
             if list2 = ( 1,3,4,6 ) and list1 = ( 3,4) ==> list2 includes list1. 
             It mimics the idea of the Boyer-Moore algorithm

=cut

sub check_embedding{
    my ( $self, $list, $bigger_list ) = @_;
    
    my @list = @$list;
    my @bigger_list = @$bigger_list;
    unless ( scalar(@list) <= scalar( @bigger_list ) ){
	$self->throw("must pass the bigger list as second argument");
    }
    
    # with a bit of preprocessing this can be made pretty quick
    # last element in list:
    my $last_in_list = $list[$#list];
    
    ############################################################
    # store the positions where this element occurs 
    # in bigger_list starting from the right (descending order)
    my @positions;
    for (my $k= $#bigger_list; $k>=0; $k-- ){
	if ( $bigger_list[$k] == $last_in_list ){
	    push( @positions, $k);
	    last;
	}
    }
    
    ############################################################
    # if it does not occur in bigger_list, return false
    unless (@positions){
	return 0;
    }
    
    ############################################################
    # else, start matching from that position:
    my $i = $#list;  
    my $j = shift @positions;
    
    while ( $j >=0 ){
	
	while ( $list[$i] == $bigger_list[$j] && $i >=0 && $j >= 0){
	    $i--;
	    $j--;
	}
	
	if ( $i == -1 ){
	    # list is embedded in bigger_list
	    return 1;
	}
	elsif ( @positions ){
	    # start looking from the next occurrence of $last_in_list in @bigger_list
	    $j = shift @positions;
	    $i = $#list;
	}
	else{
	    # no more occurrences
	    return 0;
	}
    }
    
    # sorry, we got to the end of bigger_list, and no match was found
    return 0;
}

############################################################



1;
