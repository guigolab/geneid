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
      }  
      else {
	# Start a new cluster
	$exon_cluster = new ClusterMerge::ExonCluster;
	$exon_cluster->add_Exon($exon,'EXPAND');
	$exon_cluster->strand($exon->strand);
		
	# and add it to the main_cluster feature
	push(@$cluster_list,$exon_cluster);
      }
    }
    $count++;
  }
  return $cluster_list;
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



1;
