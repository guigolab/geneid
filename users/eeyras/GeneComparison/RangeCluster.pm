package GeneComparison::RangeCluster;
use ClusterMerge::Transcript;
use ClusterMerge::ExonUtils;
use vars qw(@ISA);
use strict;

#@ISA = qw(ClusterMerge::Root);

#########################################################################

sub new{
  my ($class, @args)=@_;

  if (ref($class)){
    $class = ref($class);
  }
  my $self = {};
  bless($self,$class);
  
  return $self;
}
#########################################################################

sub type{
    my ($self, $type ) = @_;
    if ($type){
	$self->{'_type'} = $type;
    }
    return $self->{'_type'};
}

#########################################################################

=head1 Range-like methods

Methods start and end are typical for a range. We also implement the boolean
and geometrical methods for a range.

=head2 start()

  Title   : start
  Usage   : $start = $transcript_cluster->end();
  Function: get/set the start of the range covered by the cluster. This is re-calculated and set everytime
            a new transcript is added to the cluster
  Returns : a number
  Args    : optionaly allows the start to be set

=cut

sub start{
  my ($self,$start) = @_;
  if ($start){
    $self->throw( "$start is not an integer") unless $start =~/^[-+]?\d+$/;
    $self->{'_start'} = $start;
  }
  return $self->{'_start'};
}

############################################################

=head2 end()

  Title   : end
  Usage   : $end = $transcript_cluster->end();
  Function: get/set the end of the range covered by the cluster. This is re-calculated and set everytime
            a new transcript is added to the cluster
  Returns : the end of this range
  Args    : optionaly allows the end to be set
          : using $range->end($end
=cut

sub end{
  my ($self,$end) = @_;
  if ($end){
    $self->throw( "$end is not an integer") unless $end =~/^[-+]?\d+$/;
    $self->{'_end'} = $end;
  }
  return $self->{'_end'};
}

############################################################

=head2 length

  Title   : length
  Usage   : $length = $range->length();
  Function: get/set the length of this range
  Returns : the length of this range
  Args    : optionaly allows the length to be set
          : using $range->length($length)

=cut

sub length{
  my $self = shift @_;
  if (@_){
    $self->confess( ref($self)."->length() is read-only");
  }
  return ( $self->{'_end'} - $self->{'_start'} + 1 );
}

############################################################

=head2 strand

  Title   : strand
  Usage   : $strand = $transcript->strand();
  Function: get/set the strand of the transcripts in the cluster.
            The strand is set in put_Transcripts when the first transcript is added to the cluster, in that
            that method there is also a check for strand consistency everytime a new transcript is added
  Returns : the strandidness (-1, 0, +1)
  Args    : optionaly allows the strand to be set
        
=cut

sub strand{
  my ($self,$strand) = @_;
  if ($strand){
    $self->{'_strand'} = $strand;
  }
  return $self->{'_strand'};
}


############################################################

=head1 Boolean Methods

These methods return true or false. They throw an error if start and end are
not defined. They are implemented in Bio::RangeI.

 $cluster->overlaps($other_cluster) && print "Clusters overlap\n";

=head2 overlaps

  Title   : overlaps
  Usage   : if($cluster1->overlaps($cluster)) { do stuff }
  Function: tests if $cluster2 overlaps $cluster1 overlaps in the sense of genomic-range overlap,
            it does NOT test for exon overlap.
  Args    : arg #1 = a range to compare this one to (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') 
  Returns : true if the clusters overlap, false otherwise

=cut

=head2 contains

  Title   : contains
  Usage   : if($cluster1->contains($cluster2) { do stuff }
  Function: tests whether $cluster1 totally contains $cluster2 
  Args    : arg #1 = a range to compare this one to (mandatory)
	             alternatively, integer scalar to test
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : true if the argument is totaly contained within this range

=cut

=head2 equals

  Title   : equals
  Usage   : if($cluster1->equals($cluster2))
  Function: test whether the range covered by $cluster1 has the same start, end, length as the range 
            covered by $cluster2
  Args    : a range to test for equality
  Returns : true if they are describing the same range

=cut

=head1 Geometrical methods

These methods do things to the geometry of ranges, and return
Bio::RangeI compliant objects or triplets (start, stop, strand) from
which new ranges could be built. They are implemented  here and not in Bio::RangeI, since we
want them to return a new TranscriptCluster object.

=head2 overlap_extent

 Title   : overlap_extent
 Usage   : ($a_unique,$common,$b_unique) = $a->overlap_extent($b)
 Function: Provides actual amount of overlap between two different
           ranges. Implemented already in RangeI
 Example :
 Returns : array of values for 
           - the amount unique to a
           - the amount common to both
           - the amount unique to b
 Args    : 

=cut

#########################################################################

=head2 intersection

  Title   : intersection
  Usage   : $intersection_cluster = $cluster1->intersection($cluster2)
  Function: gives a cluster with the transcripts which fall entirely within the intersecting range of
            $cluster1 and $cluster2
  Args    : arg #1 = a cluster to compare this one to (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') (not implemented here)
  Returns : a TranscriptCluster object, which is empty if the intersection does not contain
            any transcript
=cut

sub intersection{
  my ($self, $cluster) = @_;

  my ($inter_start,$inter_end);
  my ($start1,$end1) = ($self->start   ,   $self->end);
  my ($start2,$end2) = ($cluster->start,$cluster->end);

  my $strand = $cluster->strand;
  
  # if clusters overlap, calculate the intersection
  if ( $self->overlaps( $cluster ) ){
      $inter_start = ($start2 > $start1) ? $start2 : $start1;
      $inter_end = ($end2 < $end1) ? $end2 : $end1;
  }
  return ($inter_start, $inter_end);
}

############################################################

=head2 union

  Title   : union
  Usage   : $union_cluster = $cluster1->union(@clusters);
  Function: returns the union of clusters 
  Args    : a TranscriptCluster or list of TranscriptClusters to find the union of
  Returns : the TranscriptCluster object containing all of the ranges

=cut

sub union{

}
############################################################

=head2 put_Range()

  function to include one or more transcripts in the cluster.
  Useful when creating a cluster. It takes as argument an array of transcripts, it returns nothing.

=cut

sub put_Ranges {
    my ($self, @new_ranges)= @_;
    
    #Get bounds of new transcripts
    my $min_start = undef;
    my $max_end   = undef;
    foreach my $range (@new_ranges) {
	#print "range is a $range\n";
	my ($start, $end) = ($range->start, $range->end);
	if (!defined($min_start) || $start < $min_start) {
	    $min_start = $start;
	}
	if (!defined($max_end) || $end > $max_end) {
	    $max_end = $end;
	}
    }
    
    # check strand consistency
    foreach my $range (@new_ranges){
	unless ( $self->strand ){
	    $self->strand( $range->strand );
	}
	if ( $self->strand != $range->strand ){
	    #$self->warn( "You're trying to put a range [$range] in a cluster of opposite strand");
	}
    }
    
    # if start is not defined, set it
    unless ( $self->start ){
	$self->start( $min_start );
    }
    
    # if end is not defined, set it
    unless ( $self->end ){
	$self->end( $max_end);
    }
    
    # extend start and end if necessary as we include more transcripts
    if ($min_start < $self->start ){
	$self->start( $min_start );
    }
    if ( $max_end > $self->end ){
	$self->end( $max_end );
    }
    
    push ( @{ $self->{'_range_array'} }, @new_ranges );
  
}

#########################################################################

sub get_Ranges {
    my $self = shift @_;
    return $self->{'_range_array'};
}

############################################################

sub print_Cluster{
    my ($self,$type) = @_;
    my @ranges = @{$self->get_Ranges};
    foreach my $t (@ranges){
	print $t->start."-".$t->end." ";
	print "\n";
    }
}

############################################################



1;
