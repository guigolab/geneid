# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

ClusterMerge::ExonCluster

=head1 CONTACT

eae@sanger.ac.uk

=cut


# Let the code begin...

package ClusterMerge::ExonCluster;

use vars qw(@ISA);
use strict;

use ClusterMerge::Root;
@ISA = qw(ClusterMerge::Root);


sub new {
  my ($class) = @_;

  if (ref($class)){
    $class = ref($class);
  }
  my $self = {};
  bless($self,$class);
  
  return $self;
}
############################################################

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

=head2 score

 Title   : score
 Usage   : $score = $feat->score()
           $feat->score($score)
 Function: get/set on score information
 Returns : float
 Args    : none if get, the new value if set


=cut

sub score {
    my ($self,$value) = @_;

    if(defined ($value) ) {
      if( $value !~ /^[+-]?\d+\.?\d*(e-\d+)?/ ) {
          $self->throw("'$value' is not a valid score");
      }
      $self->{'_gsf_score'} = $value;
  }

  return $self->{'_gsf_score'};
}

=head2 get_Exons

=cut

sub get_Exons{
  my ($self) = @_;
  
  if($self->{'_gsf_sub_array'}){
    return $self->{'_gsf_sub_array'};
  }
  else{
    return;
  }
}

=head2 add_Exon

=cut

sub add_Exon{
   my ($self,$feat,$expand) = @_;

   if( $expand eq 'EXPAND' ) {
       # if this doesn't have start/end set - forget it!
       if( !defined $self->start && !defined $self->end ) {
           $self->start($feat->start());
           $self->end($feat->end());
           $self->strand($feat->strand);
       } else {
           my ($start,$end);
           if( $feat->start < $self->start ) {
               $start = $feat->start;
           }

           if( $feat->end > $self->end ) {
               $end = $feat->end;
           }

           $self->start($start);
           $self->end($end);

       }
   } else {
       if( !$self->contains($feat) ) {
           $self->throw("$feat is not contained within parent feature, and expansion is not valid");
       }
   }

   push(@{$self->{'_gsf_sub_array'}},$feat);

}

=head2 flush_sub_SeqFeature

 Title   : flush_sub_SeqFeature
 Usage   : $sf->flush_sub_SeqFeature
 Function: Removes all sub SeqFeature
           (if you want to remove only a subset, take
            an array of them all, flush them, and add
            back only the guys you want)
 Example :
 Returns : none
 Args    : none


=cut

sub flush_sub_SeqFeature {
   my ($self) = @_;

   $self->{'_gsf_sub_array'} = []; # zap the array implicitly.
}


sub id {
    my ($self,$value) = @_;

    if (defined($value)) {
        $self->{_id} = $value;
    }

    return $self->{_id};

}


sub gffstring {
   my ($self) = @_;

   my $str;

   my $strand = "+";
   
   if ((defined $self->strand)&&($self->strand == -1)) {
     $strand = "-";
   }
   
   $str .= (defined $self->seqname)     ?   $self->seqname."\t"      :  "\t";
   $str .= (defined $self->source_tag)  ?   $self->source_tag."\t"   :  "\t";
   $str .= (defined $self->primary_tag) ?   $self->primary_tag."\t"  :  "\t";
   $str .= (defined $self->start)       ?   $self->start."\t"        :  "\t";
   $str .= (defined $self->end)         ?   $self->end."\t"          :  "\t";
   $str .= (defined $self->score)       ?   $self->score."\t"        :  "\t";
   $str .= (defined $self->strand)      ?   $strand."\t"             :  ".\t";
   $str .= (defined $self->phase)       ?   $self->phase."\t"        :  ".\t";
   eval{
     $str .= (defined $self->end_phase)   ?   $self->end_phase."\t"        :  ".\t";
   };

   return $str;
}


sub overlaps{
    my ( $self, $range ) = @_;
    if ( $self->strand == $range->strand &&
	 !( $self->start > $range->end || $self->end < $range->start ) ){
	return 1;
    }
    return 0;
}


1;

