# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=head1 NAME

ClusterMerge::Transcript - transcript object

=head1 DESCRIPTION

Contains details of coordinates of all exons that make
up a gene transcript.

Creation:

     my $tran = new Transcript();
     my $tran = new Transcript(@exons);

Manipulation:

     # Returns an array of Exon objects
     my @exons = @{$tran->get_all_Exons}     

     # Sorts exons into order (forward for + strand, reverse for - strand)
     $tran->sort()                        

=head1 CONTACT

eae@sanger.ac.uk

=cut


# Let the code begin...

package ClusterMerge::Transcript;
use vars qw(@ISA);
use strict;

use ClusterMerge::Root;
use ClusterMerge::Exon;

@ISA = qw(ClusterMerge::Root);

sub new {
  my($class,@args) = @_;

  if( ref $class ) { 
      $class = ref $class;
  }

  my $self = {};
  bless $self,$class;

  $self->{_trans_exon_array} = [];

  # set stuff in self from @args
  foreach my $a (@args) {
    $self->add_Exon($a);
  }

  return $self; # success - we hope!
}


sub id {
   my $self = shift;
   
   if( @_ ) {
      my $value = shift;
      $self->{_id} = $value;
    }
    return $self->{_id};

}

sub dbID {
   my $self = shift;
   
   if( @_ ) {
      my $value = shift;
      $self->{_dbID} = $value;
    }
    return $self->{_dbID};
}

sub stable_id {
   my $self = shift;
   
   if( @_ ) {
      my $value = shift;
      $self->{_stable_id} = $value;
    }
    return $self->{_stable_id};
}

sub type {
  my ($self, $type) = @_;

  if(defined $type) {
    $self->{'_type'} = $type;
  }

  return $self->{'_type'};
}

=head2 start

 Description: it returns the start coordinate of the lef-most exon, i.e.
              the 5prime exon in the forward strand and the 3prime exon in the reverse strand

=cut


sub start {
  my $self = shift;
  my $arg = shift;
  
  my $strand;
  my $start;
  if( defined $arg ) {
    $self->{'_start'} = $arg;
  } elsif(!  defined $self->{'_start'} ) {

    $strand = $self->start_Exon->strand();
    if( $strand == 1 ) {
      $start = $self->start_Exon->start();
    } else {
      $start = $self->end_Exon->start();
    }
    $self->{'_start'} = $start;
  }
  
  return $self->{'_start'};
}


sub end {
  my $self = shift;
  my $arg = shift;

  my $strand;
  my $end
;
  if( defined $arg ) {
    $self->{'_end'} = $arg;
  } elsif( ! defined $self->{'_end'} ) {
    $strand = $self->start_Exon->strand();
    if( $strand == 1 ) {
      $end = $self->end_Exon->end();
    } else {
      $end = $self->start_Exon->end();
    }
    $self->{'_end'} = $end;
  }
  
  return $self->{'_end'};
}

=head2 add_Exon

 Title   : add_Exon
 Usage   : $trans->add_Exon($exon)
 Returns : Nothing
 Args    :

=cut

sub add_Exon{
   my ($self,$exon) = @_;

   unless(defined $exon && ref $exon && $exon->isa("ClusterMerge::Exon") ) {
     $self->throw("[$exon] is not a ClusterMerge::Exon!");
   }

   #invalidate the start, end and strand - they may need to be recalculated
   $self->{'_start'}  = undef;
   $self->{'_end'}    = undef;
   $self->{'_strand'} = undef;

   push(@{$self->{_trans_exon_array}},$exon);
}


=head2 get_all_Exons

  Arg [1]    : none
  Example    : my @exons = @{$transcript->get_all_Exons()};
  Description: Returns an listref of the exons in this transcipr in order.
               i.e. the first exon in the listref is the 5prime most exon in 
               the transcript.
  Returntype : a list reference to ClusterMerge::Exon objects
  Exceptions : none
=cut
sub get_all_Exons {
    my ($self) = @_;
    return $self->{_trans_exon_array};
}

############################################################
=head2 get_all_Introns

  Arg [1]    : none
  Example    : my @introns = @{$transcript->get_all_Introns()};
  Description: Returns a listref of the introns, where each intron is represented
               by a transcript with the two flanking exons.
=cut

sub get_all_Introns{
    my ($self) = @_;
    my @exons = sort {$a->start <=> $b->start} @{$self->get_all_Exons};
    my @ints;
    my $intron_count = 0;
    for (my $i=0; $i<=scalar(@exons)-2; $i++){
	$intron_count++;
	my $intron = ClusterMerge::Transcript->new();
	$intron->type($self->type);
	my $intron_tag = $self->dbID."-"."intron-".$intron_count;
	$intron->dbID($intron_tag);
	$intron->add_Exon($exons[$i]);
	$intron->add_Exon($exons[$i+1]);
	push (@ints, $intron);
    }
    return \@ints;
}
	  
############################################################

=head2 coding_transcript

This method holds the transcript containing the coding sequence.
It is only a get/set method.

=cut

sub coding_transcript{
    my ($self, $tran) = @_;
    if($tran) {
	$self->{'_coding_transcript'} = $tran;
  }
    
    return $self->{'_coding_transcript'};
}

############################################################

=head2 length

my $t_length = $transcript->length
Returns the sum of the length of all the exons in
the transcript.

=cut

sub length {
    my( $self ) = @_;
    
    my $length = 0;
    foreach my $ex (@{$self->get_all_Exons}) {
        $length += $ex->length;
    }
    return $length;
}

=head2 flush_Exons

 Description : Removes all Exons from the array.


=cut

sub flush_Exons{
   my ($self,@args) = @_;
   $self->{'_exon_coord_mapper'} = undef;
   $self->{'coding_region_start'} = undef;
   $self->{'coding_region_end'} = undef;
   $self->{'_start'} = undef;
   $self->{'_end'} = undef;
   $self->{'_strand'} = undef;

   $self->{_trans_exon_array} = [];
}


=head2 sort

 Usage   : $transcript->sort()
 Function: Sorts the exon 5' to 3'

=cut

sub sort {
    my $self = shift;
    
    # Fetch all the exons
    my @exons = @{$self->get_all_Exons()};
    
    print STDERR "\n";
    foreach my $exon ( @exons ){
	print STDERR "exon: $exon\n";
    }
    

    # Empty the feature table
    $self->flush_Exons();
    
    # Now sort the exons and put back in the feature table
    my $strand = $exons[0]->strand;
    
    if ($strand == 1) {
	@exons = sort { $a->start <=> $b->start } @exons;
    } 
    elsif ($strand == -1) {
	@exons = sort { $b->start <=> $a->start } @exons;
    }
    
    foreach my $e (@exons) {
	$self->add_Exon($e);
    }
}


=head2 start_Exon

 Description: $start_exon = $transcript->start_Exon;
 Returns    : The 5prime exon in the transcript.

=cut

sub start_Exon{
    my ($self) = @_;
    
    my @exons;
   
    # Now sort the exons
    my $strand = $self->get_all_Exons->[0]->strand;
    if ( $strand == 1){
        @exons = sort { $a->start <=> $b->start } @{$self->get_all_Exons};
    } 
    elsif ($strand == -1) {
        @exons = sort { $b->start <=> $a->start } @{$self->get_all_Exons};
    }
    return $exons[0];
}


=head2 end_Exon

 Usage   : $end_exon = $transcript->end_Exon;
 Returns : The 3prime exon in the transcript.

=cut

sub end_Exon{
   my ($self,@args) = @_;
   
   my @exons = @{$self->get_all_Exons()};
   
   # Now sort the exons and put back in the feature table
   my $strand = $exons[0]->strand;
   
   if ($strand == 1) {
       @exons = sort { $a->start <=> $b->start } @exons;
   } 
   elsif ($strand == -1) {
       @exons = sort { $b->start <=> $a->start } @exons;
   }
   return $exons[$#exons];

}

=head2 description

 Title   : description
 Usage   : $obj->description($newval)
 Function: 
 Returns : value of description
 Args    : newvalue (optional)


=cut

sub description{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'description'} = $value;
    }
    return $obj->{'description'};

}

1;
