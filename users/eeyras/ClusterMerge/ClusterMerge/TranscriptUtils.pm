############################################################
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

ClusterMerge::TranscriptUtils - 

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONTACT

eae@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package ClusterMerge::TranscriptUtils;

use vars qw(@ISA);
use strict;

use ClusterMerge::Root;
use ClusterMerge::Transcript;
use ClusterMerge::Exon;

@ISA = qw(ClusterMerge::Root);

############################################################
#
# METHODS DOING THE PRINTING
#
############################################################

sub _print_SimpleTranscript{
    my ($self,$transcript,$chr_coord) = @_;
    my @exons = sort { $a->start <=> $b->start } @{$transcript->get_all_Exons};
    
    my $strand = $exons[0]->strand;

    my $id;
    if ($transcript->stable_id){
	$id = $transcript->stable_id;
    }
    elsif ( $transcript->dbID ){
	$id = $transcript->dbID;
    }
    else{
	$id = "no id";
    }
    if ( defined( $transcript->type ) ){
	$id .= " ".$transcript->type;
    }
    print "transcript ".$id.": ".$strand." ";
    

    
    my $shift = 0;
    if ( $chr_coord ){
	$shift = $exons[0]->contig->chr_start - 1;
    }
    foreach my $exon ( @exons){
	print "".($exon->start + $shift)."-".( $exon->end + $shift )." ";
    }
    print "\n";
  }

############################################################


sub _print_Transcript{
  my ($self,$transcript) = @_;
  my @exons = @{$transcript->get_all_Exons};
  my $id;
  if ($transcript->stable_id){
    $id = $transcript->stable_id;
  }
  elsif ( $transcript->dbID ){
    $id = $transcript->dbID;
  }
  else{
    $id = "no id";
  }
  if ( defined( $transcript->type ) ){
    $id .= " ".$transcript->type;
  }
  print STDERR "transcript: ".$id."\n";
  foreach my $exon ( @exons){
    print STDERR $exon->gffstring."\n";
  }
  if ( $transcript->can('translation') && $transcript->translation){
    $self->_print_Translation($transcript);
  }
}

############################################################

sub is_spliced{
  my ($self,$t) = @_;
  my @exons = @{$t->get_all_Exons};
  if ( scalar (@exons ) == 1 ){
    return 0;
  }
  elsif( scalar (@exons) > 1 ){
      
      # check that there are not funky frame shifts
      @exons = sort{ $a->start <=> $b->start } @exons;
      for(my $i=0; $i<$#exons; $i++){
	  my $intron = $exons[$i+1]->start - $exons[$i]->end - 1;
	  if ( $intron > 9 ){
	      return 1;
	  }
      }
      return 0;
  }
  else{
      return 0;
  }
}

############################################################


############################################################

=head2 _difuse_small_introns
 Function: this function is called at the level 5 comparison
           In order to simplfy things, we difuse small introns, according to the
           value of 'intron_mismatch' passe as a parameter, since COMPARISON LEVEL 4 and 5 merges introns of this size.
           inputs merge.
 WARNING:  It does not preserve translations. This was made with ests and cdnas in mind only.
 Returns : a Transcript object. A new one if the introns have been difused or
           the same one we pass in if 'intron_mismatch' is not defined, or there is no introns of this size or smaller.
=cut

sub _difuse_small_introns{
  my ($self,$tran, $intron_mismatch) = @_;

    my $verbose  = 0;

  my $modified = 0;
  
 DIFUSE:
  if ( $intron_mismatch ){
      
      my $newtran = ClusterMerge::Transcript->new();
      if ( $tran->dbID ){
	  $newtran->dbID($tran->dbID);
      }
      my @exons = sort{ $a->start <=> $b->start } @{$tran->get_all_Exons};
      my $exon_count = 0;
      my $current_exon;
      for (my $i=0; $i<=$#exons; $i++){
	  my $exon = ClusterMerge::ExonUtils->_clone_Exon( $exons[$i] );
	  if ( $i>0 ){
	      if ( $exon->start - $current_exon->end - 1 <= $intron_mismatch ){
		  $current_exon->end( $exon->end );
		  $modified++;
	      }
	      else{
		  $current_exon = $exon;
		  $newtran->add_Exon( $current_exon );
	      }
	  }
	  else{
	      $current_exon = $exon;
	      $newtran->add_Exon( $current_exon );
	  }
      }
      if ($modified){
	  if ( $verbose ){
	      print STDERR "difused transcript:\n";
	      print STDERR "before: ";
	      $self->_print_SimpleTranscript($tran);
	      print STDERR "now:    ";
	      $self->_print_SimpleTranscript($newtran);
	  }
	  return $newtran;
      }
      else{
	  return $tran;
      }
  }
  else{
      #$self->warn("LEVEL 4 invoked but no intron_mismatch value defined. Doing level 3 instead");
      return $tran;
  }
}

############################################################

sub transcript_length{
    my ($self, $t) = @_;
    my $length = 0;
    foreach my $e ( @{$t->get_all_Exons} ){
	$length += ($e->end - $e->start + 1 );
    }
    return $length;
}

############################################################

# this method aliminates redundat transcripts from a list
# according to consecutive exon overlap
# it is simplified version of the ClusterMerge idea:
# it returns a list of the transcripts that overlap with the
# rejecte ones and are at least as big as the rejected ones.

sub eliminate_redundant_transcripts{
    my ($self,$trans) =  @_;
    my @trans = sort { 
	my $result = ( scalar(@{$b->get_all_Exons}) <=> scalar(@{$a->get_all_Exons}) );
	if ($result){
	    return $result;
	}
	else{
	    return ( ClusterMerge::TranscriptUtils->transcript_length($b) <=> ClusterMerge::TranscriptUtils->transcript_length($a) );
	}
    } @$trans;
    
    my $comparator = ClusterMerge::TranscriptComparator->new(
							     -comparison_level                 => 5,
							     -intron_mismatch                  => 20,
							     );
    my @accepted_trans;
  TRANS:
    while ( @trans ){
	my $trans = shift @trans;
	my $found_embedding = 0;
	
      ACCEPTED_TRANS:
	foreach my $accepted_trans ( @accepted_trans ){
	    
	    if ( $comparator->compare( $trans, $accepted_trans ) ){
		next TRANS;
	    }
	}
	
	# if we get to this point, it means that this trans is genuine
	push( @accepted_trans, $trans );
    }
    
    return @accepted_trans;
}


############################################################

# it sorts an arrayref of ClusterMerge::Transcript objects
# ascending by their lowered coordinate and descending by their
# higher coordinate. 

# It returns an array of ClusterMerge::Transcript objects

sub sort_transcripts{
    my ($self,$transcripts) = @_;
    my @transcripts = sort { my $result = ( $self->transcript_low($a) <=> $self->transcript_low($b) );
			     if ($result){
				 return $result;
			     }
			     else{
				 return ( $self->transcript_high($b) <=> $self->transcript_high($a) );
			     }
			 } @$transcripts;
    
    return @transcripts;
    
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


1;
