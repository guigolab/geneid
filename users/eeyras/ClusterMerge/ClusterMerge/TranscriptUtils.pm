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
    print STDERR "transcript ".$id.": ";
    

    
    my $shift = 0;
    if ( $chr_coord ){
      $shift = $exons[0]->contig->chr_start - 1;
    }
    foreach my $exon ( @exons){
      print STDERR ($exon->start + $shift)."-".( $exon->end + $shift )." ";
    }
    print STDERR "\n";
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
  my $modified = 0;
  my $verbose  = 1;
  
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
1;
