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

package GeneComparison::ExonComparison;

use ClusterMerge::Transcript;
use ClusterMerge::TranscriptCluster;
use ClusterMerge::ExonUtils;
use ClusterMerge::Root;
use vars qw(@ISA);
use strict;

@ISA = qw(ClusterMerge::Root);

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

############################################################

############################################################
# methods related to EXON COMPARISON
############################################################

sub compare_exons_Gene_level{
    my ($self,$predictions,$annotations) = @_;
    
    # each array is a list of transcript clusters (ClusterMerge::TranscriptCluster objects)
    my @predictions = @$predictions;
    my @annotations = @$annotations;
    my ($exon_predictions,
	$exon_annotations,
	$found_exons,
	$wrong_exons,
	$missing_exons) = (0,0,0,0,0);
    
    if ( !@annotations && @predictions){
	my @trans;
	foreach my $transcript_cluster ( @predictions ){
	    push(@trans, @{$transcript_cluster->get_Transcripts});
	}
	my @unique_exons =  $self->get_unique_Exons(\@trans);
	
	($exon_predictions,
	 $exon_annotations,
	 $found_exons,
	 $wrong_exons,
	 $missing_exons) = (scalar(@unique_exons),0,0,scalar(@unique_exons),0);
	
    }
    elsif ( @annotations && !@predictions){
	my @trans;
	foreach my $transcript_cluster ( @annotations ){
	    push(@trans, @{$transcript_cluster->get_Transcripts});
	}
	my @unique_exons = $self->get_unique_Exons(\@trans);
	
	($exon_predictions,
	 $exon_annotations,
	 $found_exons,
	 $wrong_exons,
	 $missing_exons) = (0,scalar(@unique_exons),0,0,scalar(@unique_exons));
    }
    elsif(  @annotations && @predictions ){
	my @trans1;
	#print "eval: ".scalar(@annotations)." annotations and ".scalar(@predictions)." predictions\n";
	foreach my $transcript_cluster ( @predictions ){
	    push(@trans1, @{$transcript_cluster->get_Transcripts});
	}
	my @uniq_exon_pred = $self->get_unique_Exons(\@trans1);
	
	my @trans2;
	foreach my $transcript_cluster ( @annotations ){
	    push(@trans2, @{$transcript_cluster->get_Transcripts});
	}
	my @uniq_exon_ann = $self->get_unique_Exons(\@trans2);
	#print "eval: about to call compare_Exons\n";
	($exon_predictions,
	 $exon_annotations,
	 $found_exons,
	 $wrong_exons,
	 $missing_exons) = $self->compare_Exons(\@uniq_exon_pred, \@uniq_exon_ann);
	
	
    }
    #print  ($exon_predictions,$exon_annotations,$found_exons,$wrong_exons,$missing_exons);

    return ($exon_predictions,
	      $exon_annotations,
	      $found_exons,
	      $wrong_exons,
	      $missing_exons);
}

############################################################
# method to compare 'unique' exons from the two genes
# First it cluster the exons together
# and then it calculates the different
# measures using the clusters:

sub compare_Exons{
    my ($self,$pred_exons,$ann_exons) = @_;
    my ($exon_predictions,
	$exon_annotations,
	$found_exons,
	$wrong_exons,
	$missing_exons) = (0,0,0,0,0);
    
    $exon_predictions = scalar(@$pred_exons);
    $exon_annotations = scalar(@$ann_exons);
    my %type;
    foreach my $e (@$ann_exons){
	$type{$e} = "ann";
    }
    foreach my $e (@$pred_exons){
	$type{$e} = "pred";
    }
    my @exons = (@$pred_exons,@$ann_exons);
    #print "eval: exons = @exons\n";

    ############################################################
    # we cluster by any overlap
    my ($clusters,$exon2cluster) = ClusterMerge::ExonUtils->_cluster_Exons(@exons);
    
  EXON_CLUSTER:
    foreach my $c (@$clusters){
	my @predictions = grep { $type{$_} eq "pred" } @{$c->get_Exons};
	my @annotations = grep { $type{$_} eq "ann"  } @{$c->get_Exons};
	if ( !@annotations && @predictions){
	    $wrong_exons += scalar(@predictions);
	}
	elsif ( @annotations && !@predictions){
	    $missing_exons += scalar(@annotations);
	}
	elsif(  @annotations && @predictions ){
	  EXON_ANN:
	    foreach my $e1 (@annotations){
		foreach my $e2 (@predictions){
		    
		    ############################################################
		    # we impose exact match
		    if ( $e1->start == $e2->start && $e1->end == $e2->end ){
			$found_exons++;
			next EXON_ANN;
		    }
		}
	    }
	}
    }
    #print "compare_Exons: $exon_predictions,$exon_annotations,$found_exons,$wrong_exons,$missing_exons\n";
    return ($exon_predictions,
	    $exon_annotations,
	    $found_exons,
	    $wrong_exons,
	    $missing_exons);
}



############################################################


sub get_unique_Exons {
    my ($self,$transcripts) = @_;
    
    my @unique_exons; 
    
    # keep track of all unique exons found so far to avoid making duplicates
    # need to be very careful about translation->start_Exon and translation->end_Exon
    my @newexons;

    foreach my $tran (@$transcripts){
	#print "transcript: $tran\n";
	foreach my $exon (@{$tran->get_all_Exons}){
	    my $found;
	    
	    #print "exon: ".$exon->start."-".$exon->end."\n";
	    #always empty
	    if (@unique_exons){
	      UNIQ:
		foreach my $uni (@unique_exons) {
		    if ($uni->start  == $exon->start  
			&&
			$uni->end    == $exon->end    
			&&
			$uni->strand == $exon->strand
			) {
			$found = $exon;
			$uni->transcript_tag($uni->transcript_tag.":".$exon->transcript_tag);
			last UNIQ;
		    }
		}
	    }
	    unless (defined($found)) {
		#print "uniq-exon:  $exon:".$exon->start."-".$exon->end."\n";
		push(@unique_exons,$exon);
	    } 
	}          
    }
    #print "eval: returning ".scalar(@unique_exons)." uniq exons\n";
    return @unique_exons;
}
    
############################################################


sub get_unique_Exons_with_coding_flag {
    my ($self,$transcripts) = @_;
    
    my @unique_exons; 
    
    # keep track of all unique exons found so far to avoid making duplicates
    # need to be very careful about translation->start_Exon and translation->end_Exon
    my @newexons;

    my %label;
    foreach my $tran (@$transcripts){
	#print "transcript: $tran\n";
	
	unless ($tran->coding_transcript){
	    print "WARNING ".$tran->dbID." has no coding transcript!\n";
	}
	my @coding_exons = @{$tran->coding_transcript->get_all_Exons};
	
	foreach my $exon (@{$tran->get_all_Exons}){
	    
	    foreach my $codex ( @coding_exons ){
		if ( !( $codex->start > $exon->end || $codex->end < $exon->start ) ){
		    if ( $codex->start == $exon->start && $codex->end == $exon->end ){
			$label{$exon} = "coding";
		    }
		    else{
			$label{$exon} = "half-coding";
		    }
		}
	    }
	    $label{$exon} = "utr" unless $label{$exon};
	    
	    my $found;
	    
	    #print "exon: ".$exon->start."-".$exon->end."\n";
	    #always empty
	    if (@unique_exons){
	      UNIQ:
		foreach my $uni (@unique_exons) {
		    if ($uni->start  == $exon->start  
			&&
			$uni->end    == $exon->end    
			&&
			$uni->strand == $exon->strand
			) {
			$found = $exon;
			$uni->transcript_tag($uni->transcript_tag.":".$exon->transcript_tag);
			
			# we choose the coding to prevail, and the half-coding to prevail over utr.
			if ( $label{$exon} eq "coding" ){
			    $label{$uni} = "coding";
			}
			if( $label{$exon} eq "half-coding" && $label{$uni} eq "utr" ){
			    $label{$uni} = "half-coding";
			}

			last UNIQ;
		    }
		}
	    }
	    unless (defined($found)) {
		#print "uniq-exon:  $exon:".$exon->start."-".$exon->end."\n";
		push(@unique_exons,$exon);
	    } 
	}          
    }
    #print "eval: returning ".scalar(@unique_exons)." uniq exons\n";
    return (\@unique_exons,\%label);
}


############################################################

1;

############################################################
