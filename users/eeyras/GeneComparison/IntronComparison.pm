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
 
package GeneComparison::IntronComparison;

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
# METHODS FOR THE COMPARISON OF INTRONS
############################################################
# and intron is represented as a transcript with two exons

sub compare_Introns{
    my ($self,$pred_introns,$ann_introns) = @_;
    # $pred_introns and $ann_introns are two arrayrefs of introns
    
    my ($intron_predictions,
	$intron_annotations,
	$found_introns,
	$wrong_introns,
	$missing_introns) = (0,0,0,0,0);
    
    $intron_predictions = scalar(@$pred_introns);
    $intron_annotations = scalar(@$ann_introns);
    my %type;
    foreach my $e (@$ann_introns){
	$type{$e} = "ann";
    }
    foreach my $e (@$pred_introns){
	$type{$e} = "pred";
    }
    my @introns = (@$pred_introns,@$ann_introns);

    # we cluster introns by simple intron overlap
    my ($clusters,$intron2cluster) = $self->cluster_Introns(\@introns);
    
  EXON_CLUSTER:
    foreach my $c (@$clusters){
	my @predictions = grep { $type{$_} eq "pred" } @{$c->get_Transcripts};
	my @annotations = grep { $type{$_} eq "ann"  } @{$c->get_Transcripts};
	if ( !@annotations && @predictions){
	    $wrong_introns += scalar(@predictions);
	}
	elsif ( @annotations && !@predictions){
	    $missing_introns += scalar(@annotations);
	}
	elsif(  @annotations && @predictions ){
	  INTRON_ANN:
	    foreach my $i1 (@annotations){
		foreach my $i2 (@predictions){
		    
		    if ( $self->compare_2introns($i1,$i2,'exact') ){
			$found_introns++;
			next INTRON_ANN;
		    }
		}
	    }
	}
    }
    return ($intron_predictions,
	    $intron_annotations,
	    $found_introns,
	    $wrong_introns,
	    $missing_introns);
}





############################################################

sub compare_introns_Gene_level{
    my ($self,$predictions,$annotations) = @_;
    
    # each array is a list of transcript clusters (ClusterMerge::TranscriptCluster objects)
    my @predictions = @$predictions;
    my @annotations = @$annotations;
    my ($intron_predictions,
	$intron_annotations,
	$found_introns,
	$wrong_introns,
	$missing_introns) = (0,0,0,0,0);
    
    if ( !@annotations && @predictions){
	if ( scalar(@predictions) >= 1 ){
	    my @trans;
	    foreach my $transcript_cluster ( @predictions ){
		push(@trans, @{$transcript_cluster->get_Transcripts});
	    }
	    my @introns = $self->get_uniq_introns(@trans);
	    ($intron_predictions,
	     $intron_annotations,
	     $found_introns,
	     $wrong_introns,
	     $missing_introns) = (scalar(@introns),0,0,scalar(@introns),0);
	}
    }
    elsif ( @annotations && !@predictions){
	if ( scalar(@annotations) >= 1 ){
	    my @trans;
	    foreach my $transcript_cluster ( @annotations ){
		push(@trans, @{$transcript_cluster->get_Transcripts});
	    }
	    my @introns = $self->get_uniq_introns(@trans);
	    ($intron_predictions,
	     $intron_annotations,
	     $found_introns,
	     $wrong_introns,
	     $missing_introns) = (0,scalar(@introns),0,0,scalar(@introns));
	}
    }
    elsif(  @annotations && @predictions ){
	my @trans1;
	#print "eval: ".scalar(@annotations)." annotations and ".scalar(@predictions)." predictions\n";
	foreach my $transcript_cluster ( @predictions ){
	    push(@trans1, @{$transcript_cluster->get_Transcripts});
	}
	my @introns1 = $self->get_uniq_introns(@trans1);
	my @trans2;
	foreach my $transcript_cluster ( @annotations ){
	    push(@trans2, @{$transcript_cluster->get_Transcripts});
	}
	my @introns2 = $self->get_uniq_introns(@trans2);
	#print "eval: about to call compare_Exons\n";
	($intron_predictions,
	 $intron_annotations,
	 $found_introns,
	 $wrong_introns,
	 $missing_introns) = $self->compare_Introns(\@introns1,\@introns2);
	
	
    }
    
    return ($intron_predictions,
	    $intron_annotations,
	    $found_introns,
	    $wrong_introns,
	    $missing_introns);



}

############################################################
# method to cluster introns
#
# since each intron is in fact a transcript with two exons
# the clustering is performed according to the alignment of the
# internal boundaries of the exons

sub cluster_Introns{
    my ($self,$introns,$exact) = @_;
     
    # no point if there are no introns
    return unless ( scalar( @$introns) > 0 );
     
    # the total list of clusters:
    my $cluster_list = [];
 
    # keep track about in which cluster is each exon
    my %intron2cluster;
     
    # sort introns by start coordinate
    my @introns = sort { $self->transcript_low($a) <=> $self->transcript_low($b) } @$introns;
     
    # Create the first cluster
    my $intron_cluster = new ClusterMerge::TranscriptCluster;
     
    # Start off the cluster with the first exon
    $intron_cluster->put_Transcripts($introns[0]);
    
    $intron2cluster{$introns[0]} = $intron_cluster;                                                           
    push(@$cluster_list, $intron_cluster);
    
    # Loop over the rest of the introns
    my $count = 0;
    
  EXON:
    foreach my $intron (@introns) {
        if ($count > 0) {
	    
            # Add to cluster if overlap AND if strand matches
            if ( $self->compare_intron_to_intron_cluster( $intron, $intron_cluster , $exact) ){
                $intron_cluster->put_Transcripts($intron);
                $intron2cluster{$intron} = $intron_cluster;
            }
            else {
                # Start a new cluster
                $intron_cluster = new ClusterMerge::TranscriptCluster;
                $intron_cluster->put_Transcripts($intron);
                $intron2cluster{$intron} = $intron_cluster;
		
                # and add it to the main_cluster feature
                push(@$cluster_list,$intron_cluster);
            }
        }
        $count++;
    }
    return ($cluster_list,\%intron2cluster);
}

############################################################
# provides the unique introns from a
# a set of (possibly overlapping) transcript predictions

sub get_uniq_introns{
    my ($self,@trans) = @_;
    
    my @introns;
    foreach my $t (@trans){
	push(@introns, @{$t->get_all_Introns} );
    }
    my ($clusters,$i2c) = $self->cluster_Introns(\@introns,'exact');
    # each cluster represents one final intron
    my @final_introns;
    foreach my $cluster (@$clusters){
	my $intron = $self->get_consensus_intron($cluster);
	push(@final_introns, $intron);
    }    
    return @final_introns;
}

############################################################

sub get_uniq_introns_with_coding_flag{
    my ($self,@trans) = @_;
    
    my %label;
    my @introns;
    foreach my $tran (@trans){
	my @these_introns = @{$tran->get_all_Introns};
	
	if ($tran->coding_transcript){
	    foreach my $intron ( @these_introns ){
		foreach my $exon (@{$tran->get_all_Exons}){
		    foreach my $codex ( @{$tran->coding_transcript->get_all_Exons} ){
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
		}
	    }
	}
	else{
	    foreach my $intron ( @these_introns ){
		foreach my $exon (@{$tran->get_all_Exons}){
		    $label{$exon} = "non-coding";
		}
	    }
	}
	push(@introns, @these_introns);
    }
    
    my ($clusters,$i2c) = $self->cluster_Introns(\@introns,'exact');
    
    # each cluster represents one final intron
    my @final_introns;
    foreach my $cluster (@$clusters){
	my ($intron) = $self->get_consensus_intron_keep_flag($cluster,\%label);
	# check
        #foreach my $e (@{$intron->get_all_Exons}){
	    #print "flag: ".$label{$e}."\n";
	#}
	push(@final_introns, $intron);
    }    
    return (\@final_introns,\%label);
}




############################################################
# this method compares one intron to an intron cluster.
# It calls "compare_2introns" for the intron and for each
# intron in the cluster.

sub compare_intron_to_intron_cluster{
    my ($self,$intron, $intron_cluster, $exact) = @_;
    
    # $intron is a ClusterMerge::Transcript
    # $intron_cluster is a ClusterMerge::TranscriptCluster
    
    foreach my $i2 (@{$intron_cluster->get_Transcripts} ){
	
        if ( $self->compare_2introns($i2,$intron,$exact) ){
            return 1;
        }
    }
    return 0;
}

############################################################
# It performs the comparison between 2 introns.
# I checks whether two introns overlap:
#  exon|----|exon
# exon|-------------|exon
#
# I 'exact' is chosen, it checks that the intron
# boundaries are conserved:
#  exon|----|exon
#  exon|----|exon

# the intron boundaries (Note that the external exon
# boundaries are not considered here). If it is called
# with a flag ($exact), there is no mismatch allowed.

sub compare_2introns{
    my ($self,$i1, $i2, $exact) = @_;
    my @e1 = sort { $a->start <=> $b->start } @{$i1->get_all_Exons};
    my @e2 = sort { $a->start <=> $b->start } @{$i2->get_all_Exons};
    
    my $intron1_start = $e1[0]->end + 1;
    my $intron1_end   = $e1[1]->start - 1;
    
    my $intron2_start = $e2[0]->end + 1;
    my $intron2_end   = $e2[1]->start - 1;

    if ($exact){
	if ( $intron1_start == $intron2_start
	     &&
	     $intron1_end   == $intron2_end ){
	    return 1;
	}
	else{
	    return 0;
	}
    }
    else{
	if ( !( $intron1_start > $intron2_end
		||
		$intron1_end   < $intron2_start ) ){
	    return 1;
	}
	else{
	    return 0;
	}
    }
}

############################################################
# given a cluster of introns, represented as
# a ClusterMerge::TranscriptCluster,
# this method generates a consensus intron

sub get_consensus_intron{
    my ($self,$cluster) = @_;
    my @transcripts = @{$cluster->get_Transcripts};
    my %s1;
    my %s2;
    my %e1;
    my %e2;
    my %tid;
    my %type;
    my $target;
    my $strand;
    foreach my $t (@transcripts){
        $tid{$t->dbID}++;
        $type{$t->type}++;
	
        # we look at both exons from lower to hight coordinate
        my @e = sort {$a->start <=> $b->start} @{$t->get_all_Exons};
        $strand = $e[0]->strand unless ($strand);
        $target = $e[0]->seqname unless ($target);
        $s1{$e[0]->start}++;
        $s2{$e[1]->start}++;
        $e1{$e[0]->end}++;
        $e2{$e[1]->end}++;
    }
    
    my @starts1 = sort { $s1{$b} <=> $s1{$a} } keys %s1;
    my @starts2 = sort { $s2{$b} <=> $s2{$a} } keys %s2;
    my @ends1   = sort { $e1{$b} <=> $e1{$a} } keys %e1;
    my @ends2   = sort { $e2{$b} <=> $e2{$a} } keys %e2;
    
    # the the most common one, but if this only occurs
    # once, take the one that makes the exon shortest
    my ($s1,$s2);
    my ($e1,$e2);
    
    #  s1    e1           s2    e2
    #   ######-------------######
    #
    
    if ( $s1{$starts1[0]} > 1 ){
        $s1 = $starts1[0];
    }
    else{
        $s1 = max( keys %s1 );
    }
    
    if ( $s2{$starts2[0]} > 1 ){
        $s2 = $starts2[0];
    }
    else{
        $s2 = max( keys %s2 );
    }
    
    if ( $e1{$ends1[0]} > 1 ){
        $e1 = $ends1[0];
    }
   if ( $e1{$ends1[0]} > 1 ){
       $e1 = $ends1[0];
   }
    else{
        $e1 = min( keys %e1 );
    }
    
    if ( $e2{$ends2[0]} > 1 ){
        $e2 = $ends2[0];
    }
    else{
        $e2 = min( keys %e2 );
    }
    
    # but if the resulting exon is too short
    # take an exon of 60bp
    if ( ($e1 - $s1 + 1) < 60 ){
        $s1 = $e1 - 60 + 1;
    }
    if ( ($e2 - $s2 + 1) < 60 ){
        $e2 = $s2 + 60 - 1;
    }
    
    my $consensus = ClusterMerge::Transcript->new();
    my $tid = join ":", sort keys %tid;
    my $type = join ":" , sort keys %type;
    $consensus->dbID($tid);
    $consensus->type($type);
    
    my $exon1 = ClusterMerge::Exon->new();
    $exon1->seqname($target);
    $exon1->source_tag($type);
    $exon1->primary_tag("exon");
    $exon1->start($s1);
    $exon1->end($e1);
    $exon1->strand($strand);
    $exon1->transcript_tag($tid);
    
    my $exon2=  ClusterMerge::Exon->new();
    $exon2->seqname($target);
    $exon2->source_tag($type);
    $exon2->primary_tag("exon");
    $exon2->start($s2);
    $exon2->end($e2);
    $exon2->strand($strand);
    $exon2->transcript_tag($tid);
    
    $consensus->add_Exon($exon1);
    $consensus->add_Exon($exon2);
    
    return $consensus;
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
# given a cluster of introns, represented as
# a ClusterMerge::TranscriptCluster,
# this method generates a consensus intron

sub get_consensus_intron_keep_flag{
    my ($self,$cluster,$label) = @_;
    my @transcripts = @{$cluster->get_Transcripts};
    my %s1;
    my %s2;
    my %e1;
    my %e2;
    my %tid;
    my %type;
    my $target;
    my $strand;
    
    my $left_label;
    my $right_label;
    foreach my $t (@transcripts){
        $tid{$t->dbID}++;
        $type{$t->type}++;
	
        # we look at both exons from lower to hight coordinate
        my @e = sort {$a->start <=> $b->start} @{$t->get_all_Exons};
        $strand = $e[0]->strand unless ($strand);
        $target = $e[0]->seqname unless ($target);
        $s1{$e[0]->start}++;
        $s2{$e[1]->start}++;
        $e1{$e[0]->end}++;
        $e2{$e[1]->end}++;
    
	############################################################
	# choose a consensus label for the exons
	# In case of clash the priority is coding > half-coding > utr
	unless ($left_label){
	    $left_label = $label->{$e[0]};
	}
	if ($label->{$e[0]} eq "coding"){
	    $left_label = "coding";
	}
	elsif( $label->{$e[0]} eq "half-coding" 
	       &&
	       ($left_label eq "utr" 
		||
		$left_label eq "non-coding")
	       ){
	    $left_label = "half-coding";
	}
	unless ($right_label){
	    $right_label = $label->{$e[1]};
	}
	if ($label->{$e[1]} eq "coding"){
	    $right_label = "coding";
	}
	elsif( $label->{$e[1]} eq "half-coding" 
	       &&
	       ($right_label eq "utr" 
		||
		$right_label eq "non-coding")
	       ){
	    $right_label = "half-coding";
	}
	
    }
    
    my @starts1 = sort { $s1{$b} <=> $s1{$a} } keys %s1;
    my @starts2 = sort { $s2{$b} <=> $s2{$a} } keys %s2;
    my @ends1   = sort { $e1{$b} <=> $e1{$a} } keys %e1;
    my @ends2   = sort { $e2{$b} <=> $e2{$a} } keys %e2;
    
    # the the most common one, but if this only occurs
    # once, take the one that makes the exon shortest
    my ($s1,$s2);
    my ($e1,$e2);
    
    #  s1    e1           s2    e2
    #   ######-------------######
    #
    
    if ( $s1{$starts1[0]} > 1 ){
        $s1 = $starts1[0];
    }
    else{
        $s1 = max( keys %s1 );
    }
    
    if ( $s2{$starts2[0]} > 1 ){
        $s2 = $starts2[0];
    }
    else{
        $s2 = max( keys %s2 );
    }
    
    if ( $e1{$ends1[0]} > 1 ){
        $e1 = $ends1[0];
    }
   if ( $e1{$ends1[0]} > 1 ){
       $e1 = $ends1[0];
   }
    else{
        $e1 = min( keys %e1 );
    }
    
    if ( $e2{$ends2[0]} > 1 ){
        $e2 = $ends2[0];
    }
    else{
        $e2 = min( keys %e2 );
    }
    
    # but if the resulting exon is too short
    # take an exon of 60bp
    if ( ($e1 - $s1 + 1) < 60 ){
        $s1 = $e1 - 60 + 1;
    }
    if ( ($e2 - $s2 + 1) < 60 ){
        $e2 = $s2 + 60 - 1;
    }
    
    my $consensus = ClusterMerge::Transcript->new();
    my $tid = join ":", sort keys %tid;
    my $type = join ":" , sort keys %type;
    $consensus->dbID($tid);
    $consensus->type($type);
    
    # the 'left' exon
    my $exon1 = ClusterMerge::Exon->new();
    $exon1->seqname($target);
    $exon1->source_tag($type);
    $exon1->primary_tag("exon");
    $exon1->start($s1);
    $exon1->end($e1);
    $exon1->strand($strand);
    $exon1->transcript_tag($tid);
    $label->{$exon1} = $left_label;

    # the 'right' exon
    my $exon2=  ClusterMerge::Exon->new();
    $exon2->seqname($target);
    $exon2->source_tag($type);
    $exon2->primary_tag("exon");
    $exon2->start($s2);
    $exon2->end($e2);
    $exon2->strand($strand);
    $exon2->transcript_tag($tid);
    $label->{$exon2} = $right_label;
    
    $consensus->add_Exon($exon1);
    $consensus->add_Exon($exon2);
    
    return ($consensus,$label);
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



1;
############################################################
