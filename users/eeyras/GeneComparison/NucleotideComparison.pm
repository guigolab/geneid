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

package GeneComparison::NucleotideComparison;

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
# methods related to NUCLEOTIDE COMPARISON
############################################################

sub compare_nucleotides_Gene_level{
    my ($self,$predictions,$annotations) = @_;
    
    #print "comparing nucleotides between ".scalar(@$predictions)." predictions and ".scalar(@$annotations)." annotations\n";
    # each array is a list of transcript clusters (ClusterMerge::TranscriptCluster objts)
    my @predictions = @$predictions;
    my @annotations = @$annotations;
    my ($TP,$FP,$FN) = (0,0,0);
    
    if ( !@annotations && @predictions){
	my @projections = $self->get_Projections(\@predictions); # projections are ClusterMerge::ExonCluster objects
	$TP = 0;
	$FN = 0;
	foreach my $projection (@projections){
	    $FP += ($projection->end - $projection->start + 1 );
	}	    
    }
    elsif ( @annotations && !@predictions){
	my @projections = $self->get_Projections(\@annotations);
	$TP = 0;
	$FP = 0;
	foreach my $projection (@projections){
	    $FN += ($projection->end - $projection->start + 1 );
	}
    }
    elsif(  @annotations && @predictions ){
	my @pred_projections = $self->get_Projections(\@predictions);
	my @ann_projections  = $self->get_Projections(\@annotations);
	($TP,$FP,$FN) = $self->compare_Projections(\@pred_projections, \@ann_projections);
    }
    return ($TP,$FP,$FN);
}

############################################################

sub get_Projections{
    my ($self,$arrayref) = @_; # arrayref is an arrayref of ClusterMerge::TranscriptCluster objects
    my @trans;
    foreach my $c ( @$arrayref ){
	push( @trans, @{$c->get_Transcripts} );
    }
    return $self->get_Transcript_Projections(\@trans);
}

############################################################

sub get_Transcript_Projections{
    my ($self,$trans) = @_;   # $trans is an arrayref of ClusterMerge::Transcript objects
    my @exons;
    foreach my $t ( @$trans ){
       push( @exons, @{$t->get_all_Exons} );
   } 
    #print "eval: clustering ".scalar(@exons)." exons\n";
    my ($clusters,$exon2cluster) = ClusterMerge::ExonUtils->_cluster_Exons(@exons);
    return @$clusters;
}

############################################################
# method to compare the projections of the annotations and predictions
# These projections are ClusterMerge::ExonClusters and they function
# as ranges

sub compare_Projections{
    my ($self,$pred_proj,$ann_proj) = @_;
    
    my $tot_pred = 0;
    foreach my $pred (@$pred_proj){
	$tot_pred += ($pred->end - $pred->start + 1 );
    }
    
    my $tot_ann = 0;
    my $found = 0;
    # each $ann and $pred are ExonCluster, which function as Ranges (things with start and end)
    foreach my $ann ( @$ann_proj ){
	$tot_ann += ( $ann->end - $ann->start + 1 );
	foreach my $pred ( @$pred_proj ){
	    if ( $ann->overlaps($pred) ){
		my ($st,$en) = $ann->intersection($pred);
		$found += $en - $st + 1;
	    }
	}
    }
    my $TP = $found;
    my $FP = $tot_pred - $found;
    my $FN = $tot_ann  - $found;
    return ($TP,$FP,$FN);
}
    
1;
############################################################
