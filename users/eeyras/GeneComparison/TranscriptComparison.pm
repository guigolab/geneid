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

package GeneComparison::TranscriptComparison;

use ClusterMerge::Transcript;
use ClusterMerge::TranscriptCluster;
use ClusterMerge::ExonUtils;
use ClusterMerge::ObjectMap;
use ClusterMerge::Root;
use GeneComparison::Alineator;
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
# Method to pair-up transcripts from two lists
# THis is based on finding the best possible
# pair by the stable-marriage algorithm
# according to the correlation coefficient for 
# the nucleotide comparison (CCn)
#
# The method returns a list of lists. Each sublist
# contains the pair of transcripts, the CCn and the
# SNn and SPn.

sub pair_Transcripts{
    my ($self,$predictions,$annotations) = @_;

    my $verbose = 1;

    # cache the transcript lengths
    my %transcript_length;

    # cache the sensitivity and specificity
    my %sn;
    my %sp;
    my %TP;
    my %FP;
    my %FN;
    
    # use an ObjectMap to hold the possible pairs
    my $object_map = ClusterMerge::ObjectMap->new();
    
    # first calculate all possible overlaps
    foreach my $pred_tran ( @$predictions ){
	unless ( $transcript_length{$pred_tran} ){
	    $transcript_length{$pred_tran} = $self->transcript_length( $pred_tran );
	}
	
	foreach my $ann_tran ( @$annotations ){
	    unless ( $transcript_length{$ann_tran} ){
		$transcript_length{$ann_tran} = $self->transcript_length( $ann_tran );
	    }
	    
	    my ($overlap_number,$overlap_length) = 
		ClusterMerge::TranscriptComparator->_compare_Transcripts( $pred_tran, $ann_tran );
	    
	    print "ann: $transcript_length{$ann_tran}, pred: $transcript_length{$pred_tran}\n" if $verbose;
	    #ClusterMerge::TranscriptUtils->_print_SimpleTranscript($ann_tran);
            #ClusterMerge::TranscriptUtils->_print_SimpleTranscript($pred_tran);
	    
	    my $TP = $overlap_length;
	    my $FP = $transcript_length{$pred_tran} - $TP;
	    my $FN = $transcript_length{$ann_tran}  - $TP;
	    my $TN = $self->calculate_TrueNegatives( $pred_tran, $ann_tran);

	    $TP{$pred_tran}{$ann_tran} = $TP;
	    $FP{$pred_tran}{$ann_tran} = $FP;
	    $FN{$pred_tran}{$ann_tran} = $FN;
   
	    my $SNn = $overlap_length / $transcript_length{$ann_tran}; # $SNn = $TP/ ($TP + $FN);
	    my $SPn = $overlap_length / $transcript_length{$pred_tran};# $SPn = $TP/ ($TP + $FP);

	    # my $CCn = ( $SNn + $SPn )/2;  # look up the exact formula

	    my $CCn = (($TP*$TN) - ($FN*$FP))/sqrt(($TP+$FN)*($TN+$FP)*($TP+$FP)*($TN+$FN));
	    print "SN = $SNn\tSP = $SPn\tTP = $TP\tTN = $TN\tFP = $FP\tFN = $FN\tCCn = $CCn\n" if $verbose;
	    
    
	    #$sn{$pred_tran}{$ann_tran} = $SNn;
	    #$sp{$pred_tran}{$ann_tran} = $SPn;
	    $sn{$pred_tran}{$ann_tran} = $SNn;
	    $sp{$pred_tran}{$ann_tran} = $SPn;

    
	    # put the scores in an ObjbectMap 
	    # and apply the stable-marriage algorithm to find the
	    # best possible pairs

	    # what score to use? CCn?
	    $object_map->match($pred_tran, $ann_tran, $CCn);
	}
    }
    
    my $pairs_object_map = $object_map->stable_marriage;
    
    # get the pred_tran from each formed pair:
    my @list1 = $pairs_object_map->list1();
    
    my @pairs_list;
    foreach my $pred ( @list1 ){
	my @anns = $pairs_object_map->partners( $pred );
	push (@pairs_list, 
	      [$pred,$anns[0],$pairs_object_map->score($pred,$anns[0]),$TP{$pred}{$anns[0]},$FP{$pred}{$anns[0]}, $FN{$pred}{$anns[0]} ] );
    }
    return @pairs_list;
}

############################################################
# this method currently only calculates
# the extension of the sequence which is
# annotated as non-exonic by both prediction and annotation.
# This is in fact, whatever is no exonic between both exon sets.
# We calculate the length between the first and the last introns.
# Use padding?
sub calculate_TrueNegatives{
    my ($self,$pred_tran,$ann_tran) = @_;
    my @pred_exons = @{$pred_tran->get_all_Exons};
    my @ann_exons  = @{$ann_tran->get_all_Exons};
    my @exons = (@pred_exons,@ann_exons);
    my ($clusters,$exon2cluster) = ClusterMerge::ExonUtils->_cluster_Exons(@exons);
    
    my $padding = 0; # use padding?
    my $length = $padding;
    my @clusters = sort {$a->start <=> $b->start} @$clusters;
    for(my $i=0; $i< scalar(@clusters)-1; $i++ ){
	$length += $clusters[$i+1]->start - $clusters[$i]->end - 1;
    }
    #print "returning length $length\n";
    $length = 1000 unless ($length);
    return $length;
}

############################################################

sub transcript_length{
    my ($self,$t) = @_;
    my $length = 0;
    foreach my $e ( @{$t->get_all_Exons} ){
	$length += ($e->end - $e->start + 1 );
    }
    #print "returning length $length\n";
    return $length;
}

############################################################
# the method returns 1 if two transcripts are identical
# and 0 otherwise
sub exact_match{
    my ($self,$t1,$t2) = @_;
    my @e1 = sort{$a->start <=> $b->start} @{$t1->get_all_Exons};
    my @e2 = sort{$a->start <=> $b->start} @{$t2->get_all_Exons};
    return 0 unless (scalar(@e1) == scalar(@e2) );
    for(my $i=0; $i<scalar(@e1); $i++){
	return 0 unless ($e1[$i]->start == $e2[$i]->start
			 &&
			 $e1[$i]->end   == $e2[$i]->end
			 );
    }
    return 1;
}



############################################################

sub compare_fingerprints{
    my ($self,$t1,$t2) = @_;
    
    my $f1 = $self->get_fingerprint($t1);
    my $f2 = $self->get_fingerprint($t2);
    
    GeneComparison::Alineator->compare_fingerprints($t1->dbID,
						    scalar(@{$t1->get_all_Exons}),
						    $f1,
						    $t2->dbID,
						    scalar(@{$t2->get_all_Exons}),
						    $f2
						    );
}

############################################################

sub get_fingerprint{
    my ($self,$tran) = @_;

    my $verbose = 0;

    my $t = $tran->coding_transcript;
    
    my $id = $t->dbID;       
    my @e = @{$t->get_all_Exons};
    my $nex   = scalar(@e);
    
    my @fingerprint;

    push( @fingerprint, $id);
    push( @fingerprint, $nex);

    my $fase5=3;      #fase 5'
    my $fase3;        #fase 3'
    my $d=0;          #Es 3-fase3,lo que (nucleotidos que le faltan al ultimo codon del anterior exon)
                      # le restamos al tamanyo del siguiente exon antes de especificar su fase.

    my @exons;
    my $strand = $e[0]->strand;
    if ($strand == 1 ){
	@exons = sort {$a->start <=> $b->start} @e;
    }
    elsif ($strand == -1){
	@exons = sort {$b->start <=> $a->start} @e;
    }
    foreach my $e ( @exons){
	my $length = ($e->end - $e->start + 1) - $d;

	# phase
	if ($length%3==0){
	    $fase3= 0;
	}
	elsif ($length%3==1){
	    
	    $fase3= 1;
	}
	elsif ($length%3==2){
	    $fase3= 2;
	}
	
	if ( $e == $exons[-1] ){
	    $fase3 = 3;
	}

	#print el fingerprint de cada exon;

	my $exonfinger = $fase5.":".$fase3.":".($length+$d);
	push (@fingerprint, $exonfinger);
	
	$fase5= $fase3;       # la fase 5' de un exon es la misma que la 3' del anterior.
	$d= 3-$fase3;         # preparamos el modificador de tamaño del exon $d.
	
    }
    
    print "FINGERPRINT:\t@fingerprint\n" if $verbose;
    return \@fingerprint;

}




1;
############################################################
