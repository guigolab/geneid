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

package GeneComparison::TranscriptPair;

use vars qw(@ISA);
use strict;
use ClusterMerge::Transcript;
use ClusterMerge::TranscriptCluster;
use ClusterMerge::Root;
use ClusterMerge::GFFTools;
use Bio::Seq;
use Bio::SeqIO;
use ClusterMerge::TranscriptUtils;
use Adaptor::SequenceAdaptor;

@ISA = qw(ClusterMerge::Root);

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);
  
    my( $gap_penalty ) = $self->_rearrange([qw(
					       GAP_PENALTY
					       )], 
					   @args);
 
    if ( $gap_penalty ){
	$self->gap_penalty( $gap_penalty );
    }
    else{
	$self->gap_penalty( -100 );
    }
    
    return $self;
}

############################################################
#
# compare transcripts with blastn
#
############################################################


sub blast_isoforms{
    my ( $self,$tran1,$tran2, $dir1, $dir2, $coding_exons ) = @_;
    
    my $verbose = 1;

    my @exons1 = @{$tran1->get_all_Exons};
    my $seqname1 = $exons1[0]->seqname;
    $seqname1 = "chr".$seqname1 unless $seqname1=~/chr/;
    my $chrfile1 = $dir1."/".$seqname1.".fa";
	

    my @exons2 = @{$tran2->get_all_Exons};
    my $seqname2 = $exons2[0]->seqname;
    $seqname2 = "chr".$seqname2 unless $seqname2=~/chr/;
    my $chrfile2 = $dir2."/".$seqname2.".fa";
    # query
    my $id1 = $tran1->dbID;
    #target
    my $id2 = $tran2->dbID;
    
    print "comparing $id1 and $id2\n";
    
    my ($seq1, $seq2 );
    if ( $coding_exons ){
	my $string1 = Adaptor::SequenceAdaptor->get_transcript_seq($chrfile1, $tran1->coding_transcript);
	$seq1       = Bio::Seq->new(
				    -DISPLAY_ID => $id1,
				    -MOLTYPE    => 'dna',
				    -SEQ        => $string1,
				    );
	
	my $string2 = Adaptor::SequenceAdaptor->get_transcript_seq($chrfile2,$tran2->coding_transcript);
	$seq2       = Bio::Seq->new(
				    -DISPLAY_ID => $id2,
				    -MOLTYPE    => 'dna',
				    -SEQ        => $string2,
				    );
    }
    else{
	my $string1 = Adaptor::SequenceAdaptor->get_transcript_seq($chrfile1,$tran1);
	$seq1       = Bio::Seq->new(
				    -DISPLAY_ID => $id1,
				    -MOLTYPE    => 'dna',
				    -SEQ        => $string1,
				    );
	
	my $string2 = Adaptor::SequenceAdaptor->get_transcript_seq($chrfile2,$tran2);
	$seq2       = Bio::Seq->new(
				    -DISPLAY_ID => $id2,
				    -MOLTYPE    => 'dna',
				    -SEQ        => $string2,
				    );
    }
    my $length1 = $seq1->length;
    my $length2 = $seq2->length;
    
    ############################################################
    # create database
    my $file = 'seq_'.$$.'.fa';
    my $database = "/tmp/".$file;
    open( DB_SEQ,">$database") || die("Could not open $database $!");
    
    my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
				 '-fh'     => \*DB_SEQ);
    
    $seqout->write_seq($seq2);
    close( DB_SEQ );
    system("pressdb $database > /dev/null 2>&1");
    
    ############################################################
    # create query file
    my $query_file = '/tmp/query_'.$$.'.fa';
    open ( QUERY, ">$query_file") or die ("cannot open query file $query_file for writing");
    
    my $seqout2 = Bio::SeqIO->new('-format' => 'Fasta',
				  '-fh'     => \*QUERY);
    $seqout2->write_seq($seq1);
    close(QUERY);

    ############################################################
    # generate blast command
    
    #my $options = "-nogap W=5";
    # Ian's parameters:
    #my $options = "W=5 M=1 N=-1 Q=3 R=3";
    my $options = "W=5 -topcomboN 1";
    
    my $command = "blastn $database $query_file $options";

    open( OUT, "$command | parseblast.pl -GF |" );
    #open( OUT, "$command | ");
    
    my @featurepairs;
    while(<OUT>){
	chomp;
	my $e = ClusterMerge::GFFTools->exon_from_gff($_);
	push(@featurepairs, $e);
    }

    unlink ( $database );
    unlink ($query_file);
    
    ############################################################
    # separate by strands

    my @features;
    push( @{$features[0]}, grep { $_->strand == 1 } @featurepairs );
    push( @{$features[1]}, grep { $_->strand == -1} @featurepairs );

    my $best_score = 0;
    my $best_features;
    for (my $i=0; $i<2; $i++ ){
	unless ( $features[$i] && @{$features[$i]} ){
	    next;
	}
	
	############################################################
	# compute the score
	my $score = 0;
	foreach my $fp (sort {$a->start <=> $b->start} @{$features[$i]}) {
	    $score += $fp->score;
	    #print $fp->gffstring . "\n";
	}
	if ( $score > $best_score ){
	    $best_score    = $score;
	    $best_features = $features[$i];
	}
    }
    
    my $coverage = 0;
    my $spliced = 0;
    
    if ( $best_features ){
      
	############################################################
	# calculate coverage
	# we use query/target as in feature pairs the target=seqname and query=hseqname
	my ($query_coverage,  $query_spliced)  = 
	    $self->process_query( $best_features, $tran2 , $seq1, $coding_exons);
	my ($target_coverage, $target_spliced) = 
	    $self->process_target( $best_features, $tran1, $seq2, $coding_exons);
	$coverage = ( $query_coverage + $target_coverage )/2;
	
	if ( $query_spliced || $target_spliced ){
	    $spliced = 1;
	}
	
	############################################################
	# calculate the perc id
	#my $perc_id = 0;
	#foreach my $f ( @$best_features ){
	#    $perc_id += $f->percent_id;
	#}
	#$perc_id = sprintf "%.2f", ( $perc_id/scalar(@$best_features) );
	
	print STDERR "\tquery:$id1 coverage:$query_coverage spliced:$query_spliced\n";
	print STDERR "\ttarget:$id2 coverage:$target_coverage spliced:$target_spliced\n";
	#print STDERR "\taveraged percent id: $perc_id\n";
	$best_score = 0 if ( $target_coverage + $query_coverage < 100 );
    }
    
    $best_score = 0 if $spliced;
    return ( $best_score, $best_features );
}


############################################################
# the target
############################################################

sub process_target{
    my ($self,$feat, $tran, $seq, $coding_exons, $protein) = @_;
    
    my $transcript_length = length($seq);
    my @exons;
    
    if ( $coding_exons || $protein ){
	@exons = sort { $a->length <=> $b->length } @{$tran->coding_transcript->get_all_Exons};
    }
    else{
	@exons = sort { $a->length <=> $b->length } @{$tran->get_all_Exons};
    }
    my $min_exon_length = $exons[0]->length;
    if ( $protein && $coding_exons){
	if ( $seq=~/TAA$|TGA$|TAG$/i ){
	    $transcript_length = ($transcript_length - 3)/3;
	}
	else{
	    $transcript_length /= 3;
	}
	$min_exon_length /= 3;
    }  
    
    my $is_spliced;
    
    my @clusters;
    my @cluster_starts;
    my @cluster_ends;
    my @features = sort{ $a->start <=> $b->start} @$feat;
 
    # create the first cluster
    my $count = 0;
    my $cluster = [];
    
    # start it off with the first feature
    my $first_feat = shift( @features );
    push (@$cluster, $first_feat);
    $cluster_starts[$count] = $first_feat->start;
    $cluster_ends[  $count] = $first_feat->end;
 
    # store the list of clusters
    push(@clusters,$cluster);
    
    ############################################################
    # loop over the rest of the features
  FEATURE:
    foreach my $f ( @features ){
	if (!($f->end < $cluster_starts[$count] || $f->start > $cluster_ends[$count])) {      
	    push(@$cluster,$f);
	    
	    # re-adjust size of cluster
	    if ($f->start < $cluster_starts[$count]) {
	      $cluster_starts[$count] = $f->start;
	    }
	    if ($f->end  > $cluster_ends[$count]) {
	      $cluster_ends[$count]   = $f->end;
	    }
	  }
	else{
	    # else, start create a new cluster with this feature
	    $count++;
	    $cluster = [];
	    push (@$cluster, $f);
	    $cluster_starts[$count] = $f->start;
	    $cluster_ends[  $count] = $f->end;
	    
	    # store it in the list of clusters
	    push(@clusters,$cluster);
	}
    }
    
    ############################################################
    # check whether the transcript has one or more exons unaligned
    if ( scalar( @clusters ) == 1 ){
	$is_spliced = 0;
    }
    else{
	# compute the size of the 'gaps'
	my @gaps;
	$is_spliced = 0;
	for(my $i=0; $i<$#clusters; $i++){
	    my $gap = $cluster_starts[$i+1] - $cluster_ends[$i] - 1;
	    print STDERR "gap: $gap, min_exon_length = $min_exon_length\n";
	    if ( $gap >= $min_exon_length ){
		$is_spliced = 1;
		#print STDERR "is spliced\n";
	    }
	}
    }
    
    ############################################################
    # calculate the coverage of the transcript
    my $feature_length = 0;
    for(my $i=0; $i<=$#clusters; $i++){
	#print STDERR "target cluster $i: $cluster_starts[$i] - $cluster_ends[$i]\n";
	$feature_length += $cluster_ends[$i] - $cluster_starts[$i] + 1;
    }
    #print STDERR "feature length   : $feature_length\n";
    #print STDERR "transcript length: $transcript_length\n";
    my $coverage = sprintf "%.2f", 100*$feature_length/$transcript_length;
    #$self->print_exons_in_transcript($tran);
    return ($coverage,$is_spliced);
   
}

############################################################
# the query 
############################################################
sub process_query{
    my ($self,$feat, $tran, $seq, $coding_exons, $protein) = @_;
    
    my $transcript_length = length($seq);
    my @exons;
    
    if ( $coding_exons || $protein ){
	@exons = sort { $a->length <=> $b->length } @{$tran->coding_transcript->get_all_Exons};
    }
    else{
	@exons = sort { $a->length <=> $b->length } @{$tran->get_all_Exons};
    }
    my $min_exon_length = $exons[0]->length;
    if ( $protein && $coding_exons){
	if ( $seq=~/TAA$|TGA$|TAG$/i ){
	    $transcript_length = ($transcript_length - 3)/3;
	}
	else{
	    $transcript_length /= 3;
	}
	$min_exon_length /= 3;
    }  
    
    my $is_spliced;
    
    my @clusters;
    my @cluster_hstarts;
    my @cluster_hends;
    #my @features = sort{ $a->hstart <=> $b->hstart} @$feat;
    my @features = sort{ $a->start <=> $b->start} @$feat;
    
    # create the first cluster
    my $count = 0;
    my $cluster = [];
    
    # start it off with the first feature
    my $first_feat = shift( @features );
    push (@$cluster, $first_feat);
    #$cluster_hstarts[$count] = $first_feat->hstart;
    #$cluster_hends[  $count] = $first_feat->hend;
    
    $cluster_hstarts[$count] = $first_feat->start;
    $cluster_hends[  $count] = $first_feat->end;
    
  # store the list of clusters
    push(@clusters,$cluster);
    
    ############################################################
    # loop over the rest of the features
  FEATURE:
    foreach my $f ( @features ){
	if (!($f->end < $cluster_hstarts[$count] || $f->start > $cluster_hends[$count])) {      
	    push(@$cluster,$f);
	    
	    #if (!($f->hend < $cluster_hstarts[$count] || $f->hstart > $cluster_hends[$count])) {      
	    
	    # re-adjust size of cluster
	    #if ($f->hstart < $cluster_hstarts[$count]) {
	#	$cluster_hstarts[$count] = $f->hstart;
	#    }
	#    if ($f->hend  > $cluster_hends[$count]) {
	#	$cluster_hends[$count] = $f->hend;
	#    }
	    # re-adjust size of cluster
	    if ($f->start < $cluster_hstarts[$count]) {
		$cluster_hstarts[$count] = $f->start;
	    }
	    if ($f->end  > $cluster_hends[$count]) {
		$cluster_hends[$count] = $f->end;
	    }
	}
	else{
	    # else, start create a new cluster with this feature
	    $count++;
	    $cluster = [];
	    push (@$cluster, $f);
	    $cluster_hstarts[$count] = $f->start;
	    $cluster_hends[$count]   = $f->end;
	    #$cluster_hstarts[$count] = $f->hstart;
	    #$cluster_hends[$count]   = $f->hend;
	 
	 # store it in the list of clusters
	    push(@clusters,$cluster);
	}
    }
    
    ############################################################
    # check whether the transcript has one or more exons unaligned
    if ( scalar( @clusters ) == 1 ){
     $is_spliced = 0;
 }
 else{
   # compute the size of the 'gaps'
     my @gaps;
     $is_spliced = 0;
     for(my $i=0; $i<$#clusters; $i++){
	 my $gap = $cluster_hstarts[$i+1] - $cluster_hends[$i] - 1;
	 print STDERR "gap: $gap, min_exon_length = $min_exon_length\n";
	 if ( $gap >= $min_exon_length ){
	     $is_spliced = 1;
	     #print STDERR "is spliced\n";
	 }
   }
 }
  
 ############################################################
  # calculate the coverage of the transcript
  my $feature_length = 0;
  for(my $i=0; $i<=$#clusters; $i++){
    #print STDERR "query cluster $i: $cluster_hstarts[$i] - $cluster_hends[$i]\n";
    $feature_length += $cluster_hends[$i] - $cluster_hstarts[$i] + 1;
  }
    my $coverage = sprintf "%.2f", 100*$feature_length/$transcript_length;
    #print STDERR "coverage = $feature_length / $transcript_length = $coverage\n";
    
    #$self->print_exons_in_transcript($tran);
    
    return ($coverage,$is_spliced);
}

############################################################

1;
