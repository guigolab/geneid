#!/usr/local/bin/perl  -w

############################################################
#
# script to find novel exon assemblies in a prediction
# set compared to another prediction/annotation
#
# written by Eduardo Eyras (eae@imim.es)
#
############################################################

use strict;
use Getopt::Long;
use ClusterMerge::Transcript;
use ClusterMerge::TranscriptUtils;
use ClusterMerge::Exon;
use ClusterMerge::GFFTools;
use GeneComparison::GeneComparison;

my $gff_format = 1;
my $input;
my $output;

my $verbose = 0;
my %type;


############################################################
# deafult options

my $pred_file;
my $ann_file;
my $h;

############################################################

&GetOptions( 
	     'p:s' => \$pred_file,
	     'a:s' => \$ann_file,
	     'h:s' => \$h,                          
	     );

if ( !($pred_file && $ann_file) ||  $h ){
    print STDERR "Usage: $0 -p <predictions> -a <annotations> -h <help>\n";
    exit(0);
}

############################################################
# read predictions


open ( IN, "<$pred_file") || die("could not open input file $pred_file");

my @predictions;
my $pred_trans_tag;
my %pred_trans;

my %forward_sets;
my %reverse_sets;

my %pred_gene_sets;

my %tag_count;
while(<IN>){
    chomp;
    my @e = split;
    next unless ( $e[3] && $e[4] ); 
    my %trans_tag;
    my $exon = ClusterMerge::GFFTools->exon_from_gff ($_);
    $exon->source_tag("prediction");
    if ( $exon && $exon->start && $exon->end ){
	push( @{$pred_trans{$exon->group_tag}}, $exon );
    }
}

foreach my $tag ( keys %pred_trans ){
    my $transcript = ClusterMerge::Transcript->new();
    $transcript->dbID( $tag );
    my $seqname;
    my $strand;
    foreach my $exon ( @{ $pred_trans{$tag} } ){
	$transcript->add_Exon($exon);
	$seqname = $exon->seqname;
	$strand = $exon->strand;
    }
#    push( @predictions, $transcript );
    $type{$transcript} = "prediction";
    if ($strand == -1 ){
	push( @{$reverse_sets{$seqname}}, $transcript );
    }
    else{
	push( @{$forward_sets{$seqname}}, $transcript );
    }
}
close(IN);


############################################################
# read annotations

open ( IN, "<$ann_file") || die("could not open input file $ann_file");
my @annotations;
my $ann_trans_tag;
my %ann_trans;

my %ann_gene_sets;
while(<IN>){
    chomp;
    my @e = split;
    next unless ( $e[3] && $e[4] ); 
    my $exon = ClusterMerge::GFFTools->exon_from_gff ($_);
    $exon->source_tag("annotation");
    if ( $exon && $exon->start && $exon->end ){
	push( @{$ann_trans{$exon->group_tag}}, $exon );
    }
}


foreach my $tag ( keys %ann_trans ){
    my $transcript = ClusterMerge::Transcript->new();
    $transcript->dbID( $tag );
    my $strand;
    my $seqname;
    foreach my $exon ( @{ $ann_trans{$tag} } ){
	$transcript->add_Exon($exon);
	$seqname = $exon->seqname;
	$strand  = $exon->strand;
    }
    #push( @annotations, $transcript );
    $type{$transcript} = "annotation";
    if ($strand == -1 ){
	push( @{$reverse_sets{$seqname}}, $transcript );
    }
    else{
	push( @{$forward_sets{$seqname}}, $transcript );
    }
}
close(IN);


############################################################
# cluster transcripts

my %novel_predictions;
my %mixed_clusters;
my %missed_annotations;


############################################################
# clustering in the forward strand
#
# the clustering of transcripts is done by genomic-overlap
# i.e., a transcript is considered as a range in the genomic sequence
#
print "Forward clustering\n" if $verbose;

foreach my $seqname ( keys %forward_sets ){
    print "Looking at target: $seqname\n" if $verbose;
    if ( $forward_sets{$seqname} ){
	my @all_trans = @{$forward_sets{$seqname}};
	my @forward_clusters = GeneComparison::GeneComparison->cluster_Transcripts(\@all_trans);
	
	print scalar(  @forward_clusters )." clusters found\n" if $verbose;
	my $count = 1;
	foreach my $c ( @forward_clusters ){
	    my @predictions = grep { $type{$_} eq "prediction" } @{$c->get_Transcripts};
	    my @annotations = grep { $type{$_} eq "annotation" } @{$c->get_Transcripts};
	    
	    if ( !@annotations && @predictions){
		push( @{$novel_predictions{$seqname}}, @predictions );
		print scalar( @predictions )." novel predictions in target $seqname\n" if $verbose;
	    }
	    #elsif ( @annotations && !@predictions){
	    # push( @{$missed_annotations{$seqname}}, @annotations );
	    # print scalar( @annotations )." missed annotations in target $seqname\n" if $verbose;
	    #}
	    elsif(  @annotations && @predictions ){
		push( @{$mixed_clusters{$seqname}}, $c ); 
		if ($verbose){
		    print "Cluster $count:\n";
		    $c->print_Cluster(\%type);
		    $count++;
		}
	    }
	}
	
    }
}

############################################################
# clustering in the reverse strand
#
# the clustering of transcripts is done according to genomic overlap
#
print "*********Reverse clustering\n" if $verbose;
foreach my $seqname ( keys %reverse_sets ){
    print "Looking at target: $seqname\n" if $verbose;
    if ( $reverse_sets{$seqname} ){
	my @all_trans = @{$reverse_sets{$seqname}};
	my @reverse_clusters = GeneComparison::GeneComparison->cluster_Transcripts(\@all_trans);
	print scalar(  @reverse_clusters )." clusters found\n" if $verbose;
	my $count = 1;
	foreach my $c ( @reverse_clusters ){
	    my @predictions = grep { $type{$_} eq "prediction" } @{$c->get_Transcripts};
	    my @annotations = grep { $type{$_} eq "annotation" } @{$c->get_Transcripts};
	    
	    if ( !@annotations && @predictions){
		push( @{$novel_predictions{$seqname}}, @predictions );
	    	print scalar( @predictions )." novel predictions in target $seqname\n" if $verbose;
	    }
	    #elsif ( @annotations && !@predictions){
	    #	push( @{$missed_annotations{$seqname}}, @annotations );
	    #	print scalar( @annotations )." missed annotations in target $seqname\n" if $verbose;
	    #}
	    elsif(  @annotations && @predictions ){
		push( @{$mixed_clusters{$seqname}}, $c ); 
		if ($verbose){
		    print "Cluster $count:\n";
		    $c->print_Cluster(\%type);
		    $count++;
		}
	    }
	}
    }
}


# %novel_predictions contains all novel predictions (no exon overlap at all)
# %missed_predictions contains all missed annotations (no exon overlap at all)


my %novel_assemblies;

############################################################
# exon based comparison

# algorithm:
#
# for each cluster of transcripts
#     for each prediction in the cluster
#          cluster its exons with all the exons of the annotation
#          take the sets of consecutive exons not included in the exon clusters with annotated exons
###################

 TARGET:
    foreach my $seqname ( keys %mixed_clusters ){
	
      CLUSTER:
	foreach my $c ( @{$mixed_clusters{$seqname}} ){
	    my @predictions = grep { $type{$_} eq "prediction" } @{$c->get_Transcripts};
	    my @annotations = grep { $type{$_} eq "annotation" } @{$c->get_Transcripts};
	    my %exon_type;
	    my %exon2transcript;
	    my @all_exons;
	    my @annotated_exons;
	    
	    #print "--------------------------------------------------------------------------------\n";
	    print "Cluster with ".scalar(@predictions)." predictions and ".scalar(@annotations)." annotations\n";
	    $c->print_aligned_Cluster(\%type);
	    
	    foreach my $ann (@annotations){
		if ($verbose){
		    print "annotation:\n";
		    ClusterMerge::TranscriptUtils->_print_SimpleTranscript($ann);
		}
		foreach my $exon (@{$ann->get_all_Exons}){
		    $exon_type{$exon} = "annotation";
		    $exon2transcript{$exon} = $ann;
		    push( @annotated_exons, $exon );
		}
	    }
	    
	  PREDICTION:
	    foreach my $prediction (@predictions){
		if ($verbose){
		    print "prediction:\n";
		    ClusterMerge::TranscriptUtils->_print_SimpleTranscript($prediction);
		}
		my $id = $prediction->dbID;
		my @all_exons = @annotated_exons;
		
		my @pred_exons = sort {$a->start <=> $b->start} @{$prediction->get_all_Exons};
		my $strand = $pred_exons[0]->strand;
		print "strand=$strand\n";
		foreach my $exon ( @pred_exons ){
		    $exon2transcript{$exon} = $prediction;
		    $exon_type{$exon} = "prediction";
		    push( @all_exons, $exon );
		}
		my ($clusters,$exon2cluster)= ClusterMerge::ExonUtils->_cluster_Exons(@all_exons);
		
		my @clusters = sort { $a->start <=> $b->start } @$clusters;
		
		my %cluster_number;
		my $cluster_count = 1;
		foreach my $c ( @clusters ){
		    $cluster_number{$c} = $cluster_count;
		    $cluster_count++;
		}
		# ($cluster_count - 1) is the number of exon clusters for this transcript cluster
		
		my %exon2cluster = %{$exon2cluster};
		my @exon_assemblies;
		my $assembly;
		my $exon_count = 1;
		my %exon_number;
		my %intron_label;
		my $added_head = 0;
		my $added_tail = 0;
	      EXON:
		foreach my $exon ( @pred_exons ){
		    
		    $exon_number{$exon} = $exon_count;
		    if ( scalar( @{$exon2cluster{$exon}->get_Exons} ) <= 1 ){
			print "adding exon $exon_number{$exon}\n";
			push( @$assembly,$exon);
			$added_head = 1 if ($exon_number{$exon}==1);

			# add the first possible intersecting exon
			if ( $exon_number{$exon} > 1 && !$added_head){
			    print "head = adding exon ".($exon_number{$exon} - 1)."\n";
			    push( @$assembly, $pred_exons[ $exon_number{$exon} - 2]);
			    $added_head = 1;
			}
		    }
		    if ( scalar(  @{$exon2cluster{$exon}->get_Exons} ) > 1 && $assembly && @$assembly){
			
			# add this possible last intersecting exon
			unless ( $added_tail ){
			    push( @$assembly, $exon );
			    print "tail = adding exon $exon_number{$exon}\n";
			    $added_tail = 1;
			}
			# if the same intersecting cluster is the only one with more than one exon
			# we continue in the same assembly
			#              ###
			#  ###---###---###---###
			#
			if ( $added_tail
			     ||
			     $exon_number{$exon} == scalar(@pred_exons)
			     ||
			     !scalar( @{$exon2cluster{$pred_exons[$exon_number{$exon}]}->get_Exons} ) == 1 
			     ){
			    push( @exon_assemblies, $assembly );
			    $assembly = [];
			    $added_head = 0;
			    $added_tail = 0;
			}
		    }
		    $exon_count++;
		} # end of EXON
		# put in the last assembly:
		if ( $assembly && @$assembly){
		    push( @exon_assemblies, $assembly );
		    $assembly = [];
		}
		
		## we print out the exon-assemblies
		print scalar(@exon_assemblies)." novel exon assemblies found\n" if $verbose;
		foreach my $ass ( @exon_assemblies ){
		    ClusterMerge::ExonUtils->print_Exons($ass);
		  }
				
		# bridge:
		# within a cluster but 'bridging' (w.r.t the annotation):
		#
		#                          ###--###---### a  (annotations)
		#  ###---###                              a
		#        ###---###---###---###--###---### p  (prediction)
		#
		# they are characterized by having a common transcript between
		# any of all the previous exon-clusters and any of all the
		# the following exon clusters - To avoid a case like this:
		#
		#  ###--------------------------###---### a
		#                          ###--###---### a  (annotations)
		#  ###---###                              a
		#        ###---###---###---###--###---### p  (prediction)
		#

		# intronic assemblies:
		# 1) it can be complete (the whole predicted transcript is in an intron)
		# 2) it can be partial. It will be intronic if
		#    a) the exon assembly does not contain the first and last exons or if
		#    b) it contains either of these, they are not in the first or 
		#    c) the last exon clusters
		#
		#    ###---###---------------###     a
		#    ###---###---###---###---###     p
		# 
		#    ###---------------###---###     a
		#          ###---###---###---###     p
		#
		#    ###----###---###------###       a
		#           ###---###-###------####  p
		
		my %label;

		my %exon2ass;
	      ASSEMBLY:
		foreach my $ass ( @exon_assemblies ){
		    @$ass = sort {$a->start <=> $b->start} @$ass;
		    my $label;
		    
		    ############################################################
		    # store the exon-numbers in this assembly
		    foreach my $e ( sort {$a->start <=> $b->start } @$ass ){
			my $exonnumber;
			if ($strand == 1 ){
			    $exonnumber = $exon_number{$e};
			}
			else{
			    $exonnumber = scalar(@pred_exons) - $exon_number{$e} + 1;
			}
			$label .= $exonnumber.":";
			$exon2ass{$e} = $ass;
		    }
		    
		    ############################################################
		    # is it potentially intronic?
		    #
		    # 'intronic' means that the exon assembly are within the first and the last
		    # exons of the annotations in this transcript cluster and moreover,
		    # the assembly does not bridge non-intersecting annotations.
		    # Accordingly,  exon assembly does not contain
		    # the first and last exons from the original predicted transcripts,
		    # or if it does, they do not belong to the first or last exon-cluster,
		    # respectively.
		    #
		    # 'intronic_complete':
		    #    ###---###-----------------###     a
		    #              ###---###---###         p
		    #
		    ############################################################
		    if( ( $exon_number{$ass->[0]} != 1 
			  ||
			  ( $exon_number{$ass->[0]} == 1 
			    && 
			    $cluster_number{$exon2cluster{$ass->[0]}} != 1 
			    )
			  ||
			  ( $exon_number{$ass->[0]} == 1 
			    && 
			    $cluster_number{$exon2cluster{$ass->[0]}} == 1 
			    && 
			    scalar( @{$exon2cluster{$ass->[0]}->get_Exons} ) > 1 )
			  )
			&& 
			( $exon_number{$ass->[-1]} != scalar(@pred_exons)
			  ||
			  ( $exon_number{$ass->[-1]} == scalar(@pred_exons) 
			    && 
			    $cluster_number{$exon2cluster{$ass->[-1]}} != ($cluster_count -1 ) 
			    )
			  ||
			  ( $exon_number{$ass->[-1]} == scalar(@pred_exons) 
			    && 
			    $cluster_number{$exon2cluster{$ass->[-1]}} == ($cluster_count -1 ) 
			    &&
			    scalar( @{$exon2cluster{$ass->[-1]}->get_Exons} ) > 1 )
			  )
			){
			
			
			############################################################
			# is it actually bridging?
			# first check that there are previous and posterior exon clusters:
			
			############################################################
			# 'bridge' means that within the transcript cluster
			# the prediction bridges across two different annotations
			#
			# 'bridge_complete' means that there is bridging
			# but although there is transcript-extension overlap
			# there is no actual exon overlap, hence it can remain complete and novel:
			# 
			#                              ##--- ###----### a  (annotations)
			#  ###------###                                 a
			#        ##-----###---###---###--###----###     p  (prediction)
			#
			############################################################

			my @common_transcripts;
			if ( $cluster_number{$exon2cluster{$ass->[0]}}  >= 1
			     &&
			     $cluster_number{$exon2cluster{$ass->[-1]}} <= ($cluster_count - 1)
			     ){
			    
			    #if ($verbose){
				print "total: ".($cluster_count - 1)." clusters\n";
				print "this assembly extends from cluster  $cluster_number{$exon2cluster{$ass->[0]}} to cluster  $cluster_number{$exon2cluster{$ass->[-1]}}\n";
				#print "previous cluster = ".($cluster_number{$exon2cluster{$ass->[0]}} - 1)."\n";
				#print "next cluster     = ".($cluster_number{$exon2cluster{$ass->[-1]}} + 1)."\n";
			    #}
			    ############################################################
			    # look in all the previous exon clusters
			    my @exons_in_previous_clusters;
			    for (my $i= $cluster_number{$exon2cluster{$ass->[0]}}-1; $i>=0; $i-- ){
				push ( @exons_in_previous_clusters, @{$clusters[$i]->get_Exons} );
			    }
			    ############################################################
			    # look in all the following exon clusters
			    my @exons_in_following_clusters;
			    for (my $j= $cluster_number{$exon2cluster{$ass->[-1]}} - 1; $j< ($cluster_count - 1); $j++ ){
				push ( @exons_in_following_clusters, @{$clusters[$j]->get_Exons} );
			    }
			    
			    ############################################################
			    # find the coincident transcripts
			    my %seen;
			    foreach my $e ( @exons_in_previous_clusters ){
				next unless ($type{$exon2transcript{$e}} eq "annotation" );
				$seen{$exon2transcript{$e}} = 1;
			    }
			    foreach my $e ( @exons_in_following_clusters ){
				next unless ($type{$exon2transcript{$e}} eq "annotation" );
				$seen{$exon2transcript{$e}}++;
				if ( $seen{$exon2transcript{$e}} == 2 ){
				    push ( @common_transcripts, $exon2transcript{$e} );
				}
			    }
			}
			if ( !@common_transcripts ){
			    $label .="bridge";
			}
			else{
			    #print "not bridge as there are common transcripts:\n";
			    #foreach my $t ( @common_transcripts ){
			    #	print $t->dbID.", ";
			    #}
			    #print "\n";
			    $label .="intronic";
			}
		    }
		    ############################################################
		    # 'external' means that it is not bridging, neither intronic
		    #
		    # 
		    #                                ##----###--###   a  (annotations)
		    #  ###---###                                      a
		    #        ###--###--###      ###----###------###   p  (predictions)
		    #
		    # 'external complete' means that its extent overlaps the extent
		    # of some annotation but there is no exon overlap, and moreover,
		    # it is not bridging (it does not bridge two annotations)
		    # and it is not intronic, it has exons that are outside
		    # the introns of the annotation. eg:
		    #
		    #             ###-----###  a
		    #        ###-----###      p
		    #
		    ############################################################
		    else{
			$label .= "external";
		    }
		    ############################################################
		    # is it a complete prediction?
		    if ( scalar(@$ass) == scalar(@pred_exons) ){
			$label .="_complete";
		    }
		    
		    ############################################################
		    # print the exon-assembly in UCSC format
		    #_goldenpath_string_Exons($id,$ass,$label);
		    
		    $label{$ass} = $label;
    
		} # end of ASSEMBLY
		
		############################################################
		# for each exon add a flag with the number of bases at the 5' end which are missing
			my %fiveprime_missing;
		my @pexons = @{$prediction->get_all_Exons};
		
		# the exons from the prediction are still here: @pred_exons;
		if (  $strand == -1 ){
		    @pred_exons = sort { $b->start <=> $a->start } @pexons;
		}
		else{
		    @pred_exons = sort { $a->start <=> $b->start } @pexons;
		}
		for( my $i=0; $i<scalar( @pred_exons ); $i++ ){
		    
		    my $skipped = 0;
		    if ( $i==0){
			$fiveprime_missing{$pred_exons[$i]} = 0;
		    }
		    if ( $i>0 ){
			my $j=$i-1;
			while ( $pred_exons[$j] && !$exon2ass{$pred_exons[$j]} ){
			    $skipped += ($pred_exons[$j]->end - $pred_exons[$j]->start + 1);
			    $j--;
			}
			$fiveprime_missing{$pred_exons[$i]} = $skipped;
		    }
		}
		
		############################################################
		# we should print directly in GTF format
		foreach my $ass ( @exon_assemblies ){
		    foreach my $e ( @$ass ){
			print gtf_string( $e, $prediction, $label{$ass}, "5'missing:$fiveprime_missing{$e}" )."\n";
		    }
		}

	    

	    }     # end of PREDICTION
	}         # end of CLUSTER
    }             # end of TARGET

		

    
    

############################################################
# add whole novel predictions
#
# these are predictions that are not clustered with any annotation
#
print "WHOLE NOVEL PREDICTIONS\n" if $verbose;
foreach my $seqname ( keys %novel_predictions ){
    foreach my $t ( @{$novel_predictions{$seqname}} ){
	my $label1 = "complete";
	my $label2 = "0";
	foreach my $e ( @{$t->get_all_Exons} ){
	    print gtf_string( $e, $t, $label1, $label2 )."\n";
	}
    }
}

############################################################
# Nomenclature
############################################################
# 
# assemblies labelled as 'complete' only,
# refer to whole predictions which genomic extension does not
# overlap any annotation. Singletons from the point
# of view of clustering
#
############################################################
#
# 'intronic' means that the exon assembly are within the first and the last
# exons of the annotations in this transcript cluster and moreover,
# the assembly does not bridge non-intersecting annotations.
# Accordingly,  exon assembly does not contain
# the first and last exons from the original predicted transcripts,
# or if it does, they do not belong to the first or last exon-cluster,
# respectively.
#
# 'intronic_complete':
#
#    ###---###-----------------###     a
#              ###---###---###          p
#
############################################################
#
# 'bridging' means that within the transcript cluster
# the prediction bridges across two different annotations
#
# 'bridging_complete' means that there is bridging
# but although there is transcript-extension overlap
# there is no actual exon overlap, hence it can remain complete and novel:
# 
#                              ##--- ###----### a  (annotations)
#  ###------###                                 a
#        ##-----###---###---###--###----###     p  (prediction)
#
############################################################
#
# 'external' means that it is not bridging, neither intronic
#
# 
#                                ##----###--###   a  (annotations)
#  ###---###                                      a
#        ###--###--###      ###----###------###   p  (predictions)
#
# 'external complete' means that its extent overlaps the extent
# of some annotation but there is no exon overlap, and moreover,
# it is not bridging (it does not bridge two annotations)
# and it is not intronic, it has exons that are outside
# the introns of the annotation. eg:
#
#            ###-----###  a
#        ###-----###      p
#
############################################################


############################################################

sub exon_from_gff {
  my ($gffstring, $trans_hash) = @_;
    
  #print STDERR "gff_string: $gffstring\n";
  
  #chr6 	human_cdna	exon	1030942	1031591	100	+	0	6604
  #chr6 	human_cdna	exon	1034227	1034285	100	+	0	6604
  #chr6 	human_cdna	exon	1035280	1035339	100	+	0	6604
  my ($seqname, 
      $source, 
      $primary, 
      $start, 
      $end, 
      $score, 
      $strand, 
      $frame, 
      @group) 
      = split(/\s+/, $gffstring);
  
  $frame = 0 unless( $frame =~ /^\d+$/);
    
  my $exon = ClusterMerge::Exon->new();
  $exon->seqname($seqname);
  #$exon->source_tag($source);
  $exon->start($start);
  $exon->end($end);
  
  my $phase = ( 3 - $frame )%3;
  $exon->phase($phase);
  $exon->end_phase( ( $exon->phase + $exon->length)%3 );
  if ($score eq '.'){
      $exon->score(0);
  }
  elsif ( $score ){
      $exon->score( $score );
  }
  
  if ( $strand eq '-' ) { $exon->strand(-1); }
  elsif ( $strand eq '+' ) { $exon->strand(1); }
  elsif ( $strand eq '.' ) { $exon->strand(0); }

  ############################################################
  # warning: it parses only the first element of the group
  $trans_hash->{ $exon } = $group[0];
  return $exon;
}

############################################################


sub gtf_string{
    my ($exon,$transcript,$label1, $label2) = @_;
    my ($str,$source,$primary_tag,$score,$frame,$name,$strand);
    
    if( $exon->can('score') ) {
	$score = $exon->score();
    }
    $score = '.' unless defined $score;
    
    if( $exon->can('frame') ) {
	$frame = $exon->frame();
    }
    $frame = '.' unless defined $frame;
    
    $strand = $exon->strand();
    if(! $strand) {
	$strand = ".";
    } elsif( $strand == 1 ) {
	$strand = '+';
    } elsif ( $exon->strand == -1 ) {
	$strand = '-';
    }
    
    $name        = $exon->seqname();
    $source      = "merged";
    $primary_tag = "exon";
    
    
    # chr21	ens_and_sgp_not_twi	CDS	2251441	2251469	.	-	.	gene_id "ENSGALG00000001714xchr21_121"; transcript_id "ENSGALT00000002630xchr21_121.a";
    
    my $group_tag = "gene_id \"".$transcript->dbID."\"; ";
    $group_tag   .= "transcript_id \"".$transcript->dbID."\";";
	
    my $comments = '';
    if (defined($label1) || defined($label2) ){
	$comments .="# ";
	if ( defined($label1)){
	    $comments .= "$label1 ";
	}
	if( defined($label2) ){
	    $comments .= "$label2 ";
	}
    }

    $str = join("\t",
		$name,
		$source,
		$primary_tag,
		$exon->start(),
		$exon->end(),
		$score,
		$strand,
		$frame,
		$group_tag);

    if ( $comments ){
	$str .="\t$comments";
    }
    return $str;
}

############################################################

sub gff_string{
  my ($exon, $transcript) = @_;
  my ($str,$source,$primary_tag,$score,$frame,$name,$strand);
  
  if( $exon->can('score') ) {
    $score = $exon->score();
  }
  $score = '.' unless defined $score;
  
  if( $exon->can('frame') ) {
    $frame = $exon->frame();
  }
  $frame = '.' unless defined $frame;
  
  $strand = $exon->strand();
  if(! $strand) {
    $strand = ".";
  } elsif( $strand == 1 ) {
    $strand = '+';
  } elsif ( $exon->strand == -1 ) {
    $strand = '-';
  }
  
  $name        = $exon->seqname();
  $source      = "merged";
  $primary_tag = "exon";
  
  $str = join("\t",
	      $name,
	      $source,
	      $primary_tag,
	      $exon->start(),
	      $exon->end(),
	      $score,
	      $strand,
	      $frame);
  
  my $tag_str = $transcript->type;
  $str .= "\t$tag_str";
  return $str;
}


############################################################

sub goldenpath_string{
    my ($trans,$label) = @_;
    
    my $id = $trans->dbID;
    my @exons = @{$trans->get_all_Exons};
    _goldenpath_string_Exons($id,\@exons,$label);
}



############################################################


sub _goldenpath_string_Exons{
    my ($id,$exons,$label) = @_;
    
    # NM_014188 chr1 - 1428455 1461406 1428847 1461339 5 1428455,1430650,1431644,1451554,1461259, 1428949,1430769,1431784,1451698,1461406,
    
    my @exons = sort {$a->start <=> $b->start} @$exons;
    
    my $seqname = $exons[0]->seqname;
    my $strand  = $exons[0]->strand;
    if ( $strand == -1 ){
	$strand = "-";
    }
    else{
	$strand = "+";
    }

    my $low  = $exons[0]->start;
    my $high = $exons[-1]->end;

    my $start_string;
    my $end_string;
    foreach my $exon (@exons){
	$start_string .= $exon->start.",";
	$end_string   .= $exon->end.",";
    }
    my $exon_number = scalar(@exons);
    print "$id\t$seqname\t$strand\t$low\t$high\t$low\t$high\t$exon_number\t$start_string\t$end_string";
    if ($label){
	print "\t$label";
	}
    print "\n";
}


############################################################

sub ass_low{
    my ($ass) = @_;
    my @exons = @$ass;
    
    my $exon = shift @exons;
    my $low = $exon->start;
    foreach my $e (@exons){
	if ( $e->start < $low ){
	    $low = $e->start;
	}
    }
    return $low;
}

############################################################
