#!/usr/local/bin/perl  -w

############################################################
#
# written by Eduardo Eyras (eae@imim.es)
#
############################################################

use strict;
use Getopt::Long;
use ClusterMerge::Transcript;
use ClusterMerge::TranscriptUtils;
use ClusterMerge::TranscriptComparator;
use ClusterMerge::Exon;
use ClusterMerge::GFFTools;
use ClusterMerge::ObjectMap;
use GeneComparison::GeneComparison;
use GeneComparison::ExonComparison;
use GeneComparison::TranscriptComparison;
use GeneComparison::IntronComparison;
use GeneComparison::NucleotideComparison;


my $gff_format = 1;
my $input;
my $output;

my $short_output = 1;

my $verbose = 0;
my $verbose2 = 1;

############################################################
# There are 4 levels of comparison:
#
# gene level
# a gene is here a cluster of transcripts (according to exon-sharing)
# Sn & Sp as usual
# MG = missing genes/ real genes
# WG = overpredicted genes/ total predicted genes
# Split genes
# Joined genes
#
# transcript level
# Sn = transcripts-found/real-transcripts
# Sp = transcripts-found/predicted-transcripts
# where transcript found is each one of the paired-up transcripts
#
# MT = missing transcripts/ real transcripts
# WT = wrong transcripts/ predicted transcripts
# Split transcripts
# Joined transcripts
#
# exon-per-transcript level
# Sn & Sp of exons per each transcript pair
# ME = missing_exons/real_exons
# (not even one base overlap with the prediction)
#
# WE = wrong_exons/predicted_exons
# not even one base overlap with the annotation
#
# exon level
# Sn & Sp (of all exons)
# ME = missing_exons/real_exons
# (not even one base overlap with the prediction)
#
# WE = wrong_exons/predicted_exons
# not even one base overlap with the annotation
#
# Nucleotide-per-transcript level
# Sn & Sp of nucleotides per each transcript pair
#
# Nucleotide level
# Sn & Sp (of all found/predicted nucleotides)
#
############################################################
#
# Other two measures 'exon-connectivitity': 
# 
# How many introns are correctly predicted (at gene and transcript level)
# How many intergenic regions are correctly predicted, how many end-exon first-exon pairs 
#  are correctly found: the idea is to find how the connectivitity between exons
#  'intergenic' or 'intragenic' is found
#
############################################################


############################################################
# default options

my $pred_file = $ARGV[0];
my $ann_file  = $ARGV[1];
my $source;
my $h;
my $s;
my $gtf = 0;
my $gff = 1;
my $e;
my $p;

############################################################

&GetOptions( 
	     'h'   => \$h,                          
	     'gff' => \$gff,
	     'gtf' => \$gtf,
	     's'   => \$s,
	     'e'   => \$e,
	     'p'   => \$p,
	     );

if ( !($pred_file && $ann_file ) ||  $h ){
    print STDERR "Evaluation AT (2004)\n";
    print STDERR "program to evaluate predictions\n";
    print STDERR "taking into account alternative splicing\n";
    print STDERR "Usage: $0 predictions annotations [options]\n";
    #print STDERR "options            -gff (input in gff format - default)\n";
    #print STDERR "                   -gtf (input in gtf format)\n";
    print STDERR "                   -h   (this help)\n";
    print STDERR "                   -s   (summary results)\n";
    print STDERR "                   -e   (exact for genes and transcripts)\n";
    print STDERR "                   -p   (print details about transcript pairs)\n";
    print STDERR "\n";
    exit(0);
}

############################################################
# read predictions
############################################################
open ( IN, "<$pred_file") || die("could not open input file $pred_file");
my %tag2transcript;
my %forward_predictions;
my %reverse_predictions;
while(<IN>){
    chomp;
    my @e = split;
    next unless ( $e[3] && $e[4] ); 
    my $exon;
    #if ($gff){
    $exon = ClusterMerge::GFFTools->exon_from_gff ($_);
    #}
    #elsif($gtf){
    #$exon = ClusterMerge::GFFTools->exon_from_gtf ($_);
    #}
    $exon->source_tag("prediction");
    my $strand = $exon->strand;
    unless ( $tag2transcript{$exon->group_tag} ){
	my $trans_id = $exon->group_tag;
	$tag2transcript{$trans_id} =  ClusterMerge::Transcript->new();
	$tag2transcript{$trans_id}->dbID( $exon->group_tag );
	$tag2transcript{$trans_id}->type("prediction");
	if ($strand == 1){
	    push ( @{$forward_predictions{$exon->seqname}}, $tag2transcript{$exon->group_tag} );
	} 
	else{
	    push ( @{$reverse_predictions{$exon->seqname}}, $tag2transcript{$exon->group_tag} );
	} 
    }
    $tag2transcript{$exon->group_tag}->add_Exon($exon);
}

############################################################

############################################################
# read annotations
############################################################
open ( IN, "<$ann_file") || die("could not open input file $ann_file");
my %tag2annotation;
my %forward_annotations;
my %reverse_annotations;
while(<IN>){
    chomp;
    my @e = split;
    next unless ( $e[3] && $e[4] ); 
    my $exon;
    if ($gff){
	$exon = ClusterMerge::GFFTools->exon_from_gff ($_);
    }
    elsif($gtf){
	$exon = ClusterMerge::GFFTools->exon_from_gtf ($_);
    }
    $exon->source_tag("annotation");
    my $strand = $exon->strand;
    unless ( $tag2annotation{$exon->group_tag} ){
	my $trans_id = $exon->group_tag;
	$tag2annotation{$trans_id} =  ClusterMerge::Transcript->new();
	$tag2annotation{$trans_id}->dbID( $exon->group_tag );
	$tag2annotation{$trans_id}->type("annotation");
	if ($strand == 1){
	    push ( @{$forward_annotations{$exon->seqname}}, $tag2annotation{$exon->group_tag} );
	} 
	else{
	    push ( @{$reverse_annotations{$exon->seqname}}, $tag2annotation{$exon->group_tag} );
	} 
    }
    $tag2annotation{$exon->group_tag}->add_Exon($exon);
}
############################################################


############################################################
# Calculate TN
############################################################
# this corresponds to the sequence
# uncovered by both prediction and annotation

############################################################
# GENE LEVEL
############################################################
my %seq_names;

############################################################
# cluster predictions into genes
############################################################
my %forward_prediction_clusters;
my $predicted_genes_count;
foreach my $seqname ( keys %forward_predictions ){
    $seq_names{$seqname}=1;
    if ( $forward_predictions{$seqname} ){
	my @all_trans = @{$forward_predictions{$seqname}};
	push( @{$forward_prediction_clusters{$seqname}}, ClusterMerge::TranscriptUtils->cluster_Transcripts_into_Genes(\@all_trans,"prediction") );
	$predicted_genes_count += scalar(  @{$forward_prediction_clusters{$seqname}} );
    }
}
my %reverse_prediction_clusters;
foreach my $seqname ( keys %reverse_predictions ){
    $seq_names{$seqname}=1;
    if ( $reverse_predictions{$seqname} ){
	my @all_trans = @{$reverse_predictions{$seqname}};
	push( @{$reverse_prediction_clusters{$seqname}}, ClusterMerge::TranscriptUtils->cluster_Transcripts_into_Genes(\@all_trans,"prediction") );
	$predicted_genes_count +=  scalar(  @{$reverse_prediction_clusters{$seqname}} );
    }
}
############################################################
# cluster annotations into genes
############################################################
my %forward_annotation_clusters;
my $annotated_genes_count;
foreach my $seqname ( keys %forward_annotations ){
    if ( $forward_annotations{$seqname} ){
	$seq_names{$seqname}=1;
	my @all_trans = @{$forward_annotations{$seqname}};
	push( @{$forward_annotation_clusters{$seqname}}, 
	      ClusterMerge::TranscriptUtils->
	      cluster_Transcripts_into_Genes(\@all_trans,"annotation") );
	$annotated_genes_count += scalar( @{$forward_annotation_clusters{$seqname}} );
    }
}

my %reverse_annotation_clusters;
foreach my $seqname ( keys %reverse_annotations ){
    if ( $reverse_annotations{$seqname} ){
	$seq_names{$seqname}=1;
	my @all_trans = @{$reverse_annotations{$seqname}};
	push( @{$reverse_annotation_clusters{$seqname}}, 
	      ClusterMerge::TranscriptUtils
	      ->cluster_Transcripts_into_Genes(\@all_trans,"annotation") );
	$annotated_genes_count += scalar( @{$reverse_annotation_clusters{$seqname}} );
    }
}


############################################################
# compare the clusters:
############################################################

############################################################
# gene level
############################################################

# this counts the genes that are found by simple genomic overlap
# a stricter measure to be implemented later
my $found_exact_annotated_genes = 0;
my $total_found_annotated_genes       = 0;
my $total_found_exact_annotated_genes = 0;

# this counts the real genes with which we have no overlap at all
my $tot_missed_annotated_genes = 0;

# this counts the predicted genes that do not overlap any annotation
my $tot_overpredicted_genes    = 0;

my $total_gene_annotations = 0;
my $total_gene_predictions = 0;

my $wrong_genes  = 0;
my $missed_genes = 0;

# split_genes is the number of predictions that cover annotations
# over the number of covered annotated loci ( counting as one
# cases where the prediction covers two or more annotations )
my $split_genes  = 0;

# joined_genes is the number of annotations covered by predictions
# over the predictions that overlap (counting as one the
# predictions that overlap the same annotation
my $joined_genes = 0;

# Example:
#
#  ann genes    #######     ####   ##  ##   ##   ##
#
#  pred genes   ##  ###     ####   ######    #######
#
#  ==>    SG = 5/4   JG = 6/4

my $tot_sg_numerator   = 0;  # this counts the gene predictions for SG
my $tot_sg_denominator = 0;  # this counts the annotated loci covered
                         # by one or more predictions (the denominator of SG)
my $tot_jg_numerator   = 0;  # this counts the annotations for JG
my $tot_jg_denominator = 0;  # this counts the prediction loci
                         # covered by one or more annotations (the
                         # denominator of JG)


############################################################
# transcript level
############################################################

my $total_found_annotated_transcripts  = 0;

# this counts the transcripts which are found exactily
my $found_exact_annotated_transcripts = 0;
my $tot_found_exact_annotated_transcripts = 0;

# this counts the real transcripts to which we cannot assign a transcript prediction
my $missed_annotated_transcripts     = 0;
my $tot_missed_annotated_transcripts = 0;

# this counts the predicted transcripts to which we cannot assign an annotation
my $overpredicted_transcripts     = 0;
my $tot_overpredicted_transcripts = 0;

my $total_transcript_predictions = 0;
my $total_transcript_annotations = 0;

my $wrong_transcripts  = 0;
my $missed_transcripts = 0;
my $tot_wrong_transcripts  = 0;
my $tot_missed_transcripts = 0;

# we also want to measure the average of correct (exact) transcripts
# per gene, to know whether a raise in SNt_exact is due to
# finding more transcripts in other genes or finding more
# transcripts in the same genes (alternative transcripts)
my $average_sensitivity_exact_transcripts_per_gene = 0;
my $average_specificity_exact_transcripts_per_gene = 0;
my $average_exact_transcripts_per_gene = 0;
my $gene_counter = 0; # only counts the number of genes for which we found transcripts exactly

# this measures the annotated transcripts that correspond
# to more than one prediction
my $split_transcripts  = 0;

# this measures the predicted transcritps that correspond to
# more than one annotation
my $joined_transcripts = 0;

my $st_numerator   = 0; 
my $st_denominator = 0; 
my $jt_numerator   = 0;
my $jt_denominator = 0;
my $tot_st_numerator   = 0; 
my $tot_st_denominator = 0; 
my $tot_jt_numerator   = 0;
my $tot_jt_denominator = 0;

############################################################
# exon level
############################################################

my $tot_exon_predictions = 0;
my $tot_exon_annotations = 0;
my $tot_found_exons      = 0;
my $tot_wrong_exons      = 0;
my $tot_missing_exons    = 0;
my $total_tlevel_exon_predictions = 0;
my $total_tlevel_exon_annotations = 0;
my $total_tlevel_found_exons      = 0;
my $total_tlevel_wrong_exons      = 0;
my $total_tlevel_missing_exons    = 0;

############################################################
# nucleotide level
############################################################

my $tot_TP = 0;
my $tot_FP = 0;
my $tot_FN = 0;
my $tot_TPnt = 0;  #transcript level
my $tot_FPnt = 0;  #transcript level
my $tot_FNnt = 0;  #transcript level

############################################################
# exon connectivity (intron level)
############################################################

my $tot_intron_predictions = 0;
my $tot_intron_annotations = 0;
my $tot_found_introns      = 0;
my $tot_wrong_introns      = 0;
my $tot_missing_introns    = 0;

my $tot_tlevel_found_introns      = 0;
my $tot_tlevel_intron_annotations = 0;
my $tot_tlevel_missing_introns    = 0;
my $tot_tlevel_wrong_introns      = 0;
my $tot_tlevel_intron_predictions = 0;


## how to measure the last-exon first-exon pairs?


############################################################

my @prediction_clusters = ( \%forward_prediction_clusters, \%reverse_prediction_clusters );
my @annotation_clusters = ( \%forward_annotation_clusters, \%reverse_annotation_clusters );

############################################################
# index=0 is forward and index=1 is reverse
for (my $index = 0; $index<=1; $index++ ){
    foreach my $seqname ( keys %seq_names ){
	
	print "Looking at $seqname - index = $index\n" if $verbose2;
	my %prediction_clusters = %{$prediction_clusters[$index]};
	my %annotation_clusters = %{$annotation_clusters[$index]};
	my @clusters;
	if ( $prediction_clusters{$seqname} && @{$prediction_clusters{$seqname}} ){
	    push( @clusters,  @{$prediction_clusters{$seqname}} );
	}
	if ( $annotation_clusters{$seqname} && @{$annotation_clusters{$seqname}} ){
	    push( @clusters, @{$annotation_clusters{$seqname}} );
	}
	
	#print "eval: clusters = ".scalar(@clusters)."\n";
	next unless (@clusters);
	############################################################
	# this call compares the genomic extent of the prediction and annotation
	# clusters (genes) and groups together those that overlap
	my @clustered_clusters = 
	    GeneComparison::GeneComparison->cluster_TranscriptClusters( \@clusters );
	#print "eval: gene clusters = ".scalar(@clustered_clusters)."\n";
	
	############################################################
	# GENE LEVEL
	############################################################
	
	############################################################
	# each $c is a cluster of "ranges", a GeneComparison::RangeCluster. 
	# and represent a cluster of genes (predictions and annotations).
	# Each gene is represented by a TranscriptCluster.
	# The ranges here are transcript clusters (i.e., genes)
	
	# each $c is potentially the comparison of one prediction
	# with one annotation
	my $gene_count = 0;
      GENE_LEVEL:
	foreach my $c ( @clustered_clusters ){
	    $gene_count++;
	    print "Looking at gene: $gene_count\n" if $verbose2;
	    my $sg_numerator   = 0;  # this counts the gene predictions for SG
	    my $sg_denominator = 0;  # this counts the annotated loci covered 
	    # by one or more predictions (the denominator of SG)
	    my $jg_numerator   = 0;  # this counts the annotations for JG
	    my $jg_denominator = 0;  # this counts the prediction loci
	    # covered by one or more annotations (the denominator of JG)
	    
	    my $overpredicted_genes = 0;
	    my $missed_annotated_genes = 0;
	    $found_exact_annotated_genes = 0;

	    ############################################################
	    # @predictions and @annotations are lists of ClusterMerge::TranscriptCluster objects
	    # each element corresponding to one gene
	    my @predictions = grep { $_->type eq "prediction" } @{$c->get_Ranges};
	    my @annotations = grep { $_->type eq "annotation" } @{$c->get_Ranges};
	    
	    print "predictions: ".scalar(@predictions)." annotations: ".scalar(@annotations)."\n";
	    my $predicted_genes;
	    my $annotated_genes;
	    
	    ############################################################
	    # Exon comparison at the Gene level
	    #print "Comparing exons at gene level\n";
	    my ($exon_predictions,
		$exon_annotations,
		$found_exons,
		$wrong_exons,
		$missing_exons) = GeneComparison::ExonComparison->compare_exons_Gene_level(\@predictions,\@annotations);
	    
	    print "exon_predictions: $exon_predictions\n";
	    print "exon_annotations: $exon_annotations\n";
	    print "found_exons     : $found_exons\n";
	    print "wrong           : $wrong_exons\n";
	    print "missed          : $missing_exons\n";
	    

	    $tot_exon_predictions +=  $exon_predictions;
	    $tot_exon_annotations +=  $exon_annotations;
	    $tot_found_exons      +=  $found_exons;
	    $tot_wrong_exons      +=  $wrong_exons;
	    $tot_missing_exons    +=  $missing_exons;
	    ############################################################
	    
	    ############################################################
	    # Nucleotide comparison at the Gene level
	    my ($TP,$FP,$FN) = GeneComparison::NucleotideComparison->compare_nucleotides_Gene_level(\@predictions,\@annotations);
	    $tot_TP +=  $TP;
	    $tot_FP +=  $FP;
	    $tot_FN +=  $FN;
	    #print "nucleotide gene level\n";
	    #print "TP=$TP --> tot_TP=$tot_TP\n";
	    #print "FP=$FP --> tot_FP=$tot_FP\n";
	    #print "FN=$FN --> tot_FN=$tot_FN\n";
	    ############################################################

	    ############################################################
	    # Intron comparison at Gene level (exon connectivity at the gene level)
	    #print "Comparing introns at gene level\n";
	    my ($intron_predictions,
		$intron_annotations,
		$found_introns,
		$wrong_introns,
		$missing_introns) = GeneComparison::IntronComparison->compare_introns_Gene_level(\@predictions,\@annotations);
	    
	    $tot_intron_predictions +=  $intron_predictions;
	    $tot_intron_annotations +=  $intron_annotations;
	    $tot_found_introns      +=  $found_introns;
	    $tot_wrong_introns      +=  $wrong_introns;
	    $tot_missing_introns    +=  $missing_introns;
	    ############################################################

	    ############################################################
	    my $thisgene_found_annotated_transcripts       = 0;
	    my $thisgene_found_exact_annotated_transcripts = 0;
	    my $thisgene_missed_annotated_transcripts      = 0;
	    my $thisgene_overpredicted_transcripts         = 0;
	    my $thisgene_annotated_transcripts             = 0;
	    my $thisgene_predicted_transcripts             = 0;
	    
	    my $thisgene_tlevel_exon_predictions           = 0;
	    my $thisgene_tlevel_exon_annotations           = 0;
	    my $thisgene_tlevel_found_exons                = 0;
	    my $thisgene_tlevel_wrong_exons                = 0;
	    my $thisgene_tlevel_missing_exons              = 0;
	    
	    my $thisgene_tlevel_found_introns      = 0;
	    my $thisgene_tlevel_intron_annotations = 0;
	    my $thisgene_tlevel_missing_introns    = 0;
	    my $thisgene_tlevel_wrong_introns      = 0;
	    my $thisgene_tlevel_intron_predictions = 0;

	    my $thisgene_tlevel_TPn = 0;
	    my $thisgene_tlevel_FPn = 0;
	    my $thisgene_tlevel_FNn = 0;

	    #my $split_transcripts  = 0;
	    #my $joined_transcripts = 0;
	    #my $st_numerator   = 0; 
	    #my $st_denominator = 0; 
	    #my $jt_numerator   = 0;
	    #my $jt_denominator = 0;
	    
	    my $found_annotated_genes       = 0;

	    ############################################################
	    # if this cluster only contains predictions
	    if ( !@annotations && @predictions){
		my $temp = scalar(@predictions);
		$overpredicted_genes     += $temp;
		$tot_overpredicted_genes += $temp;
		$predicted_genes          = $temp;
		$total_gene_predictions  += $temp;
		
		# count the transcripts
		my @transcript_predictions;
		foreach my $gene_prediction ( @predictions ){
		    push( @transcript_predictions, @{$gene_prediction->get_Transcripts} );
		}
		$thisgene_predicted_transcripts     = scalar(@transcript_predictions);
		$thisgene_overpredicted_transcripts = $thisgene_predicted_transcripts;
		$total_transcript_predictions      += $thisgene_predicted_transcripts;
	    }
	    ############################################################
	    # if this cluster only contains annotations
	    elsif ( @annotations && !@predictions){
		my $temp = scalar(@annotations);
		$tot_missed_annotated_genes += $temp;
		$missed_annotated_genes     += $temp;
		$annotated_genes             = $temp;
		$total_gene_annotations     += $temp;

		# count the transcripts
		my @transcript_annotations;
		foreach my $gene_annotation ( @annotations ){
		    push( @transcript_annotations, @{$gene_annotation->get_Transcripts} );
		}
		$thisgene_annotated_transcripts        = scalar(@transcript_annotations);		
		$thisgene_missed_annotated_transcripts = $thisgene_annotated_transcripts;
		$total_transcript_annotations         += $thisgene_annotated_transcripts;
	    }
	    ############################################################
	    # if this cluster contains both predictions and annotations
	    elsif(  @annotations && @predictions ){
		
		############################################################
		# this is a cluster with both gene predictions and gene annotations:
		my $temp1 = scalar(@predictions);
		my $temp2 = scalar(@annotations);
		$total_gene_predictions     += $temp1;
		$total_gene_annotations     += $temp2;
		$total_found_annotated_genes+= $temp2;
		$found_annotated_genes       = $temp2;
		$annotated_genes             = $temp2;
		$predicted_genes             = $temp1;
		if ( @annotations == 1 && @predictions == 1 ){
		    $tot_sg_numerator++;
		    $tot_sg_denominator++;
		    $tot_jg_numerator++;
		    $tot_jg_denominator++;
		    ($sg_numerator, $sg_denominator, $jg_numerator, $jg_denominator) = (1,1,1,1);
		}
		# if we have multiple predictions/annotations for one cluster
		# we have split/joined genes
		elsif ( @annotations == 1 && @predictions > 1 ){
		    # split genes
		    $tot_sg_numerator += scalar(@predictions);
		    $tot_sg_denominator++;
		    $tot_jg_numerator++;
		    $tot_jg_denominator++;
		    ($sg_numerator,$sg_denominator,$jg_numerator,$jg_denominator) = (scalar(@predictions),1,1,1);
		}
		elsif ( @annotations > 1 && @predictions == 1 ){
		    # joined genes
		    $tot_sg_numerator++;
		    $tot_sg_denominator++;
		    $tot_jg_numerator += scalar(@annotations);
		    $tot_jg_denominator++;
		    ($sg_numerator,$sg_denominator,$jg_numerator,$jg_denominator) = (1,1,scalar(@annotations),1);
		}
		elsif ( @annotations > 1 && @predictions > 1 ){
		    # we approximate this situation by:
		    $tot_sg_numerator += scalar(@predictions);
		    $tot_sg_denominator++;
		    $tot_jg_numerator += scalar(@annotations);
		    $tot_jg_denominator++;
		    ( $sg_numerator,$sg_denominator,$jg_numerator,$jg_denominator) = (scalar(@predictions),1,scalar(@annotations),1);
		}
		else{
		    # we should not be here
		    print "ERROR\n";
		    exit(0);
		}
		
		############################################################
		# TRANSCRIPT LEVEL
		############################################################
		# transcript level measures:
		#
		# Sn = transcripts-found/real-transcripts
		# Sp = transcripts-found/predicted-transcripts
		#
		# MT = missing transcripts/ real transcripts
		# WT = wrong transcripts/ predicted transcripts
		#
		# Split transcripts
		# Joined transcripts
		
		# how are the split and joined transcripts calculated?
		# 
		
		# how to choose which transcript goes with which? best CCn?,
		my @transcript_predictions;
		foreach my $gene_prediction ( @predictions ){
		    push( @transcript_predictions, @{$gene_prediction->get_Transcripts} );
		}
		$thisgene_predicted_transcripts = scalar(@transcript_predictions);
		my @transcript_annotations;
		foreach my $gene_annotation ( @annotations ){
		    push( @transcript_annotations, @{$gene_annotation->get_Transcripts} );
		}
		$thisgene_annotated_transcripts = scalar(@transcript_annotations);
		
		# we recover a list of pairs, 
		# where each element of the list will have
		# the transcript pair and the 
		# scores: [prediction, annotation, SNn, SPn, CCn]

		my @transcript_pairs = GeneComparison::TranscriptComparison->pair_Transcripts(\@transcript_predictions,\@transcript_annotations);
		
		# To be able to look at split/joined transcripts
		# we have to match transcripts not just in pairs
		# but also to allow 1-to-many and many-to-many matchings
		# WHAT IS THE APPROPRIATE ALGORITHM?
		# idea: 
		# for each found pair:
		#   given the annotation
		#      look at the rest of the predictions and
		#      see whether there is one not already in a pair
		#      (perhaps not above some threshold) that can 'complete'
		#      the annotation and does not overlap with the first prediction
		#
		#  given the prediction in the pair, repeat 
		#  the search with the annotations
		#
		#  (each time a pair is broken we 
		#  eliminate it and subsume it in the
		#  better match.)
		#  Continue until you have exhausted the pairs
		
		# this counts the transcripts that are 
		# found by best-exon-coverage
		
		$total_transcript_annotations += scalar(@transcript_annotations);
		$total_transcript_predictions += scalar(@transcript_predictions);
		
		#this is counted below 
		#$total_found_annotated_transcripts += scalar(@transcript_pairs);
		$thisgene_missed_annotated_transcripts  =
                    scalar(@transcript_annotations) - scalar(@transcript_pairs);
		$missed_annotated_transcripts += 
		    scalar(@transcript_annotations) - scalar(@transcript_pairs);
		$thisgene_overpredicted_transcripts =
                    scalar(@transcript_predictions) - scalar(@transcript_pairs);
		$overpredicted_transcripts    += 
		    scalar(@transcript_predictions) - scalar(@transcript_pairs);
		
		# Now for each transcript pair found we can calculate the rest of the 'per-transcript'
		foreach my $pair ( @transcript_pairs ){
		    my $pred_tran = $pair->[0];
		    my $ann_tran  = $pair->[1];
		    my $CCn       = $pair->[2];
		    my $TPn       = $pair->[3];
		    my $FPn       = $pair->[4];
		    my $FNn       = $pair->[5];
		    
		    print "TRANS-PAIR\t".$pred_tran->dbID."\t".$ann_tran->dbID."\n";

		    ############################################################
		    # NUCLEOTIDE PER-TRANSCRIPT-PAIR LEVEL
		    ############################################################
		    # Nucleotide-per-transcript level measures:
		    # Sn & Sp of nucleotides per each transcript pair
		    $tot_TPnt += $TPn;
		    $tot_FPnt += $FPn;
		    $tot_FNnt += $FNn;
		    if ($verbose){
			print "adding to thisgene_tlevel:\n";
			print "thisgene_tlevel_TPn = $thisgene_tlevel_TPn <-- $TPn\n";
			print "thisgene_tlevel_FPn = $thisgene_tlevel_FPn <-- $FPn\n";
			print "thisgene_tlevel_FNn = $thisgene_tlevel_FNn <-- $FNn\n";
		    }
		    $thisgene_tlevel_TPn += $TPn;
		    $thisgene_tlevel_FPn += $FPn;
		    $thisgene_tlevel_FNn += $FNn;
		    		    
		    ############################################################
		    # INTRON PER-TRANSCRIPT-PAIR COMPARISON
		    ############################################################
		    # if this is better than the equivalent exonic
		    # comparisons, it means that we do find the correct introns
		    # even though the exon-boundaries are not correct
		    my ($thispair_intron_predictions,
			$thispair_intron_annotations,
			$thispair_found_introns,
			$thispair_wrong_introns,
			$thispair_missing_introns) = GeneComparison::IntronComparison->compare_introns_Transcript_level($pred_tran,$ann_tran);
		    
		    $thisgene_tlevel_found_introns      += $thispair_found_introns;
		    $thisgene_tlevel_intron_annotations += $thispair_intron_annotations;
		    $thisgene_tlevel_missing_introns    += $thispair_missing_introns;
		    $thisgene_tlevel_wrong_introns      += $thispair_wrong_introns;
		    $thisgene_tlevel_intron_predictions += $thispair_intron_predictions;

		    $tot_tlevel_intron_predictions      +=  $thispair_intron_predictions;
		    $tot_tlevel_intron_annotations      +=  $thispair_intron_annotations;
		    $tot_tlevel_found_introns           +=  $thispair_found_introns;
		    $tot_tlevel_wrong_introns           +=  $thispair_wrong_introns;
		    $tot_tlevel_missing_introns         +=  $thispair_missing_introns;
		    
		    ############################################################
		    # EXON PER-TRANSCRIPT-PAIR LEVEL
		    ############################################################
		    # exon-per-transcript level measures:
		    # Sn & Sp of exons per each transcript pair
		    # ME = missing_exons/real_exons (not even one base overlap with the prediction)
		    # WE = wrong_exons/predicted_exons (not even one base overlap with the annotation)
		    #print "eval: about to call compare_Exons (2)\n";
		    my ($thispair_exon_predictions,
			$thispair_exon_annotations,
			$thispair_found_exons,
			$thispair_wrong_exons,
			$thispair_missing_exons) 
			= GeneComparison::ExonComparison->compare_Exons($pred_tran->get_all_Exons,$ann_tran->get_all_Exons);
		    
		    $thisgene_tlevel_exon_predictions += $thispair_exon_predictions;
		    $thisgene_tlevel_exon_annotations += $thispair_exon_annotations;
		    $thisgene_tlevel_found_exons      += $thispair_found_exons;
		    $thisgene_tlevel_wrong_exons      += $thispair_wrong_exons;
		    $thisgene_tlevel_missing_exons    += $thispair_missing_exons;

		    $total_tlevel_exon_predictions    += $thispair_exon_predictions;
		    $total_tlevel_exon_annotations    += $thispair_exon_annotations;
		    $total_tlevel_found_exons         += $thispair_found_exons;
		    $total_tlevel_wrong_exons         += $thispair_wrong_exons;
		    $total_tlevel_missing_exons       += $thispair_missing_exons;
		    
		    $thisgene_found_annotated_transcripts++;
		    $total_found_annotated_transcripts++;
		    if ( GeneComparison::TranscriptComparison->exact_match( $pred_tran, $ann_tran ) ){
			$thisgene_found_exact_annotated_transcripts++;
			$found_exact_annotated_transcripts++;
		    }
		    
		    if ($p){
			print_TranscriptPair($pair);
		    }
		    
		}   # end of TRANSCRIPT-PAIR loop
		
		# if in this gene cluster we find exactily all the
		# the transcripts, we consider the gene(s) as exactly found
		if ( $thisgene_found_exact_annotated_transcripts == scalar(@transcript_annotations )){
		    $total_found_exact_annotated_genes += $annotated_genes;
		    $found_exact_annotated_genes = $annotated_genes;
		}

	    }  # end of the case where we have both Annotations and Predictions
	    
	    ############################################################
	    # GENE
	    # print the results for this gene/gene-pair unless we only want the summary
	    # ($SNg,$SPg,$SG,$JG,$SNt,$SPt,$ST,$JT,$SNe,$SPe,$SNet,$SPet,$SNn,$SPn,$SNnt,$SPnt)
	    unless ($s){
		
		############################################################
		# gene comparison
		############################################################
		my ($SNg,$SPg,$SNg_exact,$SPg_exact) = (0,0,0,0);
		if ( $annotated_genes ){
		    print "found annotated genes: $found_annotated_genes\n";
		    print "annotated genes       : $annotated_genes\n";
		    $SNg       = calculate($found_annotated_genes,$annotated_genes);
		    $SNg_exact = calculate($found_exact_annotated_genes,$annotated_genes);
		}
		if ( $predicted_genes ){
		    print "found annotated genes: $found_annotated_genes\n";
		    print "predicted genes       : $predicted_genes\n";
		    $SPg       = calculate($found_annotated_genes,$predicted_genes);
		    $SPg_exact = calculate($found_exact_annotated_genes,$predicted_genes);
		}
		my $WG = $overpredicted_genes;
		my $MG = $missed_annotated_genes;
		
		my ($SG,$JG) = (0,0);
		if ( $jg_denominator ){
		    $JG = calculate($jg_numerator,$jg_denominator);
		}
		if ( $sg_denominator ){
		    $SG = calculate($sg_numerator,$sg_denominator);
		}
		############################################################
		# transcript comparison
		############################################################
		my ( $SNt,$SPt,$SNt_exact,$SPt_exact,$WT,$MT,$ST,$JT) =
		    ( 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 );
		if ( $thisgene_annotated_transcripts ){
		    $SNt       = calculate($thisgene_found_annotated_transcripts,$thisgene_annotated_transcripts);
		    $SNt_exact = calculate($thisgene_found_exact_annotated_transcripts,$thisgene_annotated_transcripts);
		    $MT        = calculate($thisgene_missed_annotated_transcripts,$thisgene_annotated_transcripts);
		}
		if ($thisgene_predicted_transcripts ){
		    $SPt       = calculate($thisgene_found_annotated_transcripts,$thisgene_predicted_transcripts);
		    $SPt_exact = calculate($thisgene_found_exact_annotated_transcripts,$thisgene_predicted_transcripts); 
		    $WT        = calculate($thisgene_overpredicted_transcripts,$thisgene_predicted_transcripts);
		}
		if ( $SNt_exact && $SPt_exact ){
		    $average_sensitivity_exact_transcripts_per_gene += $SNt_exact;
		    $average_specificity_exact_transcripts_per_gene += $SPt_exact;
		    $average_exact_transcripts_per_gene             += $thisgene_found_exact_annotated_transcripts;

		    # only counts the number of genes for which we found transcripts exactly
		    $gene_counter++; 
		}

		# transcripts splits/joins
		$ST = 0;
		$JT = 0;
		
		############################################################
		# exon comparison
		############################################################
		my ($SNe,$SPe,$WE,$ME) = (0,0,0,0);
		if ( $exon_annotations){
		    $SNe  = calculate($found_exons,$exon_annotations);   # gene-level comparison
		    $ME   = calculate($missing_exons,$exon_annotations);
		}
		if ( $exon_predictions ){
		    $SPe  = calculate($found_exons,$exon_predictions);   # gene-level comparison
		    $WE   = calculate($wrong_exons,$exon_predictions);
		}
		
		my ($SNet,$SPet,$WEt,$MEt) = (0,0,0,0);
		if ( $thisgene_tlevel_exon_annotations ){
		    $SNet = calculate($thisgene_tlevel_found_exons,$thisgene_tlevel_exon_annotations);
		    $MEt  = calculate($thisgene_tlevel_missing_exons,$thisgene_tlevel_exon_annotations);
		}
		if ( $thisgene_tlevel_exon_predictions ){
		    $SPet = calculate($thisgene_tlevel_found_exons,$thisgene_tlevel_exon_predictions);
		    $WEt  = calculate($thisgene_tlevel_wrong_exons,$thisgene_tlevel_exon_predictions);
		}

		############################################################
		# intron comparison  (gene-level comparison)
		############################################################
		my ($SNi,$SPi,$WI,$MI) = (0,0,0,0);
		if ( $intron_annotations ){
		    $SNi = calculate($found_introns,$intron_annotations);   
		    $MI  = calculate($missing_introns,$intron_annotations);
		}
		if ( $intron_predictions ){
		    $SPi = calculate($found_introns,$intron_predictions);   
		    $WI  = calculate($wrong_introns,$intron_predictions);
		}
		
		my ($SNit,$SPit,$WIt,$MIt) = (0,0,0,0);


		############################################################
		# nucleotide comparison
		############################################################
		my ($SNn,$SPn,$SNnt,$SPnt) = (0,0,0,0);
		if ($TP || ( $FN && $FP ) ){
		    $SNn  = calculate($TP,($TP+$FN)); #gene-level
		    $SPn  = calculate($TP,($TP+$FP)); #gene-level
		}
		if( $thisgene_tlevel_TPn || ( $thisgene_tlevel_FNn && $thisgene_tlevel_FPn ) ){
		    $SNnt = calculate($thisgene_tlevel_TPn,($thisgene_tlevel_TPn+$thisgene_tlevel_FNn));
		    $SPnt = calculate($thisgene_tlevel_TPn,($thisgene_tlevel_TPn+$thisgene_tlevel_FPn));
		}

		############################################################
		print_evaluation($SNg,$SPg,$SNg_exact,$SPg_exact,$WG,$MG,$JG,$SG,
				 $SNt,$SPt,$SNt_exact,$SPt_exact,$WT,$MT,$ST,$JT,
				 $SNe,$SPe,$WE,$ME,$SNet,$SPet,$WEt,$MEt,
				 $SNi,$SPi,$WI,$MI,$SNit,$SPit,$WIt,$MIt,
				 $SNn,$SPn,$SNnt,$SPnt,\@predictions,\@annotations);
		
	    } # end of unless ($s)
	} # end of GENE_LEVEL
    }     # end of loop over target names
}         # end of loop over strands


############################################################
#
# totals
#
############################################################

############################################################
# gene comparison
############################################################
my $TSNg       = calculate($total_found_annotated_genes,$total_gene_annotations);
my $TSPg       = calculate($total_found_annotated_genes,$total_gene_predictions);
my $TSNg_exact = calculate($total_found_exact_annotated_genes,$total_gene_annotations);
my $TSPg_exact = calculate($total_found_exact_annotated_genes,$total_gene_predictions);
my $TWG        = calculate($tot_overpredicted_genes,$total_gene_predictions);
my $TMG        = calculate($tot_missed_annotated_genes,$total_gene_annotations);
my $TJG        = 0;
my $TSG        = 0;
$TJG           = calculate($tot_jg_numerator,$tot_jg_denominator) if ( $tot_jg_denominator && $tot_jg_numerator );
$TSG           = calculate($tot_sg_numerator,$tot_sg_denominator) if ( $tot_sg_denominator && $tot_sg_numerator );

############################################################
# transcript comparison
############################################################
my $TSNt       = calculate($total_found_annotated_transcripts,$total_transcript_annotations);
my $TSPt       = calculate($total_found_annotated_transcripts,$total_transcript_predictions);
my $TSNt_exact = calculate($found_exact_annotated_transcripts,$total_transcript_annotations);
my $TSPt_exact = calculate($found_exact_annotated_transcripts,$total_transcript_predictions);
my $TWT        = calculate($overpredicted_transcripts,$total_transcript_predictions);
my $TMT        = calculate($missed_annotated_transcripts,$total_transcript_annotations);
my $TST        = 0;
my $TJT        =  0;
$TST           = calculate($tot_st_numerator,$tot_st_denominator) if ( $tot_st_denominator && $tot_st_numerator );
$TJT           = calculate($tot_jt_numerator,$tot_jt_denominator) if ( $tot_jt_denominator && $tot_jt_numerator );

if ($gene_counter){   
    $average_sensitivity_exact_transcripts_per_gene = calculate($average_sensitivity_exact_transcripts_per_gene,$gene_counter);
    $average_specificity_exact_transcripts_per_gene = calculate($average_specificity_exact_transcripts_per_gene,$gene_counter);
    $average_exact_transcripts_per_gene             = calculate($average_exact_transcripts_per_gene,$gene_counter);
}


############################################################
# exon comparison
############################################################
my $TSPe  = calculate($tot_found_exons,$tot_exon_predictions);   # gene-level comparison
my $TSNe  = calculate($tot_found_exons,$tot_exon_annotations);   # gene-level comparison
my $TWE   = calculate($tot_wrong_exons,$tot_exon_predictions);
my $TME   = calculate($tot_missing_exons,$tot_exon_annotations);

my ($TSNet,$TMEt,$TSPet,$TWEt ) = (0,0,0,0);
if ($total_tlevel_exon_annotations){
    $TSNet = calculate($total_tlevel_found_exons,$total_tlevel_exon_annotations);
    $TMEt  = calculate($total_tlevel_missing_exons,$total_tlevel_exon_annotations);
}  
if ($total_tlevel_exon_predictions){
    $TSPet = calculate($total_tlevel_found_exons,$total_tlevel_exon_predictions);
    $TWEt  = calculate($total_tlevel_wrong_exons,$total_tlevel_exon_predictions);
}

############################################################
# intron comparison
############################################################
my ($TSNi,$TSPi,$TWI,$TMI) = (0,0,0,0);
if ( $tot_intron_annotations){
    $TSNi = calculate($tot_found_introns,$tot_intron_annotations);
    $TMI  = calculate($tot_missing_introns,$tot_intron_annotations);
}
if($tot_intron_predictions){
    $TSPi = calculate($tot_found_introns,$tot_intron_predictions);
    $TWI  = calculate($tot_wrong_introns,$tot_intron_predictions);
}
my ($TSNit,$TSPit,$TWIt,$TMIt) = (0,0,0,0);
if ( $tot_tlevel_intron_annotations ){
    $TSNit = calculate($tot_tlevel_found_introns  ,$tot_tlevel_intron_annotations);
    $TMIt  = calculate($tot_tlevel_missing_introns,$tot_tlevel_intron_annotations);
}
if ( $tot_tlevel_intron_predictions ){
    $TSPit = calculate($tot_tlevel_found_introns  ,$tot_tlevel_intron_predictions);
    $TWIt  = calculate($tot_tlevel_wrong_introns  ,$tot_tlevel_intron_predictions);
}

############################################################
# nucleotide comparison
############################################################
my $TSNn  = calculate($tot_TP,($tot_TP+$tot_FN)); #gene-level
my $TSPn  = calculate($tot_TP,($tot_TP+$tot_FP)); #gene-level

my $TSNnt = calculate($tot_TPnt,($tot_TPnt+$tot_FNnt));
my $TSPnt = calculate($tot_TPnt,($tot_TPnt+$tot_FPnt));

print "Total:\n";
my @summary 
    = ("SUMMARY\tGene-ann:$total_gene_annotations\tGene-pre:$total_gene_predictions\tGene-fou:$total_found_annotated_genes",
       "SUMMARY\tTran-ann:$total_transcript_annotations\tTran-pre:$total_transcript_predictions\tTran-fou:$found_exact_annotated_transcripts",
       "SUMMARY\tAver-SNt:$average_sensitivity_exact_transcripts_per_gene\t".
       "Aver-SPt:$average_specificity_exact_transcripts_per_gene\t".
       "Aver-exacttrans-per-gene:$average_exact_transcripts_per_gene",  
       "SUMMARY\tExon-ann:$tot_exon_annotations\tExon-pre:$tot_exon_predictions\tExon-fou:$tot_found_exons",
       "SUMMARY\tIntr-ann:$tot_intron_annotations\tIntr-pre:$tot_intron_predictions\tIntr-fou:$tot_found_introns",
       "SUMMARY\tNucl-ann:".($tot_TP+$tot_FN)."\tNucl-pre:".($tot_TP+$tot_FP)."\tNucl-fou:$tot_TP"
       );
my $summary = join "\n", @summary;
print $summary."\n";

print "SUMMARY\tTRAN-LEVEL\tnt-ann:".($tot_TPnt+$tot_FNnt)."\tnt-pred:".($tot_TPnt+$tot_FPnt)."\tnt-found:".$tot_TPnt."\n";
print "SUMMARY\tTRAN-LEVEL\tex-ann:".($total_tlevel_exon_annotations)."\tex-pred:".($total_tlevel_exon_predictions)."\tex-found:".($total_tlevel_found_exons)."\n";

print_evaluation($TSNg,$TSPg,$TSNg_exact,$TSPg_exact,$TWG,$TMG,$TJG,$TSG,
		 $TSNt,$TSPt,$TSNt_exact,$TSPt_exact,$TWT,$TMT,$TST,$TJT,
		 $TSNe,$TSPe,$TWE,$TME,$TSNet,$TSPet,$TWEt,$TMEt,
		 $TSNi,$TSPi,$TWI,$TMI,$TSNit,$TSPit,$TWIt,$TMIt,
		 $TSNn,$TSPn,$TSNnt,$TSPnt);
    

############################################################
#
# method to calculate the round-p the fractions to two decimals
#
############################################################
sub calculate{
    my ($num,$denom) = @_;
    return int(100*$num/$denom + 0.5)/100;
}

############################################################
#
# method to print the evaluation results
# it is used for both the per-gene output and the totals
#
############################################################
sub print_evaluation{
    my ($SNg,$SPg,$SNg_exact,$SPg_exact,$WG,$MG,$JG,$SG,
	$SNt,$SPt,$SNt_exact,$SPt_exact,$WT,$MT,$ST,$JT,
	$SNe,$SPe,$WE,$ME,$SNet,$SPet,$WEt,$MEt,
	$SNi,$SPi,$WI,$MI,$SNit,$SPit,$WIt,$MIt,
	$SNn,$SPn,$SNnt,$SPnt,$predictions,$annotations) = @_;

    #$SNg
    #$SPg
    #$SNg_exact
    #$SPg_exact
    #$WG
    #$MG
    #$JG
    #$SG

    #$SNt
    #$SPt
    #$SNt_exact
    #$SPt_exact
    #$WT
    #$MT
    #$ST
    #$JT

    #$SNe
    #$SPe
    #$WE
    #$ME
    #$SNet
    #$SPet
    #$WEt
    #$MEt

    #$SNi
    #$SPi
    #$WI
    #$MI
    #$SNit
    #$SPit
    #$WIt
    #$MIt,
    
    #$SNn
    #$SPn
    #$SNnt
    #$SPnt
    #$predictions
    #$annotations

    # we extract all the transcripts from a pair
    # of gene-cluster comparisons:
    my $pred_label = "pred:\t";
    my $ann_label  = "ann:\t";

    my @transcript_predictions;
    if ( $predictions && @$predictions ){
	foreach my $gene_prediction ( @$predictions ){
	    push( @transcript_predictions, @{$gene_prediction->get_Transcripts} );
	}
    }
    my @transcript_annotations;
    if ($annotations && @$annotations){
	foreach my $gene_annotation ( @$annotations ){
	    push( @transcript_annotations, @{$gene_annotation->get_Transcripts} );
	}
    }
    if ( @transcript_predictions ){
	foreach my $t ( @transcript_predictions ){
	    $pred_label .= $t->dbID.",";
	}
    }
    else{
	$pred_label .= "none";
    }
    if ( @transcript_annotations ){
	foreach my $t (@transcript_annotations ){
	    $ann_label .= $t->dbID.",";
	}
    }
    else{
	$ann_label .= "none";
    }
    
    ############################################################
    # output
    my $header;
    my $string;
    
    if ($short_output){
	$header =
	    join "\t", ("LABEL","SNg","SPg","WG","MG",
			"SNt","SPt","WT","MT",
			"SNe","SPe","WE","ME","SNet","SPet","WEt","MEt",
			"SNi","SPi","WI","MI","SNit","SPit","WIt","TMIt",
			"SNn","SPn","SNnt","SPnt");
	$string = 
	    join "\t", ("SCORE",
			$SNg,$SPg,$WG,$MG,
			$SNt_exact,$SPt_exact,$WT,$MT,
			$SNe,$SPe,$WE,$ME,
			$SNet,$SPet,$WEt,$MEt,
			$SNi,$SPi,$WI,$MI,
			$SNit,$SPit,$WIt,$MIt,
			$SNn,$SPn,$SNnt,$SPnt);
    }
    else{
	############################################################
	# exact measures
	if ( $e ){
	    $header =
		join "\t", ("LABEL","SNg_exact","SPg_exact","WG","MG","JG","SG",
			    "SNt_exact","SPt_exact","WT","MT","ST","JT",
			    "SNe","SPe","WE","ME","SNet","SPet","WEt","MEt",
			    "SNi","SPi","WI","MI","SNit","SPit","WIt","TMIt",
			    "SNn","SPn","SNnt","SPnt");
	    $string = 
		join "\t", ("SCORE",$SNg_exact, $SPg_exact, $WG, $MG, $JG, $SG, 
			    $SNt_exact, $SPt_exact, $WT,$MT,$ST,$JT,
			    $SNe,$SPe,$WE,$ME,$SNet,$SPet,$WEt,$MEt,
			    $SNi,$SPi,$WI,$MI,$SNit,$SPit,$WIt,$MIt,
			    $SNn,$SPn,$SNnt,$SPnt);
	}
	############################################################
	# approximate measures
	else{
	    $header =
		join "\t", ("LABEL","SNg","SPg","WG","MG","JG","SG",
			    "SNt","SPt","WT","MT","ST","JT",
			    "SNe","SPe","WE","ME","SNet","SPet","WEt","MEt",
			    "SNi","SPi","WI","MI","SNit","SPit","WIt","TMIt",
			    "SNn","SPn","SNnt","SPnt");
	    $string = 
		join "\t", ("SCORE",$SNg,$SPg,$WG,$MG,$JG,$SG,
			    $SNt,$SPt,$WT,$MT,$ST,$JT,
			    $SNe,$SPe,$WE,$ME,$SNet,$SPet,$WEt,$MEt,
			    $SNi,$SPi,$WI,$MI,$SNit,$SPit,$WIt,$MIt,
			    $SNn,$SPn,$SNnt,$SPnt);
	}
    }

    ############################################################
    # print output
    print "---------------------------------------------------------------\n";
    unless ($pred_label =~/none/ && $ann_label =~/none/ ){
	print "ID\t".$pred_label."\t".$ann_label."\n";
    }
    print $header."\n";
    print $string."\n";
    print "---------------------------------------------------------------\n";
}




############################################################
# method to calculate the True Negativs
# This method takes all the transcripts together, predictions
# and annotations, cluster the exons and use the clusters as
# projections. With these projections we estimate then
# how much is sequence is not covered by both predictions and annotations:
#
#    ########-----------##########
#     #########---------#######
#              |-------| 
#                 TN
#
# this method only works with transcripts that are on the same
# coordinate system

sub calculate_TN{
    my ($trans) = @_;
    my @clusters = get_Transcript_Projections($trans);
    my @fw_clusters = sort { $a->start <=> $b->start } egrep { $_->strand == 1 } @clusters;
    my @rv_clusters = sort { $a->start <=> $b->start } egrep { $_->strand == -1 } @clusters;
    
    my $fw_length = $fw_clusters[-1]->end - $fw_clusters[0]->start + 1;
    my $rv_length = $rv_clusters[-1]->end - $rv_clusters[0]->start + 1;
 
    ############################################################
    # start with TN = the total length covered between the first base
    # and the last, and then remove from this regions those
    # that correspond to either annotation, prediction or both.
    my $fw_TN = $fw_length;
    foreach my $c (@fw_clusters){
	$fw_TN -= ($c->end - $c->start + 1);
    }
    my $rv_TN = $rv_length;
    foreach my $c (@rv_clusters){
	$rv_TN -= ($c->end - $c->start + 1);
    }
    
    if ( $rv_TN < 0 || $fw_TN < 0 ){
	print "ERROR: we are getting a negative TN: TN(forward) = $fw_TN, TN(reverse) = $rv_TN\n";
    }
    return ($fw_TN,$rv_TN);
}


############################################################

sub print_TranscriptPair{
    my $pair = shift;
    my $pred_tran = $pair->[0];
    my $ann_tran  = $pair->[1];
    my $CCn       = $pair->[2];
    my $TPn       = $pair->[3];
    my $FPn       = $pair->[4];
    my $FNn       = $pair->[5];
    ClusterMerge::TranscriptUtils->_print_SimpleTranscript($pred_tran);
    ClusterMerge::TranscriptUtils->_print_SimpleTranscript($ann_tran);
    print "CCn = $CCn\n";
}
