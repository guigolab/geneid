#!/usr/local/bin/perl  -w

############################################################
#
# script to compare genes
#
# written by Eduardo Eyras (eae@imim.es)
#
############################################################

use strict;
use Getopt::Long;
#use ClusterMerge::ClusterMerge;
use ClusterMerge::Transcript;
use ClusterMerge::TranscriptComparator;
use ClusterMerge::TranscriptUtils;
use ClusterMerge::Exon;
use GeneComparison::GeneComparison;
use Utils::GFFTools;

my $gff_format = 1;
my $input;
my $output;

my $verbose = 0;
my %type;


############################################################
# deafult options

my $prediction_file;
my $annotation_file;
my $out_pred = "out_pred";
my $out_ann  = "out_ann";
my $h;

############################################################

&GetOptions( 
	     'p:s' => \$prediction_file,
	     'a:s' => \$annotation_file,
	     'outp:s' => \$out_pred,
	     'outa:s' => \$out_ann,
	     'h:s' => \$h,                          
	     );

if ( !($prediction_file && $annotation_file && $out_pred && $out_ann ) ||  $h ){
    print STDERR "Usage: $0 -p <predictions> -a <annotations> -outp <out pred> -outa <out ann> -h <help>\n";
    exit(0);
}

############################################################
# read predictions

open ( IN, "<$prediction_file") || die("could not open input file $prediction_file");
my @predictions;
my $pred_trans_tag;
my %pred_trans;
my %forward_pred_sets;
my %reverse_pred_sets;
my %pred_gene_sets;
while(<IN>){
    chomp;
    my @e = split;
    next unless ( $e[3] && $e[4] ); 
    my %trans_tag;
    my $exon = Utils::GFFTools->exon_from_gff( $_ , \%trans_tag );
    #print STDERR "trans_tag[ exon ] = $trans_tag{$exon}\n";
    if ( $exon && $exon->start && $exon->end ){
	push( @{$pred_trans{$trans_tag{$exon}}}, $exon );
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
    push( @predictions, $transcript );
    $type{$transcript} = "prediction";
    if ($strand == -1 ){
	push( @{$reverse_pred_sets{$seqname}}, $transcript );
    }
    else{
	push( @{$forward_pred_sets{$seqname}}, $transcript );
    }
}
close(IN);


############################################################
# read annotations

open ( IN, "<$annotation_file") || die("could not open input file $annotation_file");
my @annotations;
my $ann_trans_tag;
my %ann_trans;
my %forward_ann_sets;
my %reverse_ann_sets;
my %ann_gene_sets;
while(<IN>){
    chomp;
    my @e = split;
    next unless ( $e[3] && $e[4] ); 
    my %trans_tag;
    my $exon = Utils::GFFTools->exon_from_gff( $_ , \%trans_tag );
    #print STDERR "trans_tag[ exon ] = $trans_tag{$exon}\n";
    if ( $exon && $exon->start && $exon->end ){
	push( @{$ann_trans{$trans_tag{$exon}}}, $exon );
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
	$strand = $exon->strand;
    }
    push( @annotations, $transcript );
    $type{$transcript} = "annotation";
    if ($strand == -1 ){
	push( @{$reverse_ann_sets{$seqname}}, $transcript );
    }
    else{
	push( @{$forward_ann_sets{$seqname}}, $transcript );
    }
}
close(IN);


#print STDERR "transcripts created\n";
#foreach my $t ( @transcripts ){
#    print STDERR "transcript: $t\n";
#    ClusterMerge::TranscriptUtils->_print_Transcript($t);
#}


my %forward_pred_clusters;
my %reverse_pred_clusters;
my %forward_ann_clusters;
my %reverse_ann_clusters;


############################################################
# put the seqnames together
my %forward_seqnames;
foreach my $seqname ( (keys %forward_pred_sets, keys %forward_ann_sets ) ){
    $forward_seqnames{$seqname} = 1;
}
my %reverse_seqnames;
foreach my $seqname ( (keys %reverse_pred_sets, keys %reverse_ann_sets ) ){
    $reverse_seqnames{$seqname} = 1;
}


############################################################
# cluster transcripts

my %forward_clusters;
my %reverse_clusters;

my %novel_predictions;
my %missed_annotations;
my %mixed_clusters;


############################################################
# clustering in the forward strand
#
# the clustering of transcripts is done by genomic-overlap
# i.e., a transcript is considered as a range in the genomic sequence
#
print "Forward clustering\n" if $verbose;

foreach my $seqname ( keys %forward_seqnames ){
    print "Looking at target: $seqname\n" if $verbose;
    if ( $forward_pred_sets{$seqname} && $forward_ann_sets{$seqname} ){
	push( @{$forward_clusters{$seqname}}, GeneComparison::GeneComparison->cluster_Transcripts([@{$forward_pred_sets{$seqname}}, @{$forward_ann_sets{$seqname}} ] ) );
	
	print scalar(  @{$forward_clusters{$seqname}} )." clusters found\n" if $verbose;
	my $count = 1;
	foreach my $c ( @{$forward_clusters{$seqname}} ){
	    my @predictions = grep { $type{$_} eq "prediction" } @{$c->get_Transcripts};
	    my @annotations = grep { $type{$_} eq "annotation" } @{$c->get_Transcripts};
	    
	    if ( !@annotations && @predictions){
		push( @{$novel_predictions{$seqname}}, @predictions );
		print scalar( @predictions )." novel predictions in target $seqname\n" if $verbose;
	    }
	    elsif ( @annotations && !@predictions){
		push( @{$missed_annotations{$seqname}}, @annotations );
		print scalar( @annotations )." missed annotations in target $seqname\n" if $verbose;
	    }
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
    elsif( !$forward_pred_sets{$seqname} && $forward_ann_sets{$seqname} ){
	push ( @{$missed_annotations{$seqname}}, @{$forward_ann_sets{$seqname}} );
	print scalar( @{$forward_ann_sets{$seqname}} )." missed annotations in target $seqname\n" if $verbose;
    }
    elsif( $forward_pred_sets{$seqname} && !$reverse_ann_sets{$seqname} ){
	push ( @{$novel_predictions{$seqname}} , @{$forward_pred_sets{$seqname}} );
	print scalar( @{$forward_pred_sets{$seqname}} )." novel predictions in target $seqname\n" if $verbose;
    }
}

############################################################
# clustering in the reverse strand
#
# the clustering of transcripts is done according to genomic overlap
#
print "*********Reverse clustering\n" if $verbose;
foreach my $seqname ( keys %reverse_seqnames ){
    print "Looking at target: $seqname\n" if $verbose;
    if ( $reverse_pred_sets{$seqname} && $reverse_ann_sets{$seqname} ){
	push( @{$reverse_clusters{$seqname}}, GeneComparison::GeneComparison->cluster_Transcripts([@{$reverse_pred_sets{$seqname}}, @{$reverse_ann_sets{$seqname}} ] ) );
	print scalar(  @{$reverse_clusters{$seqname}} )." clusters found\n" if $verbose;
	my $count = 1;
	foreach my $c ( @{$reverse_clusters{$seqname}} ){
	    my @predictions = grep { $type{$_} eq "prediction" } @{$c->get_Transcripts};
	    my @annotations = grep { $type{$_} eq "annotation" } @{$c->get_Transcripts};
	    
	    if ( !@annotations && @predictions){
		push( @{$novel_predictions{$seqname}}, @predictions );
	    	print scalar( @predictions )." novel predictions in target $seqname\n" if $verbose;
	    }
	    elsif ( @annotations && !@predictions){
		push( @{$missed_annotations{$seqname}}, @annotations );
	   	print scalar( @annotations )." missed annotations in target $seqname\n" if $verbose;
	    }
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
    elsif( !$reverse_pred_sets{$seqname} && $reverse_ann_sets{$seqname} ){
	push ( @{$missed_annotations{$seqname}}, @{$reverse_ann_sets{$seqname}} );
	print scalar( @{$reverse_ann_sets{$seqname}} )." missed annotations in target $seqname\n" if $verbose;
    }
    elsif( $reverse_pred_sets{$seqname} && !$reverse_ann_sets{$seqname} ){
	push ( @{$novel_predictions{$seqname}} , @{$reverse_pred_sets{$seqname}} );
	print scalar( @{$reverse_pred_sets{$seqname}} )." novel predictions in target $seqname\n" if $verbose;
    }
}



############################################################
# algorithm:
#
# we make exons uniq on each set (specially important for ensembl)
# for each possible pair prediction-annotation
#   we extract the common exon-assemblies on each side
#   where common means 'overlapping'
#
# we eliminate the exon assemblies which are redundant
# (that have exons which are included in a bigger/same size exon assembly)
#
# For cases like ensembl, it is important to
# get rid of exon-assemblies which are not adding any extra informations, i.e.
# they are not redundant in the sense of the previous check, but all exons
# of this particular assembly overlap exons of a bigger assembly.

#open(OUTP,">$out_pred") or die ("cannot open $out_pred");
#open(OUTA,">$out_ann")  or die ("cannot open $out_ann");


my $comparator = ClusterMerge::TranscriptComparator->new(	
								-comparison_level         => 3,
								-exon_match               => 0,	
								-splice_mismatch          => 6,
								-intron_mismatch          => 10,
								-internal_splice_overlap  => 10,
								);

 TARGET:
    foreach my $seqname ( keys %mixed_clusters ){
	
      TRANSCRIPT_CLUSTER:
	foreach my $c ( @{$mixed_clusters{$seqname}} ){
	    
	    #$c->print_aligned_Cluster;
	    my @predictions = grep { $type{$_} eq "prediction" } @{$c->get_Transcripts};
	    my @annotations = grep { $type{$_} eq "annotation" } @{$c->get_Transcripts};

	    ############################################################
	    # if there is overlap within predictions or annotations
	    # make exons unique:
	    if ( scalar(@predictions) > scalar(GeneComparison::GeneComparison->cluster_Transcripts(\@predictions) ) ){
		#print "making shared exons in predictions unique\n";
		@predictions = @{make_shared_exons_unique(\@predictions)};
	    }
	    if ( scalar(@annotations) > scalar(GeneComparison::GeneComparison->cluster_Transcripts(\@annotations) ) ){
		#print "making shared exons in annotations unique\n";
		@annotations = @{make_shared_exons_unique(\@annotations)};
	    }

	    my %exon_number;
	    my $exon_count = 1;
	    foreach my $t (@predictions){
		foreach my $exon ( @{$t->get_all_Exons} ){
		    $exon_number{$exon} = $exon_count;
		    $exon_count++;
		}
	    }
	    $exon_count = 1;
	    foreach my $t (@annotations){
		foreach my $exon ( @{$t->get_all_Exons} ){
		    $exon_number{$exon} = $exon_count;
		    $exon_count++;
		}
	    }



	    my @pred_assemblies;
	    my @ann_assemblies;
	    my %assembly2transcript;
	  PREDICTION:
	    foreach my $pred (@predictions){
	      ANNOTATION:
		foreach my $ann (@annotations){
		    push( @ann_assemblies, $comparator->find_intersecting_exon_assemblies( $ann, $pred ) );
		    foreach my $ass ( @ann_assemblies){
			$assembly2transcript{$ass} = $ann;
		    }
		    push( @pred_assemblies,$comparator->find_intersecting_exon_assemblies( $pred, $ann ) );
		    foreach my $ass ( @pred_assemblies){
			$assembly2transcript{$ass} = $pred;
		    }
		}
	    }
	    ############################################################
	    # go to the next cluster if there are no common assemblies
	    next TRANSCRIPT_CLUSTER unless (@pred_assemblies && @ann_assemblies);
	    
	    if ($verbose){
		print "Total prediction assemblies:\n";
		foreach my $ass ( @pred_assemblies ){
		    ClusterMerge::ExonUtils->print_Exons($ass);
		  }
		print "Total annotation assemblies:\n";
		foreach my $ass ( @ann_assemblies ){
		    ClusterMerge::ExonUtils->print_Exons($ass);
		  }
	    }
	    
	    ############################################################
	    # eliminate redundand assemblies
	    # by looking at lists of exon objects which are already included in another list
	    my @accepted_ann_assemblies  = ClusterMerge::ExonUtils->eliminate_redundant_lists(\@ann_assemblies);
	    my @accepted_pred_assemblies = ClusterMerge::ExonUtils->eliminate_redundant_lists(\@pred_assemblies);
	    
	    if ($verbose){
		print "------------------------------------------------------------\n";
		print "after pruning redundant assemblies\n";
		print "Total prediction assemblies:\n";
		foreach my $ass ( @accepted_pred_assemblies ){
		    ClusterMerge::ExonUtils->print_Exons($ass);
		  }
		print "Total annotation assemblies:\n";
		foreach my $ass ( @accepted_ann_assemblies ){
		    ClusterMerge::ExonUtils->print_Exons($ass);
		  }
	    }
	    
	    ############################################################
	    # in order to avoid overlapping quasi-identical
	    # assemblies (like for instance in ensembl):
	    #
	    #E1     ####---###---###---###
	    #E2      ###---###---####-------###
	    #
	    #S1 ##---###---###---######
	    #
	    #E1 and E2 would produce assemblies
	    #       ####---###---###
	    #        ###---###---####
	    # which possibly do not add any extra information,
	    # we should perhaps choose one among these.
	    
	    ############################################################
	    # convert the assemblies into transcripts:
	    my @ann_ass_trans  = convert_assemblies_into_transcripts( \@accepted_ann_assemblies,\%assembly2transcript );
	    my @ann_pred_trans = convert_assemblies_into_transcripts(\@accepted_pred_assemblies,\%assembly2transcript );
	    
	    ############################################################
	    # cluster transcripts 
	    my @ann_ass_clusters  = GeneComparison::GeneComparison->cluster_Transcripts(\@ann_ass_trans);
	    my @pred_ass_clusters = GeneComparison::GeneComparison->cluster_Transcripts(\@ann_pred_trans);
	    
	    # now the final assemblies are ClusterMerge::Transcript objects
	    
	    ############################################################
	    # eliminate the redundant exon-assemblies 
	    # considering a fuzzy overlap between the assemblies
	    my @final_ann_ass;
	    my @final_pred_ass;
	    foreach my $c ( @ann_ass_clusters ){
		my @trans = @{$c->get_Transcripts};
		if ( scalar(@trans) == 1 ){
		    push (@final_ann_ass, $trans[0]);
		    next;
		}
		my @non_redundant = ClusterMerge::TranscriptUtils->eliminate_redundant_transcripts(\@trans);
		push ( @final_ann_ass, @non_redundant);
	    }
	    foreach my $c ( @pred_ass_clusters ){
		my @trans = @{$c->get_Transcripts};
		if ( scalar(@trans) == 1 ){
		    push (@final_pred_ass, $trans[0]);
		    next;
		}
		my @non_redundant = ClusterMerge::TranscriptUtils->eliminate_redundant_transcripts(\@trans);
		push ( @final_pred_ass, @non_redundant);
	    }
	    
	    if ($verbose){
		print "------------------------------------------------------------\n";
		print "after merging overlapping assemblies\n";
		print "Total prediction assemblies:\n";
		foreach my $ass ( @final_pred_ass ){
		    ClusterMerge::ExonUtils->print_Exons($ass->get_all_Exons);
		  }
		print "Total annotation assemblies:\n";
		foreach my $ass ( @final_ann_ass ){
		    ClusterMerge::ExonUtils->print_Exons($ass->get_all_Exons);
		  }
	    }


	    # If we want to print only overlapping CDSs:

	    ############################################################
	    # from each overlapping pair of assemblies take
	    # the intersecting CDS only:

	    my @final_preds =  ClusterMerge::TranscriptUtils->sort_transcripts(\@final_pred_ass);
	    my @final_anns  =  ClusterMerge::TranscriptUtils->sort_transcripts(\@final_ann_ass);
	    
	    my $start = 0;
	  PRED:
	    foreach my $pred ( @final_preds ){

	      ANN:
		for(my $i=$start; $i<scalar(@final_anns); $i++ ){
		    
		    ############################################################
		    # if they overlap, at this stage they are necessarily related:
		    my $ann_start =  ClusterMerge::TranscriptUtils->transcript_low( $final_anns[$i] );
		    my $ann_end   =  ClusterMerge::TranscriptUtils->transcript_low( $final_anns[$i] );
		    my $pred_start=  ClusterMerge::TranscriptUtils->transcript_low( $pred );
		    my $pred_end  =  ClusterMerge::TranscriptUtils->transcript_low( $pred );
		    
		    if ( !( $ann_end < $pred_start || $ann_start > $pred_end ) ){
			
			# they overlap, hence we should generate a consensus between the two:
			my $consensus = $comparator->get_intersecting_CDS($pred,$final_anns[$i]);
			
			my $label = '';
			if ( scalar(@{$pred->get_all_Exons})           == scalar(@{$assembly2transcript{$pred}->get_all_Exons}) 
			     ||  
			     scalar(@{$final_anns[$i]->get_all_Exons}) == scalar(@{$assembly2transcript{$final_anns[$i]}->get_all_Exons}) 
			     ){
			    $label = "\tcomplete";
			}
			
			print goldenpath_string($consensus,$label);
			
			# jump to the next PRED and compare directly with the next ANN
			$start = $i+1;
			next PRED;
		    }
		}
	    }
	    
	    
	    # if we wanted to print separate exon-assemblies on each set:
	    
	    #foreach my $ass ( @final_pred_ass ){
	#	my $label ='';
	#	my @exons = @{$ass->get_all_Exons};
	#	foreach my $e ( @exons ){
	#	    $label .= $exon_number{$e}.":";
	#	}
	#	if ( scalar(@exons) == scalar(@{$assembly2transcript{$ass}->get_all_Exons}) ){
	#	    $label = "\tcomplete";
	#	}
	#	else{
	#	    $label = "\tincomp";
	#	}
	#	print OUTP goldenpath_string($ass,$label);
	#    }
	#    foreach my $ass ( @final_ann_ass ){
	#	my $label ='';
	#	my @exons = @{$ass->get_all_Exons};
	#	foreach my $e ( @exons ){
	#	    $label .= $exon_number{$e}.":";
	#	}
	#	if ( scalar(@exons) == scalar(@{$assembly2transcript{$ass}->get_all_Exons}) ){
	#	    $label = "\tcomplete";
	#	}
	#	else{
	#	    $label = "\tincomp";
	#	}
	#	print OUTA goldenpath_string($ass,$label);
	#    }
	    

	    
	} # end of TRANSCRIPT_CLUSTER
    } # end of TARGET

#close(OUTP);
#close(OUTA);

############################################################
# add whole novel predictions
#
# these are predictions that are not clustered with any annotation
#
#print "WHOLE NOVEL PREDICTIONS\n" if $verbose;
#foreach my $seqname ( keys %novel_predictions ){
#    foreach my $t ( @{$novel_predictions{$seqname}} ){
#	#ClusterMerge::TranscriptUtils->_print_SimpleTranscript($t);
#	goldenpath_string($t,"complete");
#    }
#}


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
  
  my $tag_str = $transcript->type || $transcript->dbID;
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

    my $string = "$id\t$seqname\t$strand\t$low\t$high\t$low\t$high\t$exon_number\t$start_string\t$end_string";
    if ($label){
	$string .= "\t$label";
    }
    $string .= "\n";
}



############################################################


sub make_shared_exons_unique{
  my ( $transcripts ) = @_;
  
  my @unique_Exons; 
  
  # keep track of all unique exons found so far to avoid making duplicates
  # need to be very careful about translation->start_Exon and translation->end_Exon
  
  foreach my $tran (@$transcripts) {
      my @newexons;
      foreach my $exon (@{$tran->get_all_Exons}) {
	  my $found;
	  #always empty
	UNIQUE_EXON:
	  foreach my $uni (@unique_Exons) {
	      if ($uni->start     == $exon->start  &&
		  $uni->end       == $exon->end    &&
		  $uni->strand    == $exon->strand &&
		  $uni->phase     == $exon->phase  &&
		  $uni->end_phase == $exon->end_phase
		  ) {
		  $found = $uni;
		  last UNIQUE_EXON;
	      }
	  }
	  if (defined($found)) {
	      push(@newexons,$found);
	  } 
	  else {
	      push(@newexons,$exon);
	      push(@unique_Exons, $exon);
	  }
      }          
      $tran->flush_Exons;
      foreach my $exon (@newexons) {
	  $tran->add_Exon($exon);
      }
  }
  return $transcripts;
}

############################################################


sub convert_assemblies_into_transcripts{
    my ($assemblies, $ass2tra ) = @_;
    
    my @transcripts;
    foreach my $ass ( @$assemblies ){
	my $trans = ClusterMerge::Transcript->new();
	
	$ass2tra->{$trans} = $ass2tra->{$ass};
	$trans->dbID( $ass2tra->{$ass}->dbID );
	foreach my $e ( @$ass ){
	    $trans->add_Exon($e);
	}
	push( @transcripts, $trans);
    }
    
    return @transcripts;
}
    
