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

my $prediction;
my $annotation;
my $h;

############################################################

&GetOptions( 
	     'p:s' => \$prediction,
	     'a:s' => \$annotation,
	     'h:s' => \$h,                          
	     );

if ( !($prediction || $annotation) ||  $h ){
    print STDERR "Usage: $0 -p <predictions> -a <annotations> -h <help>\n";
    exit(0);
}


print "predictions: $prediction\tannotation: $annotation\n";

############################################################
# read predictions

open ( IN, "<$prediction") || die("could not open input file $prediction");
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

open ( IN, "<$annotation") || die("could not open input file $annotation");
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
my %all_seqnames;

my %forward_seqnames;
foreach my $seqname ( (keys %forward_pred_sets, keys %forward_ann_sets ) ){
    $forward_seqnames{$seqname} = 1;
    $all_seqnames{$seqname} = 1;
}
my %reverse_seqnames;
foreach my $seqname ( (keys %reverse_pred_sets, keys %reverse_ann_sets ) ){
    $reverse_seqnames{$seqname} = 1;
    $all_seqnames{$seqname} = 1;
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

my ($total_pred_solo, $total_ann_solo, $total_pred_overlaps, $total_ann_overlaps) = (0,0,0,0);
print "############################################################\n";
print "Forward clustering\n";
foreach my $seqname ( keys %forward_seqnames ){
    print "Looking at target: $seqname\n" if $verbose;

    my $pred_overlaps= 0;
    my $ann_overlaps = 0;
    my $pred_solo    = 0;
    my $ann_solo     = 0;
    if ( $forward_pred_sets{$seqname} && $forward_ann_sets{$seqname} ){
	push( @{$forward_clusters{$seqname}}, GeneComparison::GeneComparison->cluster_Transcripts([@{$forward_pred_sets{$seqname}}, @{$forward_ann_sets{$seqname}} ] ) );
	
	print scalar(  @{$forward_clusters{$seqname}} )." clusters found\n" if $verbose;
	my $count = 1;
	foreach my $c ( @{$forward_clusters{$seqname}} ){
	    my @predictions = grep { $type{$_} eq "prediction" } @{$c->get_Transcripts};
	    my @annotations = grep { $type{$_} eq "annotation" } @{$c->get_Transcripts};
	    
	    if ( !@annotations && @predictions){
		#push( @{$novel_predictions{$seqname}}, @predictions );
		#print scalar( @predictions )." novel predictions in target $seqname\n";
		$pred_solo += scalar(@predictions);
	    }
	    elsif ( @annotations && !@predictions){
		#push( @{$missed_annotations{$seqname}}, @annotations );
		#print scalar( @annotations )." missed annotations in target $seqname\n";
		$ann_solo += scalar(@annotations);
	    }
	    elsif(  @annotations && @predictions ){
		push( @{$mixed_clusters{$seqname}}, $c ); 
		#print "target:$seqname\t"."annotations:".scalar( @annotations )."\tpredictions:".scalar( @predictions )."\n";
		$pred_overlaps += scalar(@predictions);
		$ann_overlaps  += scalar(@annotations);
		#if ($verbose){
		#    print "Cluster $count:\n";
		#    $c->print_Cluster(\%type);
		#    $count++;
		#}
	    }
	}
	
    }
    elsif( !$forward_pred_sets{$seqname} && $forward_ann_sets{$seqname} ){
	#push ( @{$missed_annotations{$seqname}}, @{$forward_ann_sets{$seqname}} );
	$ann_solo += scalar(@{$forward_ann_sets{$seqname}});
	#print scalar( @{$forward_ann_sets{$seqname}} )." missed annotations in target $seqname\n" if $verbose;
    }
    elsif( $forward_pred_sets{$seqname} && !$reverse_ann_sets{$seqname} ){
	#push ( @{$novel_predictions{$seqname}} , @{$forward_pred_sets{$seqname}} );
	$pred_solo += scalar( @{$forward_pred_sets{$seqname}});
	#print scalar( @{$forward_pred_sets{$seqname}} )." novel predictions in target $seqname\n" if $verbose;
    }
    print "target:$seqname\t\tann_solo:$ann_solo\tpred_solo:$pred_solo\tann_overl:$ann_overlaps\tpred_overl:$pred_overlaps\n";
    $total_pred_solo     += $pred_solo;
    $total_ann_solo      +=  $ann_solo;
    $total_pred_overlaps +=  $pred_overlaps;
    $total_ann_overlaps  +=  $ann_overlaps;
}

print "------------------------------------------------------------\n";
print "TOTAL-FORWARD\t\ttotal_ann_solo:$total_ann_solo\ttotal_pred_solo:$total_pred_solo\ttotal_ann_overl:$total_ann_overlaps\ttotal_pred_overl:$total_pred_overlaps\n";
print "------------------------------------------------------------\n";


($total_pred_solo, $total_ann_solo, $total_pred_overlaps, $total_ann_overlaps) = (0,0,0,0);

############################################################
# clustering in the reverse strand
#
# the clustering of transcripts is done according to genomic overlap
#
print "############################################################\n";
print "Reverse clustering\n";
foreach my $seqname ( keys %reverse_seqnames ){
    print "Looking at target: $seqname\n" if $verbose;
    my $pred_overlaps= 0;
    my $ann_overlaps = 0;
    my $pred_solo    = 0;
    my $ann_solo     = 0;
    if ( $reverse_pred_sets{$seqname} && $reverse_ann_sets{$seqname} ){
	push( @{$reverse_clusters{$seqname}}, GeneComparison::GeneComparison->cluster_Transcripts([@{$reverse_pred_sets{$seqname}}, @{$reverse_ann_sets{$seqname}} ] ) );
	print scalar(  @{$reverse_clusters{$seqname}} )." clusters found\n" if $verbose;
	my $count = 1;
	foreach my $c ( @{$reverse_clusters{$seqname}} ){
	    my @predictions = grep { $type{$_} eq "prediction" } @{$c->get_Transcripts};
	    my @annotations = grep { $type{$_} eq "annotation" } @{$c->get_Transcripts};
	    
	    if ( !@annotations && @predictions){
		#push( @{$novel_predictions{$seqname}}, @predictions );
	    	#print scalar( @predictions )." novel predictions in target $seqname\n" if $verbose;
		$pred_solo += scalar(@predictions);
	    }
	    elsif ( @annotations && !@predictions){
		#push( @{$missed_annotations{$seqname}}, @annotations );
	   	#print scalar( @annotations )." missed annotations in target $seqname\n" if $verbose;
		$ann_solo += scalar(@annotations);

	    }
	    elsif(  @annotations && @predictions ){
		$pred_overlaps += scalar(@predictions);
		$ann_overlaps  += scalar(@annotations);
		#push( @{$mixed_clusters{$seqname}}, $c ); 
		#print "target:$seqname\t"."annotations:".scalar( @annotations )."\tpredictions:".scalar( @predictions )."\n";
		#if ($verbose){
		#    print "Cluster $count:\n";
		#    $c->print_Cluster(\%type);
		#    $count++;
		#}
	    }
	}
    }
    elsif( !$reverse_pred_sets{$seqname} && $reverse_ann_sets{$seqname} ){
	#push ( @{$missed_annotations{$seqname}}, @{$reverse_ann_sets{$seqname}} );
	$ann_solo += scalar( @{$reverse_ann_sets{$seqname}});
	#print scalar( @{$reverse_ann_sets{$seqname}} )." missed annotations in target $seqname\n";
    }
    elsif( $reverse_pred_sets{$seqname} && !$reverse_ann_sets{$seqname} ){
	$pred_solo += scalar(@{$reverse_pred_sets{$seqname}});
	#push ( @{$novel_predictions{$seqname}} , @{$reverse_pred_sets{$seqname}} );
	#print scalar( @{$reverse_pred_sets{$seqname}} )." novel predictions in target $seqname\n";
    }
    print "target:$seqname\t\tann_solo:$ann_solo\tpred_solo:$pred_solo\tann_overl:$ann_overlaps\tpred_overl:$pred_overlaps\n";
    $total_pred_solo     +=  $pred_solo;
    $total_ann_solo      +=  $ann_solo;
    $total_pred_overlaps +=  $pred_overlaps;
    $total_ann_overlaps  +=  $ann_overlaps;
}

print "------------------------------------------------------------\n";
print "TOTAL-REVERSE\t\ttotal_ann_solo:$total_ann_solo\ttotal_pred_solo:$total_pred_solo\ttotal_ann_overl:$total_ann_overlaps\ttotal_pred_overl:$total_pred_overlaps\n";
print "------------------------------------------------------------\n";



############################################################
# %novel_predictions contains all novel predictions (no exon overlap at all)
# %missed_predictions contains all missed annotations (no exon overlap at all)
# %mixed_clusters contains all clusters with predictions and annotations

my %novel_assemblies;

($total_pred_solo, $total_ann_solo, $total_pred_overlaps, $total_ann_overlaps) = (0,0,0,0);

############################################################
# strandless clustering:

print "############################################################\n";
print "Strandless Clustering\n";
my %all_clusters;

foreach my $seqname ( keys %all_seqnames ){
    my $pred_overlaps= 0;
    my $ann_overlaps = 0;
    my $pred_solo    = 0;
    my $ann_solo     = 0;
    my @all_transcripts;
    if ( $reverse_pred_sets{$seqname} ){
	push( @all_transcripts, @{$reverse_pred_sets{$seqname}} );
    }
    if ( $forward_pred_sets{$seqname} ){
	push( @all_transcripts, @{$forward_pred_sets{$seqname}} );
    }
    if ( $reverse_ann_sets{$seqname} ){
	push( @all_transcripts, @{$reverse_ann_sets{$seqname}} );
    }
    if ( $forward_ann_sets{$seqname} ){
	push( @all_transcripts, @{$forward_ann_sets{$seqname}} );
    }
    push( @{$all_clusters{$seqname}}, GeneComparison::GeneComparison->cluster_Transcripts(\@all_transcripts) );

    foreach my $c ( @{$all_clusters{$seqname}} ){
	my @prediction = grep { $type{$_} eq "prediction" } @{$c->get_Transcripts};
	my @annotation = grep { $type{$_} eq "annotation" } @{$c->get_Transcripts};
	
	if ( !@annotation && @prediction){
	    $pred_solo += scalar(@prediction);
	}
	elsif ( @annotation && !@prediction){
	    $ann_solo += scalar(@annotation);
	}
	elsif(  @annotation && @prediction ){
	    $pred_overlaps += scalar(@prediction);
	    $ann_overlaps  += scalar(@annotation);
	}
    }
    print "target:$seqname\t\tann_solo:$ann_solo\tpred_solo:$pred_solo\tann_overl:$ann_overlaps\tpred_overl:$pred_overlaps\n";
    $total_pred_solo     +=  $pred_solo;
    $total_ann_solo      +=  $ann_solo;
    $total_pred_overlaps +=  $pred_overlaps;
    $total_ann_overlaps  +=  $ann_overlaps;
}

print "------------------------------------------------------------\n";
print "TOTAL-STRANDLESS\t\ttotal_ann_solo:$total_ann_solo\ttotal_pred_solo:$total_pred_solo\ttotal_ann_overl:$total_ann_overlaps\ttotal_pred_overl:$total_pred_overlaps\n";
print "------------------------------------------------------------\n";
							   
								   




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
    print "$id\t$seqname\t$strand\t$low\t$high\t$low\t$high\t$exon_number\t$start_string\t$end_string";
    if ($label){
	print "\t$label";
	}
    print "\n";
}
