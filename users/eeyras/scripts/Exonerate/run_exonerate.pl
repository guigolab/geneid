#!/usr/local/bin/perl -w

use strict;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Getopt::Long;
use Bio::SeqIO;

my $query;
my $target;
my $query_type  = 'dna';
my $target_type = 'dna';

my $min_coverage = 90;
my $min_perc_id  = 97;

my $verbose = 1;

&GetOptions( 
	     'target:s' => \$target,
	     'query:s'  => \$query,
	     'target_type:s' => \$target_type,
	     'query_type:s'  => \$query_type,
	     );

unless ( $target && $query ){
    print STDERR "Usage: $0 -target <target_file> -query <query_file>\n";
    exit(0);
}

############################################################
# get the length of the query sequences
############################################################
my %length;
my %coverage;
my %perc_id;
my %query_id;
my %target_id;
my $seqio= Bio::SeqIO->new(
			   -format => "Fasta",
			   -file   => "$query",
			  );

while( my $seq = $seqio->next_seq() ){
    if (defined($seq) && !($seq->seq eq '') && !($seq->display_id eq '') ){
	$length{$seq->display_id} = $seq->length;
    }
}


############################################################
# run exonerate
############################################################


my $options = " --softmasktarget  --score 500 --fsmmemory 800 ".
    " --saturatethreshold 100 --dnawordlen 14 --exhaustive FALSE --percent 90 ".
    " --model est2genome --ryo \"RESULT: %S %pi %g %V\\n\" ";

# around 60 cdnas with parameters:
# --softmasktarget  --score 500 --fsmmemory 800
# --saturatethreshold 100 --dnawordlen 14 --exhaustive FALSE
# --model est2genome --ryo
#
#real	6m5.496s  (real time)
#user	2m24.980s (user CPU)
#sys	0m9.260s  (system CPU)

# with exhaustive and percent 90 seems to get out of memory


my $command = "exonerate ".
    $options." --querytype $query_type --targettype $target_type --query $query --target $target/* ";
  
$command .= " | ";
  
print STDOUT "running exonerate: $command\n";
  
open( EXO, $command ) || die("Error running exonerate $!");
  
# system calls return 0 (true in Unix) if they succeed

############################################################
# store each alignment as a transcript with supporting features
my @transcripts;

############################################################
# parse results - avoid writing to disk the output
while (<EXO>){
    
    print STDOUT $_ if $verbose;
    
    ############################################################
    # the output is of the format:
    #
    # --ryo "RESULT: %S %pi %V\n"
    # 
    # It shows the alignments in "sugar" + percent_id + gene orientation + "vulgar blocks" format. 
    #
    # Sugar contains 9 fields
    # ( <qy_id> <qy_start> <qy_len> <qy_strand> <tg_id> <tg_start> <tg_len> <tg_strand> <score> ), 
    # 
    # The vulgar (Verbose Useful Labelled Gapped Alignment Report) blocks are a series 
    # of <label, query_length, target_length> triplets. The label may be one of the following: 
    #
    # M    Match 
    # G    Gap 
    # C    Codon gap 
    # N    Non-equivalenced region 
    # 5/F  5' splice site 
    # 3/T  3' splice site 
    # I    Intron 
    #
    # example:
    # RESULT: AW793782.1 25 250 + 10_NT_035057.1-141075 104447 126318 + 652 82.88 M 30 30 5 0 2 I 0 21645 3 0 2 M 35 35 G 1 0 M 16 16 G 1 0 M 134 134 G 1 0 M 7 7
    # 
    # This gives rise to:
    # M 30 30  ---> exon
    # 5 0 2 
    # I 0 21645 --> intron
    # 3 0 2 
    # M 35 35    \
    # G 1 0      |
    # M 16 16    |
    # G 1 0      |-> exon
    # M 134 134  |
    # G 1 0      |
    # M 7 7     /
    # 
    chomp;
    next unless ( /RESULT/ );
    my ( $tag, $q_id, $q_start, $q_end, $q_strand, $t_id, $t_start, $t_end, $t_strand, $score, $perc_id, $gene_orientation, @blocks) = split;
    
    # the SUGAR 'start' coordinates are 1 less than the actual position on the sequence
    $q_start++;
    $t_start++;
    
    next unless ( $tag && $tag eq 'RESULT:' );
   
    ############################################################
    # initialize the feature
    my (%query, %target);
    my (%rev_query, %rev_target);    
    
    ( $query{score}  , $target{score} )  = ($score,$score);
    ( $query{percent}, $target{percent}) = ( $perc_id, $perc_id );
    ( $query{source} , $target{source} ) = ('exonerate','exonerate');
    ( $query{name}   , $target{name})    = ( $q_id, $t_id);
 
    if ( $q_strand eq '+' ){ $q_strand = 1; }
    else{ $q_strand = -1; }
    if ( $t_strand eq '+' ){ $t_strand = 1; }
    else{ $t_strand = -1; }
    
    ( $query{strand} , $target{strand})  = ( $q_strand, $t_strand );
    

    ############################################################
    # coordinates are corrected later according to strand
    $query{start}  = $q_start;
    $target{start} = $t_start;
        
    ############################################################
    # orientation is not used at the moment
    print STDOUT "gene_orientation: $gene_orientation\n" if $verbose;

    if ( $gene_orientation eq '+' ){   $gene_orientation = 1 }
    elsif( $gene_orientation eq '-' ){ $gene_orientation = -1 }
    else{ $gene_orientation = 0 }
    
    my $transcript = Bio::EnsEMBL::Transcript->new();      
    $query_id{ $transcript } = $q_id;
    $target_id{ $transcript } = $t_id;
    $perc_id{ $transcript } = $perc_id;
    my $exon = Bio::EnsEMBL::Exon->new();
    my @features;
    my $in_exon = 0;
    my $target_gap_length = 0;
    
    ############################################################
    # exons are delimited by M - I
    # supporting features are delimited by M - G  and  M - I
  TRIAD:
    for( my $i=0; $i<=$#blocks; $i+=3){
      
	# do not look at splice sites now
	if ( $blocks[$i] eq '5' || $blocks[$i] eq 'F' 
	     || 
	     $blocks[$i] eq '3' || $blocks[$i] eq 'T'
	     ){
	    
	    # the re-set of the coordinates is done in the intron triad
	    next TRIAD;
	}
	
	############################################################
	# match
	if ( $blocks[$i] eq 'M' ){
	    if ( $in_exon == 0 ){
		$in_exon = 1;
		
		# start a new exon
		$exon = Bio::EnsEMBL::Exon->new();
		$exon->seqname( $t_id );
		$exon->start( $target{start} );
		$exon->strand( $gene_orientation );
		$exon->phase( 0 );
		$exon->end_phase( 0 );
	    }
	    # for every match we increase the end coordinate of the current exon
	    
	    if ( $exon->end ){
		#print STDOUT "shifting exon end to ". ($exon->end +  $blocks[$i+2] )."\n";
		$exon->end( $exon->end +  $blocks[$i+2] );
	    }
	    else{
		#print STDOUT "setting exon end = ". ($exon->start +  $blocks[$i+2] - 1 )."\n";
		$exon->end( $exon->start +  $blocks[$i+2] - 1 );
	    }
	    #start a new feature pair
	    $query{end}  = $query{start}  + $blocks[$i+1] - 1;
	    $target{end} = $target{start} + $blocks[$i+2] - 1;
	    
	    if ($verbose){
		print STDOUT "FEATURE:\n";
		#print STDOUT "query length: ".$length{$q_id}."\n";
		print STDOUT "query : $query{start} -  $query{end}\n";
		print STDOUT "target: $target{start} -  $target{end}\n";
	    }
	    
	    ############################################################
	    # if the query has been inverted we count from the end and
	    # invert start and end to inforce start < end
	    
	    %rev_query  = %query;
	    %rev_target = %target;
	    if ( $q_strand == -1 ){
		$rev_query{end}   = $length{$q_id} - $query{start} + 1;
		$rev_query{start} = $length{$q_id} - $query{end} + 1;
		
		#print STDOUT "rev_query{end}   = ".($length{$q_id})." - ".$query{start}." + 1\n" if $verbose;
		#print STDOUT "rev_query{start} = ".($length{$q_id})." - ".$query{start}." + 1\n" if $verbose;
		#my $feature_pair = create_FeaturePair(\%rev_target, \%rev_query);	
		
		#print STDOUT "adding feature: ".$feature_pair->gffstring."\n" if $verbose;
		#push( @features, $feature_pair);
	    }
	    else{	      	    
		#my $feature_pair = create_FeaturePair(\%target, \%query);
		#print STDOUT "adding feature: ".$feature_pair->gffstring."\n" if $verbose;
		#push( @features, $feature_pair);
	    }
	    
	    
	    # re-set the start:
	    $query{start}  = $query{end}  + 1;
	    $target{start} = $target{end} + 1;
	    #if ($verbose){
	    #  print STDOUT "Re-SET POINTERS:\n";
	    #  print STDOUT "query : $query{start} -  $query{end}\n";
	    #  print STDOUT "target: $target{start} -  $target{end}\n";
	    #}
	    
      }
      ############################################################
      # gap 
	if ( $blocks[$i] eq 'G' ){
	    print STDOUT "in_exon: $in_exon GAP: $blocks[$i+1] $blocks[$i+2]\n" if $verbose;
	    if ( $in_exon ){
		# keep the same exon
		
		# if the gap is in the query, we move the target
		if ( $blocks[$i+2] ){
		    $target{start} += $blocks[$i+2];
		    # also move the exon:
		    print STDOUT "moving exon end from ". $exon->end. " to ". ( $exon->end + $blocks[$i+2] )."\n" if $verbose;
		    $exon->end( $exon->end + $blocks[$i+2]); 
		}
		# if the gap is in the target, we move the query
		if ( $blocks[$i+1] ){
		    $query{start} += $blocks[$i+1];
		    $target_gap_length += $blocks[$i+1];
		}
		#if ($verbose){
		#  print STDOUT "GAP:\n";
		#  print STDOUT "query : $query{start} -  $query{end}\n";
		#	print STDOUT "target: $target{start} -  $target{end}\n";
		#  }
	    }
	}
	############################################################
	# intron
	if( $blocks[$i] eq 'I' ){
	    if ( $in_exon ){
		# emit the current exon
		if ($verbose){
		    print STDOUT "EXON: ".$exon->start."-".$exon->end."\n";
		}
		
		$transcript->add_Exon($exon);
		
		# add the supporting features
		#my $supp_feature;
		#eval{
		#    $supp_feature = Bio::EnsEMBL::DnaDnaAlignFeature->new( -features => \@features);
		#};
		#print STDOUT "intron: adding evidence : ".$supp_feature->gffstring."\n" if $verbose;
		#$exon->add_supporting_features( $supp_feature );
		
		$in_exon = 0;
		
		@features = ();
		
		# reset the start in the target only
		my $intron_length = $blocks[$i+2];
		if ( $i>=3 && 
		     ( $blocks[$i-3] eq '3' || $blocks[$i-3] eq 'T' )
		     || 
		     ( $blocks[$i-3] eq '5' || $blocks[$i-3] eq 'F' ) 
		     ){
		    $intron_length += 2;
		}
		if ( ( $blocks[$i+3] eq '5' || $blocks[$i+3] eq 'F' )
		     || 
		     ( $blocks[$i+3] eq '3' || $blocks[$i+3] eq 'T' )
		     ){
		    $intron_length += 2;
		}
		
		#print STDOUT "intron length = $intron_length\n" if $verbose;
		$target{start} += $intron_length;
		print STDOUT "new target starts at $target{start}\n" if $verbose;
	    }
	}
    } # end of TRIAD
    
    ############################################################
    # emit the last exon and the last set of supporting features
    # and add the supporting features
    #print STDOUT "created features:\n";
    #foreach my $f (@features){
    #print STDOUT $f->gffstring."\n";
    #}
    
    if ($verbose){
	print STDOUT "EXON: ".$exon->start."-".$exon->end."\n";
    }


#if ( scalar(@features) ){
#	my $supp_feature;
#	eval{
#	    $supp_feature = Bio::EnsEMBL::DnaDnaAlignFeature->new( -features => \@features);
#	};
#	if ($@){
#	    print STDOUT $@."\n";
#	}
#	#print STDOUT "outside: adding evidence : ".$supp_feature->gffstring."\n" if $verbose;
#	$exon->add_supporting_features( $supp_feature );
    #   }
    #   else{
#	#print STDOUT "No more features to add" if $verbose;
#    }
    $transcript->add_Exon($exon);
    
    
    ############################################################
    # compute coverage
    # q_start reported by the sugar/cigar lines is one less 
    # as exonerate counts between lines
    my $aligned_length = $q_end - ($q_start + 1) + 1;
    
    my $coverage = sprintf "%.2f", 100 * ( $aligned_length - $target_gap_length ) / $length{$q_id};
    print STDOUT "coverage = ( $aligned_length - $target_gap_length ) / ".$length{$q_id}." =  $coverage\n" if $verbose;
    
    $coverage{$transcript} = $coverage;
    
    ############################################################
    # print GFF of this transcript
    printGFF($transcript);
    
    ############################################################
    # check splice sites
    my $chr_location = $target."/".$t_id.".fa";
    check_splice_sites( $transcript, $chr_location );


    ############################################################
    # lower bound to avoid unnecessary processing
    if ( $coverage > 60 && $perc_id > 60 ){
	print STDOUT "****************************** KEEP **********************************\n";
	push( @transcripts, $transcript );
    }
    
} # end of while loop

close(EXO) || die("couldn't close pipe ");  

############################################################
# filter alignments
my @alignments = filter_alignments(@transcripts);

# test #
print STDOUT "Results:\n";
foreach my $t ( @transcripts ){
    print_SimpleTranscript( $t);
}







############################################################
# create feature pairs
############################################################

sub create_FeaturePair {
    my ($feat1, $feat2) = @_;
    #create analysis object
    my $analysis_obj = new Bio::EnsEMBL::Analysis
                        (   -db              => $feat2->{db},
                            -db_version      => $feat2->{db_version},
                            -program         => $feat2->{program},
                            -program_version => $feat2->{p_version},
                            -gff_source      => $feat2->{source},
                            -gff_feature     => $feat2->{primary},
                            -logic_name      => $feat2->{logic_name} );
    
    #create and fill Bio::EnsEMBL::Seqfeature objects
    my $seqfeature1 = new Bio::EnsEMBL::SeqFeature
                        (   -seqname        => $feat1->{name},
                            -start          => $feat1->{start},
                            -end            => $feat1->{end},
                            -strand         => $feat1->{strand},
                            -score          => $feat1->{score},
                            -percent_id     => $feat1->{percent},
                            -p_value        => $feat1->{p},
                            -analysis       => $analysis_obj);
    
    my $seqfeature2 = new Bio::EnsEMBL::SeqFeature
                        (   -seqname        => $feat2->{name},
                            -start          => $feat2->{start},
                            -end            => $feat2->{end},
                            -strand         => $feat2->{strand},
                            -score          => $feat2->{score},
                            -percent_id     => $feat2->{percent},
                            -p_value        => $feat2->{p},
                            -analysis       => $analysis_obj);
    #create featurepair
    my $fp = Bio::EnsEMBL::FeaturePair->new  (  -feature1 => $seqfeature1,
                                                -feature2 => $seqfeature2 ) ;

    #print "Feature pair " . $fp->gffstring . "\n";

    return $fp;
}


############################################################
# method for printing the transcript result
############################################################

sub print_SimpleTranscript{
    my ($transcript,$chr_coord) = @_;
    my @exons = sort { $a->start <=> $b->start } @{$transcript->get_all_Exons};
    
    print STDOUT "query:".$query_id{$transcript}.
	" coverage:".$coverage{$transcript}.
	" perc_id:".$perc_id{$transcript}.
	" target:".$target_id{$transcript}." ";
    
    my $shift = 0;
    if ( $chr_coord ){
      $shift = $exons[0]->contig->chr_start - 1;
    }
    foreach my $exon ( @exons){
	print STDOUT ($exon->start + $shift)."-".( $exon->end + $shift )." ";
    }
    print STDOUT "\n";
}


############################################################
# method for GFF printing
############################################################

sub printGFF{
    my ($transcript,$chr_coord) = @_;
    my @exons = sort { $a->start <=> $b->start } @{$transcript->get_all_Exons};
    
    my $shift = 0;
    if ( $chr_coord ){
      $shift = $exons[0]->contig->chr_start - 1;
    }
    # exon coordinates may be shifted : print STDOUT ($exon->start + $shift)."-".( $exon->end + $shift )." ";
    
    foreach my $exon ( @exons ){
	my $strand = "+";
	if ( $exon->strand == -1 ){
	    $strand = "-";
	}
	
	print STDOUT 
	    $target_id{$transcript}."\t".
	    "cDNA"."\t".
	    "exon"."\t".
	    $exon->start."\t".
	    $exon->end."\t".
	    $coverage{$transcript}."\t".
	    $strand."\t".
	    "."."\t".
	    $query_id{$transcript}."\t".
	    "# perc_id:".$perc_id{$transcript}."\n";
	
    }	
}


############################################################
# best in genome
############################################################

sub filter_alignments{
  my (@transcripts) = @_;
  
  # results are Bio::EnsEMBL::Transcripts with exons and supp_features
  
  my @good_matches;

  my %matches;

  ############################################################
  # bin transcripts according to the identifier
  ############################################################
  foreach my $transcript (@transcripts ){
      push ( @{$matches{$query_id{$transcript}}}, $transcript );
  }
  
  my %matches_sorted_by_coverage;
  my %selected_matches;
  
 RNA:
  foreach my $rna_id ( keys( %matches ) ){
    
      @{$matches_sorted_by_coverage{$rna_id}} = 
	sort { my $result = ( $coverage{$b} <=> $coverage{$a} );
	       if ( $result){
		   return $result;
	       }
	       else{
		   my $result2 = ( scalar(@{$b->get_all_Exons}) <=> scalar(@{$a->get_all_Exons}) );
		   if ( $result2 ){
		       return $result2;
		   }
		   else{
		       return ( $perc_id{$b} <=> $perc_id{$a} );
		   }
	       }
	   }   @{$matches{$rna_id}} ;
      
      
      my $count = 0;
      my $is_spliced = 0;
      my $max_score;
      my $perc_id_of_best = 0;
      my $best_has_been_seen = 0;
      
      print STDOUT "####################\n";
      print STDOUT "Matches for $rna_id:\n";
      
    TRANSCRIPT:
      foreach my $transcript ( @{$matches_sorted_by_coverage{$rna_id}} ){
	  $count++;
	  unless ($max_score){
	      $max_score = $coverage{$transcript};
	  }
	  unless ( $perc_id_of_best ){
	      $perc_id_of_best = $perc_id{$transcript};
	  }
	 	  
	  my $score   = $coverage{$transcript};
	  my $perc_id = $perc_id{$transcript};
	  
	  my @exons  = sort { $a->start <=> $b->start } @{$transcript->get_all_Exons};
	  my $start  = $exons[0]->start;
	  my $end    = $exons[$#exons]->end;
	  my $strand = $exons[0]->strand;
	  my $seqname= $exons[0]->seqname;
	  $seqname   =~ s/\.\d+-\d+$//;
	  my $extent = $seqname.".".$start."-".$end;
	  
	  ############################################################
	  # put flag is the first one is spliced
	  if ( $count == 1 && is_spliced( $transcript ) ){
	      $is_spliced = 1;
	  }

	  my $label;
	  if ( $count == 1 ){
	      $label = 'best_match';
	  }
	  elsif ( $count > 1 
		  && $is_spliced 
		  && !is_spliced( $transcript )
		  ){
	      $label = 'potential_processed_pseudogene';
	  }
	  else{
	      $label = $count;
	  }
	 
	  my $accept;
	
	  ############################################################
	  # FILTERING:
	  my $best_in_genome  = 1;
	  my $reject_potential_pseudos = 1;

	  if ($best_in_genome){
	      if ( ( $score  == $max_score && 
		     $score >= $min_coverage && 
		     $perc_id >= $min_perc_id
		     )
		   ||
		   ( $score == $max_score &&
		     $score >= (1 + 5/100)*$min_coverage &&
		     $perc_id >= ( 1 - 3/100)*$min_coverage
		     )
		   ){
		  if ( $reject_potential_pseudos
		       && $count > 1 
		       && $is_spliced 
		       && !is_spliced( $transcript )
		       ){
		      $accept = 'NO';
		  }
		  else{
		      $accept = 'YES';
		      push( @good_matches, $transcript);
		  }
	      }
	      else{
		  $accept = 'NO';
	      }
	      print STDOUT "match:$rna_id coverage:$score perc_id:$perc_id extent:$extent strand:$strand comment:$label accept:$accept\n";
	      
	      print STDOUT "--------------------\n";
	      
	  }
	  else{
	      ############################################################
	      # we keep anything which is 
	      # within the 2% of the best score
	      # with score >= $min_coverage and percent_id >= $min_perc_id
	      if ( ( $score >= (0.98*$max_score) && 
		     $score >= $min_coverage && 
		     $perc_id >= $min_coverage )
		   ||
		   ( $score >= (0.98*$max_score) &&
		     $score >= (1 + 5/100)*$min_coverage &&
		     $perc_id >= ( 1 - 3/100)*$min_perc_id
		     )
		   ){
	  
		  ############################################################
		  # non-best matches are kept only if they are not unspliced with the
		  # best match being spliced - otherwise they could be processed pseudogenes
		  if ( $reject_potential_pseudos
		       && $count > 1 
		       && $is_spliced 
		       && !is_spliced( $transcript )
		       ){
		      $accept = 'NO';
		  }
		  else{
		      $accept = 'YES';
		      push( @good_matches, $transcript);
		  }
	      }
	      else{
		  $accept = 'NO';
	      }
	      print STDOUT "match:$rna_id coverage:$score perc_id:$perc_id extent:$extent strand:$strand comment:$label accept:$accept\n";
	      
	      print STDOUT "--------------------\n";
	  }
      }
  }
  
  return @good_matches;
}


############################################################
# check whether a transcript is spliced
############################################################


sub is_spliced{
  my ($t) = @_;
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
#
#We want introns of the form:
#    
#    ...###GT...AG###...   ...###AT...AC###...   ...###GC...AG###...
#    
#if we see introns like these:
#    
#    ...###CT...AC###...   ...###GT...AT###...   ...###CT...GC###...
#
#we need to set the strand to the opposite. This can happen when 
#an est/cdna is annotated backwards in the db, if blat reverse 
#complement it to map it, it will find exactily the same exon 
#sequence of an homolog annotated forward, but in the opposite 
#strand. As blat does not reconfirm splice sites like est2genome,
#we need to do it ourselves. Exonerate will do this work for you.
############################################################

sub check_splice_sites{
  my ($transcript, $chr_location) = @_;

  my $verbose = 1;
     
  my $strand = $transcript->start_Exon->strand;
  my @exons;
  if ( $strand == 1 ){
      @exons = sort { $a->start <=> $b->start } @{$transcript->get_all_Exons};
  }
  else{
      @exons = sort { $b->start <=> $a->start } @{$transcript->get_all_Exons};
  }
  my $introns  = scalar(@exons) - 1 ; 
  if ( $introns <= 0 ){
      return $transcript;
  }
  
  my $correct  = 0;
  my $wrong    = 0;
  my $other    = 0;
  
  if ($strand == 1 ){
      
    INTRON:
      for (my $i=0; $i<$#exons; $i++ ){
	  my $upstream_exon   = $exons[$i];
	  my $downstream_exon = $exons[$i+1];
	  
	  my $upstream_start = ($upstream_exon->end     + 1);
	  my $upstream_end   = ($upstream_exon->end     + 2);      
	  my $downstream_start = $downstream_exon->start - 2;
	  my $downstream_end   = $downstream_exon->start - 1;
	  
	  #print STDOUT "upstream $upstream_site, downstream: $downstream_site\n";
	  ## good pairs of upstream-downstream intron sites:
	  ## ..###GT...AG###...   ...###AT...AC###...   ...###GC...AG###.
	  
	  ## bad  pairs of upstream-downstream intron sites (they imply wrong strand)
	  ##...###CT...AC###...   ...###GT...AT###...   ...###CT...GC###...
	  
	  my $upstream_site   = 
	      get_chr_subseq($chr_location, $upstream_start, $upstream_end, $strand );
	  my $downstream_site = 
	      get_chr_subseq($chr_location, $downstream_start, $downstream_end, $strand );
	  
	  unless ( $upstream_site && $downstream_site ){
	      print STDOUT "problems retrieving sequence for splice sites\n";
	      next INTRON;
	  }
	  
	  print STDOUT "strand: + upstream (".
	      ($upstream_start)."-".($upstream_end).") = $upstream_site, downstream ".
	      ($downstream_start)."-".($downstream_end).") = $downstream_site\n" if $verbose;
	  
	  if (  ($upstream_site eq 'GT' && $downstream_site eq 'AG') ||
		($upstream_site eq 'AT' && $downstream_site eq 'AC') ||
		($upstream_site eq 'GC' && $downstream_site eq 'AG') ){
	      $correct++;
	  }
	  elsif (  ($upstream_site eq 'CT' && $downstream_site eq 'AC') ||
		   ($upstream_site eq 'GT' && $downstream_site eq 'AT') ||
		   ($upstream_site eq 'CT' && $downstream_site eq 'GC') ){
	      $wrong++;
	  }
	  else{
	      $other++;
	  }
      } # end of INTRON
  }
  elsif ( $strand == -1 ){
      
      #  example:
      #                                  ------CT...AC---... 
      #  transcript in reverse strand -> ######GA...TG###... 
      # we calculate AC in the slice and the revcomp to get GT == good site
      
    INTRON:
      for (my $i=0; $i<$#exons; $i++ ){
	  my $upstream_exon   = $exons[$i];
	  my $downstream_exon = $exons[$i+1];
	  my $up_site;
	  my $down_site;
	  
	  my $up_start   = $upstream_exon->start - 2;
	  my $up_end     = $upstream_exon->start - 1;
	  my $down_start = $downstream_exon->end + 1;
	  my $down_end   = $downstream_exon->end + 2;
	 
	  my $upstream_site   = 
	      get_chr_subseq($chr_location, $up_start, $up_end, $strand );
	  my $downstream_site = 
	      get_chr_subseq($chr_location, $down_start, $down_end, $strand );
	  
	  unless ( $upstream_site && $downstream_site ){
	      print STDOUT "problems retrieving sequence for splice sites\n";
	      next INTRON;
	  }
	  
	  print STDOUT "strand: + upstream (".
	      ($up_start)."-".($up_end).") = $upstream_site, downstream ".
	      ($down_start)."-".($down_end).") = $downstream_site\n" if $verbose;
	  
	  
	  #print STDOUT "strand: - upstream $upstream_site, downstream: $downstream_site\n";
	  if (  ($upstream_site eq 'GT' && $downstream_site eq 'AG') ||
		($upstream_site eq 'AT' && $downstream_site eq 'AC') ||
		($upstream_site eq 'GC' && $downstream_site eq 'AG') ){
	      $correct++;
	  }
	  elsif (  ($upstream_site eq 'CT' && $downstream_site eq 'AC') ||
		   ($upstream_site eq 'GT' && $downstream_site eq 'AT') ||
		   ($upstream_site eq 'CT' && $downstream_site eq 'GC') ){
	      $wrong++;
	  }
	  else{
	      $other++;
      }
	  
      } # end of INTRON
  }
  unless ( $introns == $other + $correct + $wrong ){
      print STDOUT "STRANGE: introns:  $introns, correct: $correct, wrong: $wrong, other: $other\n";
  }
  if ( $wrong > $correct ){
      print STDOUT "Perhaps this transcript should change strand\n" if $verbose;
  }
}

############################################################

sub get_chr_subseq{
  my ( $chr_location, $start, $end, $strand ) = @_;

  my $command = "chr_subseq $chr_location $start $end |";
 
  #print STDOUT "command: $command\n";
  open( SEQ, $command ) || die("Error running chr_subseq within ExonerateToGenes");
  my $seq = uc <SEQ>;
  chomp $seq;
  close( SEQ );
  
  if ( length($seq) != 2 ){
      print STDOUT "WRONG: asking for chr_subseq $chr_location $start $end and got = $seq\n";
  }
  if ( $strand == 1 ){
      return $seq;
  }
  else{
      ( my $revcomp_seq = reverse( $seq ) ) =~ tr/ACGTacgt/TGCAtgca/;
      return $revcomp_seq;
  }
}


############################################################
#
#    this method changes the strand of the exons
#
sub change_strand{
    my ($transcript) = @_;
    my $original_strand = $transcript->start_Exon->strand;
    my $new_strand      = (-1)*$original_strand;
    foreach my $exon (@{$transcript->get_all_Exons}){
	$exon->strand($new_strand);
	#foreach my $evi ( @{$exon->get_all_supporting_features} ){
	#    $evi->strand($new_strand);
	#    $evi->hstrand( $evi->hstrand*(-1) );
	#}
    }
    $transcript->sort;
    return $transcript;
}

############################################################
























































