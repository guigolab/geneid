#!/usr/local/bin/perl -w

use strict;

use ClusterMerge::GTFTools;
use ClusterMerge::GFFTools;
use GeneComparison::TranscriptComparison;
use Adaptor::SequenceAdaptor;
use Exonerate::ExonerateParser;
use GeneComparison::IntronComparison;
use GeneComparison::SpliceSiteComparison;
use Bio::Seq;
use Bio::SeqIO;


# mouse gene id to be aligned:
my $mouse_gene_id = $ARGV[0];

# file with the mouse genes in GTF format
my $mouse_file = "/home/ug/eeyras/Projects/ATD/AltSplicingEvaluation/ComparativePredictions/mouse.genes.gtf"; # version NCBIm33

my $verbose = 1;

# where the human sequence is
my $human_seq_file = "/seq/genomes/H.sapiens/golden_path_200307/chromFa/chr10.fa";

# where the mouse_seq is
my $mouse_seq_dir  = "/seq/genomes/M.musculus/golden_path_200405mm5/chromFa/";

# we use Ensembl gene build for the NCBI m33 mouse assembly (freeze May 27, 2004, strain C57BL/6J).

unless ( $mouse_gene_id ){
    print STDERR "Usage: $0 mouse.gene.id\n";
    exit(0);
}

print "getting gene $mouse_gene_id from $mouse_file\n" if $verbose;
my ($mouse_genes,$gene_hash) = ClusterMerge::GTFTools->get_gene_from_GTF($mouse_file,$mouse_gene_id);

my $mouse_gene = $mouse_genes->[0];
unless ($mouse_gene){
    print STDERR "Could not find gene $mouse_gene_id in $mouse_file\n";
    exit(0);
}

my @transcripts = @{$mouse_gene->get_Transcripts};

foreach my $t_in_mouse (@transcripts){
  
    ############################################################
    # don't run if the mouse transcript is non-coding or
    # it is in a random sequence:
    next unless ($t_in_mouse->coding_transcript);
    my @mexons = @{$t_in_mouse->get_all_Exons};
    my $num_exons = scalar(@mexons);
    next if ( $mexons[0]->seqname =~/random/ );

    # get rid of frameshifts:
    
    $t_in_mouse = ClusterMerge::TranscriptUtils->_difuse_small_introns($t_in_mouse, 10);
    @mexons = @{$t_in_mouse->get_all_Exons};
    
    my $label = '';
    unless ( $num_exons == scalar(@mexons) ){
	$label = "\tframeshift";
    }
    # mouse GFF
    print "TESTING\t$mouse_gene_id\t".$t_in_mouse->dbID."\n";
    foreach my $e ( @mexons ){
	print "MOUSEGFF\t".$e->gff_string."\t$mouse_gene_id".$label."\n";
    }


    ############################################################
    # RUN-1
    print "############################################################\n" if $verbose;
    print "RUN 1\n" if $verbose;
    my $options = " -m est2genome ";
    my ($perc_id,$cover,$t_in_human,$q_start,$q_end, $insertions,$deletions)
	= get_best_alignment($t_in_mouse,$mouse_seq_dir,$human_seq_file, $options);
    
    print "cover   = $cover\n";
    print "perc_id = $perc_id\n";
    print "q_start = $q_start\n";
    print "q_end   = $q_end\n";
    print "insertions = ".values(%$insertions)."\n";
    print "deletions  = ".values(%$deletions)."\n";
    
    my ($avescore,$maxscore,$comment) = check_splice_site_conservation($t_in_human,$t_in_mouse);
    
    my @hexons = @{$t_in_human->get_all_Exons};
    
    print "SCORE\tAVERAGE\t$avescore\tMAX\t$maxscore".
	"\t$mouse_gene_id\t".$t_in_mouse->dbID.
	"\tmexons:\t".scalar(@mexons)."\thexons:\t".scalar(@hexons).
	"\tcover:\t".$cover."\tperc_id\t".$perc_id.
	"\t$comment\t".
	"\n";
    
    # full transcript:
    foreach my $e ( @hexons ){
	unless ( $e->transcript_tag ){
	    $e->transcript_tag( $t_in_mouse->dbID."\t".$mouse_gene_id );
	}
	print "HUMANGFF\t".$e->gff_string.
	    "\tAVERAGE\t$avescore\tMAX\t$maxscore".
	    "\tmexons:\t".scalar(@mexons)."\thexons:\t".scalar(@hexons).
	    "\tcover:\t".$cover."\tperc_id\t".$perc_id.
	    "\t$comment\t".
	    "\n";
    }
}

############################################################

sub check_splice_site_conservation{
    my ($ht,$mt) = @_;

    my @hexons = sort { $a->start <=> $b->start} @{$ht->get_all_Exons};
    my @mexons = sort { $a->start <=> $b->start} @{$mt->get_all_Exons};
    
    if ( $hexons[0]->strand == -1 ){
	@hexons = sort { $b->start <=> $a->start} @{$ht->get_all_Exons};
    }
    if ( $mexons[0]->strand == -1 ){
	@mexons = sort { $b->start <=> $a->start} @{$mt->get_all_Exons};
    }
    
    my $score_matrix = "/home/ug/eeyras/Projects/ATD/AltSplicingEvaluation/ComparativePredictions/Splice.Sites.Conservation/orthosites.validated.hm.smatrices";
    
    my $mousechr;
    if ( $mexons[0]->seqname =~/chr/ ){
	$mousechr = $mexons[0]->seqname.".fa";
    }
    else{
	$mousechr = "chr".$mexons[0]->seqname.".fa";
    }
    my $human_seqfile = $human_seq_file;
    my $mouse_seqfile = $mouse_seq_dir."/".$mousechr;
    
    my $average_score = 0;
    my $max_score     = -9999;

    my $comment ='';
    if (scalar(@hexons)==scalar(@mexons)){
	
	if ( scalar(@hexons) == 1 ){
	    
	    $comment = "single-exon-conserved";
	}
	else{
	    for(my $i=0; $i<scalar(@mexons)-1; $i++ ){
		
		# donors
		my $md = Adaptor::SequenceAdaptor->get_donor_seq($mouse_seqfile,$mexons[$i],3,7);
		my $hd = Adaptor::SequenceAdaptor->get_donor_seq($human_seqfile,$hexons[$i],3,7);
		
		# acceptors
		my $ma = Adaptor::SequenceAdaptor->get_acceptor_seq($mouse_seqfile,$mexons[$i+1],3,20);
		my $ha = Adaptor::SequenceAdaptor->get_acceptor_seq($human_seqfile,$hexons[$i+1],3,20);
		
		print "intron $i\n";
		print "md:$md\tma:$ma\n";
		print "hd:$hd\tha:$ha\n";
	    
		my ($donor_score,$acc_score) = GeneComparison::SpliceSiteComparison->score_orthologous_splice_sites($score_matrix,$hd,$ha,$md,$ma);
		#print "SCORE: d-score:$donor_score a-score:$acc_score\n";
		$average_score += ( $donor_score + $acc_score );
		
		if ( ($donor_score + $acc_score) > $max_score ){
		    $max_score = ( $donor_score + $acc_score );
		}
	    }
	    $comment = "conserved";
	    $average_score = $average_score/(scalar(@hexons)-1);
	}
    }
    elsif( scalar(@hexons) == 1 && scalar(@mexons)>1 ){
	$comment = "split-in-mouse";
    }
    elsif( scalar(@hexons)>1 && scalar(@mexons) == 1 ){
	$comment = "split-in-human";
    }
    else{	
	my @matrix;
	my $initial = 0;
	# get a simple greedy alignment:
	for(my $i=0; $i<scalar(@mexons)-1; $i++ ){
	    
	    my $md = Adaptor::SequenceAdaptor->get_donor_seq($mouse_seqfile,$mexons[$i],3,7);
	    my $ma = Adaptor::SequenceAdaptor->get_acceptor_seq($mouse_seqfile,$mexons[$i+1],3,20);
	    
	    for(my $j=$initial; $j<scalar(@hexons)-1; $j++){

		my $hd = Adaptor::SequenceAdaptor->get_donor_seq($human_seqfile,$hexons[$i],3,7);
		my $ha = Adaptor::SequenceAdaptor->get_acceptor_seq($human_seqfile,$hexons[$i+1],3,20);
		
		#print "intron $i with intron $j\n";
		#print "md:$md\tma:$ma\n";
		#print "hd:$hd\tha:$ha\n";
		my ($donor_score,$acc_score) = GeneComparison::SpliceSiteComparison->score_orthologous_splice_sites($score_matrix,$hd,$ha,$md,$ma);
		# mouse-human matrix
		$matrix[$i][$j] = ( $donor_score + $acc_score );
	    }   
	}
	
	# the one with least introns should at least have all its introns corresponding to the other one
	#                   i
	#  mouse    ######------####
	#            / \          \
	#  human  ###---###-------#####
	#                     j
	#
	# For each intron in the one with least introns
	# we choose the best pair (available) in the ortholog
	
	
	if( scalar(@hexons) < scalar(@mexons) ){
	    for(my $j=0; $j<scalar(@hexons)-1; $j++){
		my $max = max_from_column(\@matrix,$j);
		$average_score += $max;
		if ( $max > $max_score ){
		    $max_score = $max;
		}
	    }
	    $comment       = "intron-lost-in-human";
	    $average_score = $average_score/( scalar(@hexons) - 1);
	}
	else{
	    for(my $i=0; $i<scalar(@mexons)-1; $i++){
		my $max = max_from_row(\@matrix,$i);
		$average_score += $max;
		if ( $max > $max_score ){
		    $max_score = $max;
		}
	    }
	    $comment       = "intron-gained-in-human";
	    $average_score = $average_score/(scalar(@mexons) - 1);
	}
    }
    return ($average_score,$max_score,$comment);
}

############################################################

sub max_from_row{
    my ($matrix,$row) = @_;
    my @list;
    my $j=0;
    while ( $matrix->[$row]->[$j] ){
	push (@list,$matrix->[$row]->[$j]);
	$j++;
    }
    return max( @list );
}

############################################################
sub max_from_column{
    my ($matrix,$column) = @_;
    my @list;
    my $i=0;
    while ( $matrix->[$i]->[$column] ){
	push (@list,$matrix->[$i]->[$column]);
	$i++;
    }
    return max( @list );
}

############################################################

sub max{
    my ($max,@rest) = @_;
    foreach my $r (@rest){
	$max = $r if ( $r > $max );
    }
    return $max;
}

############################################################

sub get_best_alignment{
    my ($mouse_t,$mouse_seq_dir,$human_seq_file,$options) = @_;
    
    my $mouse_chr;
    my @exons = @{$mouse_t->get_all_Exons};
    my $seqname = $exons[0]->seqname;
    if ( $seqname =~/chr/ ){
	$mouse_chr = $seqname;
    }
    else{
	$mouse_chr = "chr".$seqname;
    }
    my $mouse_seq_file = $mouse_seq_dir."/".$mouse_chr.".fa";
    
    ############################################################
    # query
    #print "getting transcript sequence from $human_seq_file\n" if $verbose;
    my $tran_string = Adaptor::SequenceAdaptor->get_transcript_seq($mouse_seq_file,$mouse_t);
    my $tran_seq    = Bio::Seq->new(
				    -DISPLAY_ID => $mouse_t->dbID,
				    -MOLTYPE    => 'dna',
				    -SEQ        => $tran_string,
				    );
    
    #print "tran seq     = ".$tran_string."\n";
    #print "query length = ".$tran_seq->length."\n";
    
    print $mouse_t->dbID."\t".$tran_string."\n";
    my $query_file = '/tmp/query_'.$$.'.fa';
    open ( QUERY, ">$query_file") or die ("cannot open query file $query_file for writing");
    my $query_seqout = Bio::SeqIO->new('-format' => 'Fasta',
				       '-fh'     => \*QUERY);
    $query_seqout->write_seq($tran_seq);
    close(QUERY);
    
    ############################################################
    # target
    my $target = $human_seq_file;

    #print "getting genomic sequence\n" if $verbose;
    #my $strand       = 1;
    #my $syntenic_dna = Adaptor::SequenceAdaptor->get_seq_from_multifasta($syntenic_file,$syntenic);
    #my $target_seq   = Bio::Seq->new(
#				     -DISPLAY_ID => $syntenic,
#				     -MOLTYPE    => 'dna',
#				     -SEQ        => $syntenic_dna,
#				     );
#    my $target = "/tmp/target_".$$.".fa";
#    open( TARGET,">$target") || die("Could not open $target $!");
#    my $target_seqout = Bio::SeqIO->new('-format' => 'Fasta',
#					'-fh'     => \*TARGET);
#    $target_seqout->write_seq($target_seq);
#    close(TARGET);
    
    ############################################################
    # generate exonerate command
        
    #print "running exonerate\n" if $verbose;
    $options .=   " --terminalrangeint 24 --terminalrangeext 24 --joinrangeint 24 --joinrangeext 24 --spanrangeint 24 --spanrangeext 24 ";
    $options .=   " --showvulgar FALSE ";
    $options .=   " --showalignment TRUE ";
    $options .=   " --bestn 1 ";
    #$options .=   "  -w 3 ";
    $options .=   " --softmasktarget TRUE ";
    $options .=   " --dnahspdropoff 30 ";
    $options .=   " --ryo \"RESULT: %S %ql %tl %pi %V\\n\"";
    my $command = "exonerate -t $target -q $query_file $options";
    
    open( OUT, "$command | ");
    #RESULT: ENST00000335466.2:37369656-37369745.1 0 90 + ENSMUST00000046463.2:80047601-80047687.1 0 87 + 380 90 87 96.55 M 56 56 G 3 0 M 31 31
    my $perc_id = 0;
    my $cover1  = 0;
    my $cover2  = 0;

    my $transcript;
    my $insertions;
    my $deletions;

    my $qq_start = 0;
    my $qq_end   = 0;
    while(<OUT>){
        chomp;
	next if /exonerate/;
	print "EXO\t".$_."\n" if $verbose;
	next unless /RESULT/;
		
	($transcript, $insertions, $deletions) 
	    = Exonerate::ExonerateParser->parse_output($_) unless $transcript;
	
	#print "exonerate result:\n";
	#ClusterMerge::TranscriptUtils->_print_SimpleTranscript($transcript);
	my ( $tag, $q_id, $q_start, $q_end, $q_strand, $t_id, $t_start, $t_end, $t_strand, $score,  $q_length, $t_length, $p_id, @blocks) = split;
	
	if ( $p_id > $perc_id ){
	    $perc_id = $p_id;
	}
	if ( 100*($q_end - $q_start)/$q_length > $cover1 ){
	    $cover1 = 100*($q_end - $q_start)/$q_length;
	}
	if ( 100*($t_end - $t_start)/$t_length > $cover2 ){
	    $cover2 = 100*($t_end - $t_start)/$t_length;
	}
	#if ($verbose){
	#    my $s = join "\t",($q_id,$t_id,$perc_id,$cover1,$cover2);
	#    print $s."\n";
	#}
	
	unless ($qq_start && $qq_end ){
	    $qq_start = $q_start + 1;
	    $qq_end   = $q_end;
	}

    }
    
    unlink ($query_file);

    #print "returning query coverage = $cover1\n" if $verbose;
    #print "returning insertions :".values(%$insertions)."\n";
    #print "returning deletions  :".values(%$deletions)."\n";
    return ($perc_id,$cover1,$transcript,$qq_start,$qq_end, $insertions,$deletions);
}




############################################################

sub check_intron_conservation{
    my ($human,$mouse,$insertions,$deletions,$q_start,$q_end) = @_;
    
    
    ############################################################
    # check whether the intron is conserved:
    # 
    #  1        p         p+1           L
    #  ##########-----------#############
    #                 |
    #                \|/
    #         s               e     s>=1 && e<=L
    #         #################     (s,e) = piece of the transcript aligned
    #
    #    ####-##-##---------###-###### query transcript
    #    ##-########################## genomic
    #                 |
    #                \|/  we must take into account the gaps
    #
    #     #########---------##########     
    #        l1                 l2
    #
    # if 
    # we have aligned the donor site of the first exon, 
    # i.e., (l1 + gaps_in_genomic - gaps_in_query + s - 1) == p
    # AND the second exon is aligned 
    # i.e., (e > p+1)
    # then we know the exon-exon boundaries (positions) are conserved.
    
    my %insertions = %$insertions;  # gap in the genomic (G 2 0), etc...
    my %deletions  = %$deletions;   # gap in the exon (G 0 1), etc....

    my @human_exons = sort {$a->start <=> $b->start} @{$human->get_all_Exons};
    my @mouse_exons = sort {$a->start <=> $b->start} @{$mouse->get_all_Exons};
    
    my $left_length  = $human_exons[0]->end  - $human_exons[0]->start + 1;
    my $right_length = $human_exons[-1]->end - $human_exons[-1]->start + 1;

    my $cut_point = $left_length;
    
    # initial value:
    my $intron_conserved = 0;
    
    my $exon_index;
    if (scalar(@mouse_exons)>1){
	
	############################################################
	for(my $i=0; $i<scalar(@mouse_exons) - 1; $i++){
	    
	    #print "******************8 IN LOOP\n";

	    $exon_index = $i;
	    #print "exon index = ".$exon_index."\n";
	    my $left_mouse_exon  = $mouse_exons[$i];
	    my $right_mouse_exon = $mouse_exons[$i+1];
	    
	    ############################################################
	    # which mouse exon I have to look at depends on whether
	    # upon cross-species alignment the strand is flipped or not
	    my $insertions = $insertions{$left_mouse_exon};
	    my $deletions  = $deletions{$left_mouse_exon};
	    my $left_exon_length = ($left_mouse_exon->end - $left_mouse_exon->start + 1);
	    
	    if ($left_mouse_exon->strand != $human_exons[0]->strand){
		$cut_point = $right_length;
	    }
	    print "-------------------------------\n";
	    print "l1        = ".$left_exon_length."\n";
	    print "cut_point = ".$cut_point."\n";
	    print "insertions= ".$insertions."\n";
	    print "deletions = ".$deletions."\n";
	    print "q start   = ".$q_start."\n";
	    print "q end     = ".$q_end."\n";
	    print "-------------------------------\n";
	    
	    if ( abs( ($left_exon_length 
		       + $insertions
		       - $deletions 
		       + $q_start - 1 ) 
		      - $cut_point ) < 2
		 &&
		 $q_end > $cut_point+1
		 ){
		$intron_conserved = 1;
	    }
	    
	    if ( $intron_conserved ){
		#$exon_index = $i;
		#print "before last: exon index = ".$exon_index."\n";
		   
		last;
	    }
	    else{
		$q_start = $q_start + $left_exon_length - 1 - $deletions;
	    }
	}
    }
    #if (defined($exon_index)){
#	print "returning exon index $exon_index\n";
#    }
#    else{
#	print "No exon index, number of mouse exons = ".scalar(@mouse_exons)."\n";
#    }
    return ($intron_conserved,$exon_index);
}

############################################################


