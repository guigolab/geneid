#!/usr/local/bin/perl -w

use strict;

use ClusterMerge::GTFTools;
use ClusterMerge::GFFTools;
use GeneComparison::TranscriptComparison;
use Adaptor::SequenceAdaptor;
use Exonerate::ExonerateParser;
use GeneComparison::IntronComparison;
use Bio::Seq;
use Bio::SeqIO;


# file with the correspondance of syntenic regions
my $synteny_file  = $ARGV[0];

my $species = "mouse";

my $verbose = 1;

# where the syntenic sequences are:
#my $mouse_seq_dir = "/projects/encode/data/mm5/";
my $mouse_seq_dir = "/projects/encode/data/msa/genome-test.cse.ucsc.edu/encode/downloads/msa/OCT-2004/";

unless ( $synteny_file ){
    print STDERR "Usage: $0 synteny.file\n";
    exit(0);
}

############################################################
# read synteny
my %synteny;
open (IN,"<$synteny_file") or die("cannot open file $synteny_file");
while(<IN>){
    chomp;
    my ($human,$mouse) = split;
    push (@{$synteny{$human}}, $mouse );
}


my $intron_file = "/projects/encode/computational_predictions/10regions_777introns.gff";

my @introns = ClusterMerge::GFFTools->get_transcripts_from_GFF($intron_file);

my $cut_point;
my $human_seq_dir = "/projects/encode/data/hg16/";
my $human_seq_file;

foreach my $intron (@introns){
    
    my @exons = sort {$a->start <=> $b->start} @{$intron->get_all_Exons};
    next if ( $exons[0]->start < 0 );
     
    my $encode_region = $exons[0]->seqname;
    
    # where the encode sequences are:
    $human_seq_file = $human_seq_dir."/".$encode_region."/".$encode_region.".fa";

    # $transcript is the transcript in mouse
    # $intron is the transcript in human
 
    ############################################################
    # RUN-1
    print "############################################################\n" if $verbose;
    print "RUN 1\n" if $verbose;
    my $options = " -m est2genome ";
    my ($cover,$perc_id,$transcript,$syntenic,$q_start,$q_end,$insertions,$deletions) 
	= get_best_alignment($intron,$encode_region,\%synteny, $options);
    
    print "cover   = $cover\n";
    print "perc_id = $perc_id\n";
    print "syntenic= $syntenic\n";
    print "q_start = $q_start\n";
    print "q_end   = $q_end\n";
    print "insertions = ".values(%$insertions)."\n";
    print "deletions  = ".values(%$deletions)."\n";
    
    unless ($transcript){
	print "NO alignment found for ".$intron->dbID."\n";
	next;
    }
    
    my ($intron_conserved, $mouse_index)
	= check_intron_conservation($intron,$transcript,$insertions,$deletions,$q_start,$q_end);
        
    my $intron_label = "NOT-CONSERVED";
    if ( $intron_conserved == 0 ){
	print "intron not conserved - trying with --forcegtag\n";
	
	############################################################
	# RUN-2
	print "############################################################\n" if $verbose;
	print "RUN 2\n" if $verbose;
	my $options2 = " -m est2genome --forcegtag ";
	my ($cover2,$perc_id2,$transcript2,$syntenic2,$q_start2,$q_end2,$insertions2,$deletions2) 
	    = get_best_alignment($intron,$encode_region,\%synteny,$options2);
	
	my ($intron_conserved2,$mouse_index2);
	if ( $transcript2 ){
	    ($intron_conserved2,$mouse_index2) = check_intron_conservation($intron,$transcript2,$insertions2,$deletions2,$q_start2,$q_end2);
	}
	
	if ( $intron_conserved2 ){ 
	    $intron_label = "CONSERVED";
	    ($cover,$perc_id,$transcript,$syntenic) =
		($cover2,$perc_id2,$transcript2,$syntenic2);
	    $mouse_index = $mouse_index2;
	}
	else{    
	    print "we go back to the first case\n" if $verbose;
	    #    ############################################################
	    #    # RUN-3
	    #    print "############################################################\n" if $verbose;
	    #    print "RUN 3\n" if $verbose;
	    #    my $options3 = " -m coding2genome ";
	    #    my ($cover3,$perc_id3,$transcript3,$syntenic3,$q_start3,$q_end3,$insertions3,$deletions3) 
	    #	= get_best_alignment($intron,$encode_region,\%synteny,$options3);
	    #    
	    #    my $intron_conserved3 
	    #	= check_intron_conservation($intron,$transcript3,$insertions3,$deletions3,$q_start3,$q_end3);
	    #    if ($intron_conserved3 == 1 ){
	    #	($cover,$perc_id,$transcript,$syntenic) =
	    #	    ($cover3,$perc_id3,$transcript3,$syntenic3);
	    #	
	    #    }
	    #    else{
	    #	print "intron conservation not found after 3 runs\n"; 
	    #    }
	}
	
    }
    else{
	$intron_label = "CONSERVED";
    }
    
    ############################################################
    # are both exons fully recovered
    my $label;
    if ( $cover > 0 && $cover < 100 ){
	$label = "PARTIAL";
    }
    elsif( $cover == 100 ){
	$label = "FOUND";
    }
    elsif( $cover == 0 ){
	$label = "MISSING";
    }
    
    ############################################################
    # info about the human intron
    my @human_exons  = sort {$a->start <=> $b->start} @{$intron->get_all_Exons};
    my $human_strand = $human_exons[0]->strand;
    my @h_exons;
    foreach my $he (@human_exons){
	my $st = $he->start."-".$he->end;
	push (@h_exons, $st );
    }
    my $h_exons = join ",",@h_exons;
    my @h = ($encode_region,
	     $intron->dbID,
	     $human_exons[0]->source_tag,
	     "exons:",
	     $h_exons,
	     "h_strand:",
	     $human_strand);
    
    ############################################################
    # info about the mouse transcript
    my @mouse_exons = sort {$a->start <=> $b->start} @{$transcript->get_all_Exons};
    foreach my $e (@mouse_exons){
	print $e->gff_string."\n";
    }
    my @m_exons;
    if (defined $mouse_index){
	for (my $i=$mouse_index; $i<=$mouse_index+1; $i++ ){
	    my $me = $mouse_exons[$i];
	    my $st = $me->start."-".$me->end;
	    push (@m_exons, $st );
	}
    }
    else{
	foreach my $me (@mouse_exons){
	    my $st = $me->start."-".$me->end;
	    push (@m_exons, $st );
	}
    }
    my $m_exons = join ",",@m_exons;
    my $mouse_strand = $mouse_exons[0]->strand;
    my @m = ($syntenic,
	     $transcript->dbID,
	     "exons:",
	     scalar(@mouse_exons),
	     $m_exons,
	     "m_strand",
	     $mouse_strand);
    
    ############################################################
    # number of exons conserved?
    my $conserved;
    if ( scalar(@mouse_exons) == scalar(@human_exons) ){
	$conserved = "enumber-con";
    }
    else{
	$conserved = "enumber-not-con";
    }

    ############################################################
    # print out the result
    my $s = join "\t",("ALIGNMENT",@h,@m,$label,"cover:",$cover,"perc_id:",$perc_id,"intron-pos:",$intron_label,"exon-number:",$conserved);
    print $s."\n";
}



############################################################

sub get_best_alignment{
    my ($intron,$encode_region,$synteny,$options) = @_;

    my %synteny = %$synteny;


    my %cover;
    foreach my $syntenic ( @{$synteny{$encode_region}} ){
	
	my $syntenic_seq_file 
	    = $mouse_seq_dir."/".$encode_region."/".$species.".".$encode_region.".fa";
	
	print "syntenic file: ".$syntenic_seq_file."\n";
	
	############################################################
	# map the intron (exon-pair) to the syntenic region:
	my ($perc_id, $cover, $transcript, $q_start, $q_end, $insertions, $deletions) 
	    = exonerate_transcript_to_syntenic($intron,$syntenic,$syntenic_seq_file,$options);
	
	next unless $transcript;
	push( @{$cover{$cover}}, 
	      [$perc_id,$transcript,$syntenic,$q_start,$q_end,$insertions,$deletions]);
    }
    
    # we store all the result but we only print the best(s)
    # this will keep all the 100% coverage alignments
    # and all the completely missing cases
    # This can be useful when comparing the conservation of
    # exons in orthologous isoforms
    my @ranking = sort {$b <=> $a } keys %cover;

    my $max_cover = $ranking[0];
    
    return unless $cover{$ranking[0]};
    my $result = $cover{$ranking[0]}->[0];
    
    return $max_cover, @{$result};
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

sub exonerate_transcript_to_syntenic{
    my ($tran,$syntenic,$syntenic_file,$options) = @_;

    ############################################################
    # query
    #print "getting transcript sequence from $human_seq_file\n" if $verbose;
    my $tran_string = Adaptor::SequenceAdaptor->get_transcript_seq($human_seq_file,$tran);
    my $tran_seq    = Bio::Seq->new(
				    -DISPLAY_ID => $tran->dbID,
				    -MOLTYPE    => 'dna',
				    -SEQ        => $tran_string,
				    );
    
    #print "tran seq     = ".$tran_string."\n";
    #print "query length = ".$tran_seq->length."\n";
    
    print $tran->dbID."\t".$tran_string."\n";
    my $query_file = '/tmp/query_'.$$.'.fa';
    open ( QUERY, ">$query_file") or die ("cannot open query file $query_file for writing");
    my $query_seqout = Bio::SeqIO->new('-format' => 'Fasta',
				       '-fh'     => \*QUERY);
    $query_seqout->write_seq($tran_seq);
    close(QUERY);

    ############################################################
    # target
    #print "getting genomic sequence\n" if $verbose;
    my $strand       = 1;
    my $syntenic_dna = Adaptor::SequenceAdaptor->get_seq_from_multifasta($syntenic_file,$syntenic);
    my $target_seq   = Bio::Seq->new(
				     -DISPLAY_ID => $syntenic,
				     -MOLTYPE    => 'dna',
				     -SEQ        => $syntenic_dna,
				     );

    my $target = "/tmp/target_".$$.".fa";
    open( TARGET,">$target") || die("Could not open $target $!");
    my $target_seqout = Bio::SeqIO->new('-format' => 'Fasta',
					'-fh'     => \*TARGET);
    $target_seqout->write_seq($target_seq);
    close(TARGET);
    
    ############################################################
    # generate exonerate command
        
    #print "running exonerate\n" if $verbose;

    $options .=   " --terminalrangeint 24 --terminalrangeext 24 --joinrangeint 24 --joinrangeext 24 --spanrangeint 24 --spanrangeext 24 ";
    $options .=   " --showvulgar FALSE ";
    $options .=   " --showalignment TRUE ";
    $options .=   " --bestn 1 ";
    $options .=   " -w 3 ";
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
	print $_."\n" if $verbose;
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
    
    unlink ( $target );
    unlink ($query_file);

    #print "returning query coverage = $cover1\n" if $verbose;
    #print "returning insertions :".values(%$insertions)."\n";
    #print "returning deletions  :".values(%$deletions)."\n";
    return ($perc_id,$cover1, $transcript,$qq_start,$qq_end, $insertions,$deletions);
}



