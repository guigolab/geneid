#!/usr/local/bin/perl -w

use strict;

use ClusterMerge::GTFTools;
use GeneComparison::TranscriptComparison;
use Adaptor::SequenceAdaptor;
use Bio::Seq;
use Bio::SeqIO;


######
# script to map introns from one species
# into the syntenic region in another species

# we need to calculate first for
# each gene (or gene region) the
# syntenic region in teh other species

# the prototype here is the set of vega genes
# that we want to map to the mouse syntenic regions

my $human_id = $ARGV[0];
my $mouse_region = $ARGV[1];

my $human_dir = "/seq/genomes/H.sapiens/golden_path_200307/chromFa/";
my $mouse_dir = "/seq/genomes/M.musculus/golden_path_200405mm5/chromFa/";

# provide a file with the GTF for the human genes
my $human_file = "/home/ug/eeyras/Projects/ATD/ComparativeGenomics/HumanMouse/AllExons/noNT.human.gtf";

unless ( $human_id && $mouse_region ){
    print STDERR "Usage: $0 human.id mouse-region (region.start.end.strand)\n";
    exit(0);
}

my ($human_genes,$human_gene_hash) = ClusterMerge::GTFTools->get_gene_from_GTF($human_file,$human_id);

my $human_gene = $human_gene_hash->{$human_id};

unless ($human_gene){
    print print "CANNOT FIND $human_id\n";
    exit(0);
}


# get all the transcripts from this geen and pass the
# list of transcripts to IntronComparison::get_uniq_introns



#my ($exons1, $e2t1, $pos1) = $gene1->get_Exons_and_positions;
#my ($exons2, $e2t2, $pos2) = $gene2->get_Exons_and_positions;
#my ($exons1, $pos1) = $gene1->get_unique_Exons_and_positions;
#my ($exons2, $pos2) = $gene2->get_unique_Exons_and_positions;

# gene1 is human
#my ($exons1, $pos1, $cod_flag1) = $gene1->get_unique_Exons_and_positions_and_coding_flag;
# for human we use the transcripts;
my ($exons1, $e2t1) = $gene1->get_Exons;

# gene2 is mouse
my ($exons2, $pos2, $cod_flag2) = $gene2->get_unique_Exons_and_positions_and_coding_flag;

my @exons1 = sort { $a->start <=> $b->start } @$exons1;
my $seq1    = $exons1[0]->seqname;
my $start1  = $exons1[0]->start;
my $end1    = $exons1[-1]->end;
my $strand1 = $exons1[0]->strand;
my $region1 = $seq1.":".$start1."-".$end1.".".$strand1;

my @exons2 = sort { $a->start <=> $b->start } @$exons2;
my $seq2    = $exons2[0]->seqname;
my $start2  = $exons2[0]->start;
my $end2    = $exons2[-1]->end;
my $strand2 = $exons2[0]->strand;
my $region2 = $seq2.":".$start2."-".$end2.".".$strand2;

############################################################
# check whether all exons in gene2 are also in gene1:

my @found;
my @missing;
my %adj_list;
my %cover;
my %perc_id;

foreach my $e2 (@exons2){
    
    # we compare against each transcript in gene1, and choose the one(s) that contains ht most
    # of exon2
    my %cover;
    foreach my $tran (@{$gene1->get_Transcripts}){
	
	my ($perc_id, $cover, $exon_tran, $exon_locus, $tran_id) = exonerate_e2t($e2,$tran,$dir2,$dir1);
	
	# let's print out all of the results
	# we will filter them afterwards

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
	
	my $e_id = $exon_tran."\t"."exon:\t".$exon_locus."\t".$pos2->{$e2}."\t".$cod_flag2->{$e2};
	
	my $s = join "\t", ($label,"region2:",$region2,"region1:",$region1,$id2,$id1,$e_id,$tran_id,"cover:",$cover,"perc_id:",$perc_id);
	push( @{$cover{$cover}}, $s );
    }
    
    # we store all the result but we only print the best(s)
    # this will keep all the 100% coverage alignments
    # and all the completely missing cases
    # This can be useful when comparing the conservation of
    # exons in orthologous isoforms
    my @ranking = sort {$b <=> $a } keys %cover;
    foreach my $result ( @{$cover{$ranking[0]}} ){
	print $result."\n";
    }
}


############################################################

sub exonerate_e2t{
    my ( $exon,$tran,$exon_dir,$tran_dir) = @_;

    my $verbose = 0;
    
    my $exon_seqname = $exon->seqname;
    $exon_seqname    = "chr".$exon_seqname unless $exon_seqname=~/chr/;
    my $exon_chrfile = $exon_dir."/".$exon_seqname.".fa";

    my @exons2       = sort {$a->start <=> $b->end} @{$tran->get_all_Exons};
    my $tran_seqname = $exons2[0]->seqname;
    $tran_seqname    = "chr".$tran_seqname unless $tran_seqname=~/chr/;
    my $tran_chrfile = $tran_dir."/".$tran_seqname.".fa";
    
    # query
    my $exon_tran  = $exon->transcript_tag;
    my $exon_locus = $exon->start."-".$exon->end.".".$exon->strand;
    my $exon_id = $exon_tran.";".$exon_locus;
    #target
    my $tran_id = $tran->dbID.":".$exons2[0]->start."-".$exons2[-1]->end.".".$exons2[0]->strand;
    
    print "comparing $exon_id and $tran_id\n" if $verbose;
    
    my ($exon_seq, $tran_seq );
    my $exon_string = Adaptor::SequenceAdaptor->get_exon_seq($exon_chrfile,$exon);
    $exon_seq       = Bio::Seq->new(
				    -DISPLAY_ID => $exon_id,
				    -MOLTYPE    => 'dna',
				    -SEQ        => $exon_string,
				    );
    
    my $tran_string = Adaptor::SequenceAdaptor->get_transcript_seq($tran_chrfile,$tran);
    $tran_seq       = Bio::Seq->new(
				    -DISPLAY_ID => $tran_id,
				    -MOLTYPE    => 'dna',
				    -SEQ        => $tran_string,
				    );
    
    ############################################################
    # target
    
    my $target = "/tmp/target_".$$.".fa";
    open( DB_SEQ,">$target") || die("Could not open $target $!");
    
    my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
                                 '-fh'     => \*DB_SEQ);
    
    $seqout->write_seq($tran_seq);
    close( DB_SEQ );
        
    ############################################################
    # create query file
    my $query_file = '/tmp/query_'.$$.'.fa';
    open ( QUERY, ">$query_file") or die ("cannot open query file $query_file for writing");
    
    my $seqout2 = Bio::SeqIO->new('-format' => 'Fasta',
                                  '-fh'     => \*QUERY);
    $seqout2->write_seq($exon_seq);
    close(QUERY);
    
    ############################################################
    # generate exonerate command
        
    my $options = " --showvulgar FALSE ";
    $options .=   " --showalignment TRUE ";
    $options .=   " --bestn 1 ";
    $options .=   " -m affine:local ";
    $options .=   " -w 3 ";
    $options .=   " --dnahspdropoff 30 ";
    $options .=   " --ryo \"RESULT: %S %ql %tl %pi %V\\n\"";
    my $command = "exonerate -t $target -q $query_file $options";
    
    open( OUT, "$command | ");
    #RESULT: ENST00000335466.2:37369656-37369745.1 0 90 + ENSMUST00000046463.2:80047601-80047687.1 0 87 + 380 90 87 96.55 M 56 56 G 3 0 M 31 31
    my $perc_id = 0;
    my $cover1  = 0;
    my $cover2  = 0;
    while(<OUT>){
        chomp;
	next if /exonerate/;
	print $_."\n" if $verbose;
	next unless /RESULT/;
	
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
	if ($verbose){
	    my $s = join "\t",($q_id,$t_id,$perc_id,$cover1,$cover2);
	    print $s."\n";
	}
    }
    
    unlink ( $target );
    unlink ($query_file);



    #my $cover = int(100*($cover1 + $cover2)/2)/100;
 
    print "returning query coverage = $cover1\n" if $verbose;
    return ($perc_id,$cover1, $exon_tran, $exon_locus,$tran_id);
}



