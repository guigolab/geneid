#!/usr/local/bin/perl -w

use strict;

use ClusterMerge::GTFTools;
use GeneComparison::TranscriptComparison;
use Adaptor::SequenceAdaptor;
use Bio::Seq;
use Bio::SeqIO;


my $human_dir = "/seq/genomes/H.sapiens/golden_path_200307/chromFa/";
my $mouse_dir = "/seq/genomes/M.musculus/golden_path_200405mm5/chromFa/";
my $padding = 5000;

     
#job_10.out:MISSING      15:79889032-79932571.-1 22:37157993-37209228.-1 ENSMUSG00000022429.3    ENSG00000100206.1       ENSMUST00000023065.3:ENSMUST00000080992.1       79929004-79929048.-1 internal
    
while(<>){
    chomp;
    my ($label,
	$mouse_region,
	$human_region,
	$mouse_gene_id,
	$human_gene_id,
	$mouse_transcripts,
	$exon,
	$exon_label) = split;

    
    my ($perc_id,$cover) = align_exon_to_region($exon,$mouse_transcripts,$mouse_region,$mouse_dir,$human_region,$human_dir);
    
    my $s = join "\t", ("ALIGNMENT",$mouse_region,$mouse_gene_id,$mouse_transcripts,$exon,$exon_label,$human_region,$human_gene_id,"cover:",$cover,"percid:",$perc_id);
    print $s."\n";
	
}

    
sub align_exon_to_region{
    my ($exon,$mouse_transcripts,$mouse_region,$mouse_dir,$human_region,$human_dir) = @_;
    
    my $verbose = 1;

    $exon =~/(\d+)\-(\d+)\.(1|-1)/;
    my $e_start  = $1;
    my $e_end    = $2;
    my $e_strand = $3;

    $mouse_region =~/(\S+):(\d+)\-(\d+)\.(1|-1)/;
    my $m_chr   = $1;
    my $m_start = $2;
    my $m_end   = $3;
    my $m_strand= $4;

    $human_region =~/(\S+):(\d+)\-(\d+)\.(1|-1)/;
    my $h_chr   = $1;
    my $h_start = $2;
    my $h_end   = $3;
    my $h_strand= $4; 

    # query
    my $exon_seqname = "chr".$m_chr unless $m_chr=~/chr/;
    my $exon_chr     = $mouse_dir."/".$exon_seqname.".fa";
    my $id1          = $mouse_transcripts."_".$e_start."-".$e_end.".".$e_strand;
    my $string1      = Adaptor::SequenceAdaptor->get_subseq( $exon_chr, $e_start, $e_end, $e_strand);
    my $exon_seq     = Bio::Seq->new(
				     -DISPLAY_ID => $id1,
				     -MOLTYPE    => 'dna',
				     -SEQ        => $string1,
				     );
    
    # target
    my $syntenic = "chr".$h_chr unless $h_chr=~/chr/;
    my $syn_chr  = $human_dir."/".$syntenic.".fa";    
    my $id2      = $h_chr.".".($h_start - $padding)."-".($h_end + $padding).".".$h_strand;
    my $string2  = Adaptor::SequenceAdaptor->get_subseq($syn_chr,$h_start - $padding ,$h_end + $padding ,$h_strand);
    my $target   = Bio::Seq->new(
				   -DISPLAY_ID => $id2,
				   -MOLTYPE    => 'dna',
				   -SEQ        => $string2,
				   );
    
    
    print "comparing $id1 and $id2\n" if $verbose;
    exonerate($exon_seq,$target);
}

sub exonerate{
    my ( $query,$target ) = @_;

    my $verbose = 1;

    ############################################################
    # target
    my $file = 'seq_'.$$.'.fa';
    my $target_file = "/tmp/".$file;
    open( DB_SEQ,">$target_file") || die("Could not open $target_file $!");
    
    my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
                                 '-fh'     => \*DB_SEQ);
    
    $seqout->write_seq($target);
    close( DB_SEQ );
        
    ############################################################
    # create query file
    my $query_file = '/tmp/query_'.$$.'.fa';
    open ( QUERY, ">$query_file") or die ("cannot open query file $query_file for writing");
    
    my $seqout2 = Bio::SeqIO->new('-format' => 'Fasta',
                                  '-fh'     => \*QUERY);
    $seqout2->write_seq($query);
    close(QUERY);
    
    ############################################################
    # generate exonerate command
    
    
    my $options = " --showvulgar FALSE ";
    $options .=   " --showalignment TRUE ";
    $options .=  " --bestn 1 ";
    #$options .=   " -m ungapped:trans ";
    #$options .=  " -m coding2coding ";
    $options .=  " -m est2genome ";
    #$options .=  " -m affine:local ";
    $options .=  " --dnawordlen 3 ";
    $options .=  " --dnahspdropoff 30 ";
    $options .=  " --gapextend -3 ";
    #$options .=  " -m affine:local ";
    $options .=   " --ryo \"RESULT: %S %ql %tl %pi %V\\n\"";
    my $command = "exonerate -t $target_file -q $query_file $options";
    
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
	
	if ($verbose){
	    my $s = join "\t",($q_id,$t_id,$perc_id,$cover1,$cover2);
	    print $s."\n";
	}
    }
    
    unlink ( $target );
    unlink ($query_file);
    
    return ($perc_id,$cover1);
}



