#!/usr/local/bin/perl -w

use strict;

use ClusterMerge::GTFTools;
use GeneComparison::TranscriptComparison;
use Adaptor::SequenceAdaptor;
use Exonerate::ExonerateParser;
use Bio::Seq;
use Bio::SeqIO;


my $human_dir = "/seq/genomes/H.sapiens/golden_path_200307/chromFa/";
my $mouse_dir = "/seq/genomes/M.musculus/golden_path_200405mm5/chromFa/";
my $padding = 5000;

#MISSING region2:        7:128791202-128798650.-1        region1:        11:410653-421869.-1     ENSMUSG00000054662.2    ENSG00000185101.1       ENSMUST00000067836.2:ENSMUST00000047589.3 exon:    128793063-128793210.-1  internal        coding  ENST00000332826.1:410653-421869.-1      cover:  0       perc_id:        0

#PARTIAL region2:        X:62986285-63058865.-1  region1:        X:148575278-148707529.-1        ENSMUSG00000035776.3    ENSG00000102181.6       ENSMUST00000080035.1    exon:   62986292-62988969.-1       internal        half-coding     ENST00000309284.5:148575278-148707647.-1        cover:  15.1232262882748        perc_id:        72.35

    
my $verbose = 0;

############################################################
# we cache the chromosome lengths
my %chr_length;

my %seen;
while(<>){
    chomp;
    my ($label,
	$mr,
	$mouse_region,
	$hr,
	$human_region,
	$mouse_gene_id,
	$human_gene_id,
	$mouse_transcripts,
	$e,
	$exon,
	$position,
	$coding_type,
	$othologous_transcript) = split;
    
    next if $seen{$mouse_gene_id}{$mouse_region}{$exon}{$human_gene_id};

    my ($perc_id,$cover, $transcript) = align_exon_to_region($exon,$mouse_transcripts,$mouse_region,$mouse_dir,$human_region,$human_dir);
    $seen{$mouse_gene_id}{$mouse_region}{$exon}{$human_gene_id}++;

    if ($transcript){
	
	if ($verbose){
	    foreach my $e (@{$transcript->get_all_Exons}){
		print $e->gff_string."\n";
	    }
	}
    
    
my @remaped = remap_exons($transcript, $human_region);
	foreach my $ex (@remaped){
	    print $ex->gff_string."\n";
	}
    }
    

    my $s = join "\t", ("ALIGNMENT",$label,$mouse_region,$mouse_gene_id,$mouse_transcripts,"exon:",$exon,$position,$coding_type,$human_region,$human_gene_id,"cover:",$cover,"percid:",$perc_id);
    print $s."\n";
}

    
sub align_exon_to_region{
    my ($exon,$mouse_transcripts,$mouse_region,$mouse_dir,$human_region,$human_dir) = @_;
    
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
    
    unless( $chr_length{$syn_chr} ){
	$chr_length{$syn_chr} = Adaptor::SequenceAdaptor->get_seq_length($syn_chr);
    }

    ############################################################
    # check that we do not go over the edge
    my ($target_start,$target_end);
    if ( $h_start - $padding < 1 ){
	$target_start = 1;
    }
    else{
	$target_start = $h_start - $padding;
    }
    if ( $h_end + $padding > $chr_length{$syn_chr} ){
	$target_end = $chr_length{$syn_chr};
    }
    else{
	$target_end   = $h_end + $padding;
    }
    my $string2  = Adaptor::SequenceAdaptor->get_subseq($syn_chr, $target_start, $target_end, $h_strand);
    my $target   = Bio::Seq->new(
				 -DISPLAY_ID => $id2,
				 -MOLTYPE    => 'dna',
				 -SEQ        => $string2,
				 );
    print "comparing $id1 and $id2\n" if $verbose;
    return exonerate($exon_seq,$target);
}

sub exonerate{
    my ( $query,$target ) = @_;


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
    my $transcript;
    while(<OUT>){
        chomp;
	next if /exonerate/;
	print $_."\n" if $verbose;
	next unless /RESULT/;
	
	$transcript = Exonerate::ExonerateParser->parse_output($_) unless $transcript;
	
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
    
    return ($perc_id,$cover1,$transcript);
}


############################################################

sub remap_exons{
    my ($transcript,$human_region) = @_;
    
    #11.435280-486398.1      cDNA    exon    45141   45292   336     +       .       ENSMUST00000026568.2_128845211-128846031.1.22 ...
    
    
    $human_region =~/(\S+):(\d+)\-(\d+)\.(1|-1)/;
    my $chr          = $1;
    my $slice_start  = $2 - $padding;
    my $slice_end    = $3 + $padding;
    my $slice_strand = $4; 

    print "human chr=$chr start=$slice_start end=$slice_end strand=$slice_strand\n" if $verbose;
    
    my @exons = @{$transcript->get_all_Exons};
    foreach my $e ( @exons ){
	$e->source_tag("remapped");
	if ( $e->strand == 1 
	     &&
	     $slice_strand == 1 
	     ){
	    #print "slice_start = $slice_start\n";
	    #print "exon start  = ".$e->start."\n";
	    #print "exon end    = ".$e->end."\n";
	    
	    $e->start( $e->start + $slice_start - 1 );
	    $e->end(   $e->end   + $slice_start - 1 );
	    #print " new start  = ".$e->start."\n";
	    #print " new end    = ".$e->end."\n";


	}
	elsif( $e->strand == 1 
	    &&
	    $slice_strand == -1
	    ){
	    my $slice_length = $slice_end - $slice_start + 1;
	    my $end   = $e->end;
	    my $start = $e->start;

	    $e->start( $slice_start + ($slice_length - $end   + 1) - 1 );
	    $e->end(   $slice_start + ($slice_length - $start + 1) - 1 );
	    $e->strand(-1);
	}
	elsif ( $e->strand == -1 
		&&
		$slice_strand == 1 
		){
	    $e->start( $e->start + $slice_start - 1 );
	    $e->end(   $e->end   + $slice_start - 1 );
	}
	elsif ($e->strand == -1 
	       &&
	       $slice_strand == -1
	       ){
	    my $slice_length = $slice_end - $slice_start + 1;
	    my $end   = $e->end;
	    my $start = $e->start;
	    
	    $e->start( $slice_start + ($slice_length - $end   + 1) - 1 );
	    $e->end(   $slice_start + ($slice_length - $start + 1) - 1 );
	    $e->strand(-1);
	}
    }
    return @exons;

}
    
