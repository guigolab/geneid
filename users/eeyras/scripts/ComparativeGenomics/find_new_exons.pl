#!/usr/local/bin/perl -w

use strict;

use ClusterMerge::GTFTools;
use GeneComparison::TranscriptComparison;
use Adaptor::SequenceAdaptor;
use Bio::Seq;
use Bio::SeqIO;


my $human_file = $ARGV[0];
my $mouse_file = $ARGV[1];
my $orthology_file = $ARGV[2];

my $dir1 = "/seq/genomes/H.sapiens/golden_path_200307/chromFa/";
my $dir2 = "/seq/genomes/M.musculus/golden_path_200405mm5/chromFa/";

unless ( $human_file && $mouse_file && $orthology_file ){
    print STDERR "Usage: $0 human.gtf mouse.gtf human.mouse.orthology\n";
    exit(0);
}

my ($human_genes,$human_gene_hash) = ClusterMerge::GTFTools->get_genes_from_GTF($human_file);
my ($mouse_genes,$mouse_gene_hash) = ClusterMerge::GTFTools->get_genes_from_GTF($mouse_file);


# get the orthologous pairs:
open (ORT,"<$orthology_file") or die("cannot open $orthology_file");
my %t2g;

my %e2t;

while(<ORT>){
    chomp;
    my @e = split;
    my $id1 = $e[0];
    my $id2 = $e[1];
    
    next unless ($id1=~/ENS/ && $id2=~/ENS/);
    
    my $gene1 = $human_gene_hash->{$id1};
    my $gene2 = $mouse_gene_hash->{$id2};
    
    next unless ($gene1 && $gene2);
    
    print "ORTHOLOGOUS\t".$gene2->dbID."\t".$gene1->dbID."\n";
    #my ($exons1, $e2t1, $pos1) = $gene1->get_Exons_and_positions;
    #my ($exons2, $e2t2, $pos2) = $gene2->get_Exons_and_positions;
    my ($exons1, $pos1) = $gene1->get_unique_Exons_and_positions;
    my ($exons2, $pos2) = $gene2->get_unique_Exons_and_positions;
    
    my @exons1 = sort { $a->start <=> $b->start } @$exons1;
    my @exons2 = sort { $a->start <=> $b->start } @$exons2;
    
    my $seq1    = $exons1[0]->seqname;
    my $start1  = $exons1[0]->start;
    my $end1    = $exons1[-1]->end;
    my $strand1 = $exons1[0]->strand;
    my $region1 = $seq1.":".$start1."-".$end1.".".$strand1;

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
	
	foreach my $e1 (@exons1){
	    
	    my ($perc_id, $cover) = exonerate_Exons($e1,$e2,$dir1,$dir2);
	    if ($cover >= 50 ){
		push ( @{$adj_list{$e2}}, $e1 );
		$cover{$e1}{$e2} = $cover;
		$perc_id{$e1}{$e2} = $perc_id;
	    }
	}
	
	my $e_id2 = $e2->transcript_tag."\t".$e2->start."-".$e2->end.".".$e2->strand."\t".$pos2->{$e2};
	
	if ( $adj_list{$e2} && @{$adj_list{$e2}} ){
	    
	    my $label;
	    $label = "PAIR" if (  @{$adj_list{$e2}} == 1 );
	    $label = "MULTI-PAIR" if (  @{$adj_list{$e2}} > 1 );
	    
	    foreach my $pair ( @{$adj_list{$e2}} ){
		
		my $e_id1 = $pair->transcript_tag."\t".$pair->start."-".$pair->end.".".$pair->strand."\t".$pos1->{$pair};
		my $s = join "\t", ($label,$region2,$region1,$id2,$id1,$e_id2,$e_id1,$cover{$pair}{$e2},$perc_id{$pair}{$e2});
		print $s."\n";
		push( @found, $e2 );
	    }
	}
	unless ( $adj_list{$e2} && scalar @{$adj_list{$e2}} ){
	    my $s = join "\t", ("MISSING",$region2,$region1,$id2,$id1,$e_id2);
	    print "$s\n";
	    push (@missing, $e2);
	}
	
    }
    #print "total exons: ".scalar(@$exons2)."\n";
    #print "found      : ".scalar(@found)."\n";
    #print "missing    : ".scalar(@missing)."\n";
    

    
    
}


############################################################

sub cluster{
    my ($exons1, $exons2, $score) = @_;


}

############################################################
    
sub blast_Exons{
    my ( $exon1,$exon2, $dir1, $dir2) = @_;

    my $verbose = 0;
    
    my $seqname1 = $exon1->seqname;
    $seqname1 = "chr".$seqname1 unless $seqname1=~/chr/;
    my $chrfile1 = $dir1."/".$seqname1.".fa";

    my $seqname2 = $exon2->seqname;
    $seqname2 = "chr".$seqname2 unless $seqname2=~/chr/;
    my $chrfile2 = $dir2."/".$seqname2.".fa";
    
    # query
    my $id1 = $exon1->transcript_tag.":".$exon1->start."-".$exon1->end.":".$exon1->strand;
    #target
    my $id2 = $exon2->transcript_tag.":".$exon2->start."-".$exon2->end.":".$exon2->strand;
    
    print "comparing $id1 and $id2\n" if $verbose;
    
    my ($seq1, $seq2 );
    my $string1 = Adaptor::SequenceAdaptor->get_exon_seq($chrfile1,$exon1);
    $seq1       = Bio::Seq->new(
				-DISPLAY_ID => $id1,
				-MOLTYPE    => 'dna',
				-SEQ        => $string1,
				);
    
    my $string2 = Adaptor::SequenceAdaptor->get_exon_seq($chrfile2,$exon2);
    $seq2       = Bio::Seq->new(
				-DISPLAY_ID => $id2,
				-MOLTYPE    => 'dna',
				-SEQ        => $string2,
				);
    
    my $length1 = $seq1->length;
    my $length2 = $seq2->length;
    
    ############################################################
    # create database
    my $file = 'seq_'.$$.'.fa';
    my $database = "/tmp/".$file;
    open( DB_SEQ,">$database") || die("Could not open $database $!");
    
    my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
                                 '-fh'     => \*DB_SEQ);
    
    $seqout->write_seq($seq2);
    close( DB_SEQ );
    system("pressdb $database > /dev/null 2>&1");
    
    ############################################################
    # create query file
    my $query_file = '/tmp/query_'.$$.'.fa';
    open ( QUERY, ">$query_file") or die ("cannot open query file $query_file for writing");
    
    my $seqout2 = Bio::SeqIO->new('-format' => 'Fasta',
                                  '-fh'     => \*QUERY);
    $seqout2->write_seq($seq1);
    close(QUERY);
    
    ############################################################
    # generate blast command

    #my $options = "-nogap W=5";
    # Ian's parameters:
    #my $options = "W=5 M=1 N=-1 Q=3 R=3";
    my $options = "W=5 -topcomboN 1";
    
    my $command = "blastn $database $query_file $options";
    
    #open( OUT, "$command | parseblast.pl -GF |" );
    open( OUT, "$command | ");
    
    # ENST00000335466.2_37369656-37371283_1 BLASTN hsp 1 514 992 + . ENSMUST00000046463.2_80047601-80049069_1; Start 1; End 511; Strand +; Frame .; 
    while(<OUT>){
        print $_;
    }
    
    unlink ( $database );
    unlink ($query_file);
}


############################################################

sub exonerate_Exons{
    my ( $exon1,$exon2, $dir1, $dir2) = @_;

    my $verbose = 1;
    
    my $seqname1 = $exon1->seqname;
    $seqname1 = "chr".$seqname1 unless $seqname1=~/chr/;
    my $chrfile1 = $dir1."/".$seqname1.".fa";

    my $seqname2 = $exon2->seqname;
    $seqname2 = "chr".$seqname2 unless $seqname2=~/chr/;
    my $chrfile2 = $dir2."/".$seqname2.".fa";
    
    # query
    my $id1 = $exon1->transcript_tag.":".$exon1->start."-".$exon1->end.".".$exon1->strand;
    #target
    my $id2 = $exon2->transcript_tag.":".$exon2->start."-".$exon2->end.".".$exon2->strand;
    
    print "comparing $id1 and $id2\n" if $verbose;
    
    my ($seq1, $seq2 );
    my $string1 = Adaptor::SequenceAdaptor->get_exon_seq($chrfile1,$exon1);
    $seq1       = Bio::Seq->new(
				-DISPLAY_ID => $id1,
				-MOLTYPE    => 'dna',
				-SEQ        => $string1,
				);
    
    my $string2 = Adaptor::SequenceAdaptor->get_exon_seq($chrfile2,$exon2);
    $seq2       = Bio::Seq->new(
				-DISPLAY_ID => $id2,
				-MOLTYPE    => 'dna',
				-SEQ        => $string2,
				);
    
    my $length1 = $seq1->length;
    my $length2 = $seq2->length;
    
    ############################################################
    # target
    my $file = 'seq_'.$$.'.fa';
    my $target = "/tmp/".$file;
    open( DB_SEQ,">$target") || die("Could not open $target $!");
    
    my $seqout = Bio::SeqIO->new('-format' => 'Fasta',
                                 '-fh'     => \*DB_SEQ);
    
    $seqout->write_seq($seq2);
    close( DB_SEQ );
        
    ############################################################
    # create query file
    my $query_file = '/tmp/query_'.$$.'.fa';
    open ( QUERY, ">$query_file") or die ("cannot open query file $query_file for writing");
    
    my $seqout2 = Bio::SeqIO->new('-format' => 'Fasta',
                                  '-fh'     => \*QUERY);
    $seqout2->write_seq($seq1);
    close(QUERY);
    
    ############################################################
    # generate exonerate command
    
    
    my $options = " --showvulgar FALSE ";
    $options .=   " --showalignment TRUE ";
    $options .=  " --bestn 1 ";
    #$options .=   " -m ungapped:trans ";
    #$options .=  " -m coding2coding ";
    $options .=  " -m affine:local ";
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

    my $cover = int(100*($cover1 + $cover2)/2)/100;
    return ($perc_id,$cover);
}



