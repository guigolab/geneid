#!/usr/local/bin/perl -w

use strict;

use ClusterMerge::GTFTools;
use GeneComparison::GenePair;
use GeneComparison::Alineator;
use GeneComparison::TranscriptComparison;


my $human_id   = $ARGV[0];
my $human_file = $ARGV[1];
my $mouse_id   = $ARGV[2];
my $mouse_file = $ARGV[3];


my $dir1 = "/seq/genomes/H.sapiens/golden_path_200307/chromFa/";
my $dir2 = "/seq/genomes/M.musculus/golden_path_200405mm5/chromFa/";


unless ( $human_id && $human_file && $mouse_id && $mouse_file ){
    print STDERR "Usage: $0 human.id human.gtf mouse.id mouse.gtf\n";
    exit(0);
}


my ($human_genes,$human_gene_hash) = ClusterMerge::GTFTools->get_genes_from_GTF($human_file);
my ($mouse_genes,$mouse_gene_hash) = ClusterMerge::GTFTools->get_genes_from_GTF($mouse_file);


my $gene1 = $human_gene_hash->{$human_id};
my $gene2 = $mouse_gene_hash->{$mouse_id};
    

compare($gene1,$gene2,$dir1,$dir2,1);



sub compare{
    my ($human_gene,$mouse_gene, $dir1, $dir2, $coding_exons) = @_;
    
    my @human_transcripts = @{$human_gene->get_Transcripts};
    my @mouse_transcripts = @{$mouse_gene->get_Transcripts};
    
    # for each mouse gene, see whether you find them all in human:
    my @transcript_matches;
    foreach my $mouse_t ( @mouse_transcripts ){
	foreach my $human_t ( @human_transcripts ){
	    
	    my $m_t = $mouse_t->coding_transcript;
	    my $h_t = $human_t->coding_transcript;
	    
	    
	    GeneComparison::TranscriptComparison->compare_fingerprints( $human_t, $mouse_t );
	    
	}
    }
    
}
