#!/usr/local/bin/perl -w

use strict;

use ClusterMerge::GTFTools;
use GeneComparison::TranscriptComparison;

my $human_id = $ARGV[0];
my $mouse_id = $ARGV[1];

unless ( $human_id && $mouse_id ){
    print STDERR "Usage: $0 human_id mouse_id\n";
    exit(0);
}


my $human_file = "/home/ug/eeyras/Projects/ATD/ComparativeGenomics/HumanMouse/AllExons/noNT.human.gtf";
my $mouse_file = "/home/ug/eeyras/Projects/ATD/ComparativeGenomics/HumanMouse/AllExons/noNT.mouse.gtf";

my ($human_genes,$human_gene_hash) = ClusterMerge::GTFTools->get_gene_from_GTF($human_file, $human_id);
my ($mouse_genes,$mouse_gene_hash) = ClusterMerge::GTFTools->get_gene_from_GTF($mouse_file, $mouse_id);

# get the orthologous pairs:

my $object_map = ClusterMerge::ObjectMap->new();

my %t2g;
my %alignment;
my %fingerprint;

my $gene1 = $human_gene_hash->{$human_id};
my $gene2 = $mouse_gene_hash->{$mouse_id};

unless ($gene1 && $gene2){
    print "gene $human_id not found" unless ($gene1);
    print "gene $mouse_id not found" unless ($gene2);
    exit(0);
}

foreach my $t1 ( @{$gene1->get_Transcripts} ){
    
    $t2g{$t1} = $gene1->dbID;
    foreach my $t2 ( @{$gene2->get_Transcripts} ){
	
	$t2g{$t2} = $gene2->dbID;
	
	print "here\n";
	my ($score,$finger1,$finger2,$alignment1,$alignment2,$gaps1,$gaps2) = GeneComparison::TranscriptComparison->compare_fingerprints($t1,$t2);
	
	my $s = int(100*$score)/100;
	$fingerprint{$t1} = $finger1;
	$fingerprint{$t2} = $finger2;
	
	$alignment{$t1} = $alignment1;
	$alignment{$t2} = $alignment2;
	$object_map->match($t1, $t2, $s );
    }
}

my $best_pairs_object = $object_map->stable_marriage;

############################################################
# pairs created:
my $pair_count = scalar($best_pairs_object->list1);
print "pairs created: ".$pair_count."\n";
foreach my $element1 ( $best_pairs_object->list1 ){
    foreach my $partner ( $best_pairs_object->partners( $element1 ) ){
	# there should be only one
	my $id1 = $element1->dbID;
	my $id2 = $partner->dbID;
	print "PAIR\t$t2g{$element1}\t$id1\t$t2g{$partner}\t$id2\t".$best_pairs_object->score( $element1, $partner )."\n";
	print "FING\t".fingerprint_string($fingerprint{$element1})."\n";
	print "FING\t".fingerprint_string($fingerprint{$partner} )."\n";
	print "ALIG\t@{$alignment{$element1}}\n";
	print "ALIG\t@{$alignment{$partner}}\n";
    }
}



sub fingerprint_string{
    my ($f) = @_;
    my $s = join "\t", @$f;
    return $s;
}
