#!/usr/local/bin/perl -w

use strict;

use ClusterMerge::GFFTools;
use Adaptor::SequenceAdaptor;
use GeneComparison::IntronComparison;
use Bio::Seq;
use Bio::SeqIO;


# human encode region
my $encode_region = $ARGV[0];
unless ( $encode_region ){
    print STDERR "Usage: $0 encode.region\n";
    exit(0);
}

# where the encode sequences are:
my $human_seq_dir = "/projects/encode/data/hg16/";
my $human_seq_file = $human_seq_dir."/".$encode_region."/".$encode_region.".fa";

# where the VEGA annotations are
my $human_vega_dir = "/projects/encode/vega/gff/";



############################################################
# read vega genes for these region
my $vega_file = $human_vega_dir."/".$encode_region."_transl_VEGA_all_exons.gff";

#/projects/encode/data/hg16/ENr333/gff

my @transcripts;
if ( -e $vega_file && `ls $vega_file` ){
    @transcripts = ClusterMerge::GFFTools->get_transcripts_from_GFF($vega_file);
}
else{
    print "$encode_region is not a vega region\n";
    exit(0);
}

# we would like to label each exon 
# whether it is coding/utr/half-coding

my ($introns,$flag) = GeneComparison::IntronComparison->get_uniq_introns_with_coding_flag(@transcripts);

my $cut_point;

my $seq_file = $encode_region.".vega.introns.fa";
my $gff_file = $encode_region.".vega.introns.gff";

open ( OUT, ">$seq_file") or die ("cannot open query file $seq_file");
my $query_seqout = Bio::SeqIO->new('-format' => 'Fasta',
				   '-fh'     => \*OUT);

open ( GFF_OUT, ">$gff_file") or die ("cannot open query file $gff_file");

foreach my $intron (@$introns){
    
    my $tran_string = Adaptor::SequenceAdaptor->get_transcript_seq($human_seq_file,$intron);
    next unless $tran_string;  

    my @exons = sort {$a->start <=> $b->start} @{$intron->get_all_Exons};
    foreach my $exon (@exons){
	print GFF_OUT $exon->gff_string."\t".$flag->{$exon}."\n";
    }
    
    my $tran_seq    = Bio::Seq->new(
				    -DISPLAY_ID => $intron->dbID,
				    -MOLTYPE    => 'dna',
				    -SEQ        => $tran_string,
				    );
    
    
    $query_seqout->write_seq($tran_seq);
}
close(GFF_OUT);
close(OUT);

