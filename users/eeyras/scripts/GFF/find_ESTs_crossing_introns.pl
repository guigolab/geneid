#!/usr/local/bin/perl -w

use strict;
use ClusterMerge::Exon;
use ClusterMerge::ExonUtils;
use ClusterMerge::GFFTools;

# this program reads two files. One with the GFF with the gene structures of target genes (e.g. annotations or predictions)
# and one GFF with the blast results of comparing ESTs with these genes.
# The aim is to find out which ESTs cross an intron in the target genes, like this:
#
#
#           EST             ##----###
#   target gene   ####----####----####----####
#


my $verbose = 0;

if (scalar(@ARGV) < 2) {
    print STDERR "$0 target.genes.gff  blast.features.gff\n";
    print STDERR "\n";
    print STDERR "This program reads two files. One with the GFF with the gene structures of target genes\n";
    print STDERR "and one GFF with the blast results of comparing ESTs with these genes\n";
    print STDERR "and finds out which ESTs cross an intron in the target genes\n";
    exit(0);
}

my $target_file = $ARGV[0];
my $blast_file  = $ARGV[1]; 

############################################################


my @targets;
my %trans;
open (ANN,"<$target_file") or die("cannot open $target_file");
while (<ANN>){
    chomp;
    my @e = split;
    next unless ( $e[3] && $e[4] ); 
    my $exon = ClusterMerge::GFFTools->exon_from_gff( $_);
    $exon->source_tag("target");
    if ( $exon && $exon->start && $exon->end ){
	push( @{$trans{$exon->group_tag}}, $exon );
    }
}

############################################################
# hold the exon number for each position within the transcript:
my %exon_position;

foreach my $tag ( keys %trans ){
    my $transcript = ClusterMerge::Transcript->new();
    $transcript->dbID( $tag );
    my $seqname;
    my $strand;
    my @exons = sort { $a->start <=> $b->start } @{ $trans{$tag} };
    my $start = 1;
    my $exon_count = 1;
    foreach my $exon ( @exons ){
	print "exon: starts=$start - ends=".($start + $exon->length - 1)." length= ".$exon->length."\n" if $verbose;
	
	for (my $i=$start; $i<= ($start + $exon->length - 1) ; $i++){
	    $exon_position{$exon->group_tag}->[$i] = $exon_count;
	    #print $exon->group_tag." pos=$i exon=$exon_count\n";
	}
	$exon_count++;
	$start += $exon->length;
    }
}
close(ANN);

############################################################

# the blast input is like this: ( from parseblast.pl -F -s )

#chr15_51	BLASTN	hsp	1	233	161	+	.	Target "gb_BG106491.1_BG106491";	Start 401;	End 629;	Strand +;	Frame .;	E_value 4e-84;	P_sum .0;	Aln_Sore 161;	Bit_Score 319;	Idn_Score  95.63 (219/229);	Gaps 1.72 (Q:0|S:4);	Lengths Q:233|S:229|T:233
#chr15_51	BLASTN	hsp	1	190	158	+	.	Target "gb_BU171451.1_BU171451";	Start 386;	End 575;	Strand +;	Frame .;	E_value 3e-82;	P_sum .0;	Aln_Sore 158;	Bit_Score 313;	Idn_Score  95.79 (182/190);	Gaps 0.00 (Q:0|S:0);	Lengths Q:190|S:190|T:190


open (BLAST,"<$blast_file") or die("cannot open $blast_file");

my %hits;
while (<BLAST>){
    chomp;
    my @e = split;
    my $target = $e[0];
    my $start  = $e[3];
    my $end    = $e[4];

    my $where_starts = $exon_position{$target}->[$start];
    my $where_ends   = $exon_position{$target}->[$end];
    if ( $where_starts != $where_ends ){
	print $_."\n";
    }
}
close(BLAST);

############################################################
