package ClusterMerge::GFFTools;

use vars qw(@ISA);
use strict;

use ClusterMerge::Root;
use ClusterMerge::Transcript;
use ClusterMerge::Exon;

@ISA = qw(ClusterMerge::Root);




############################################################
# read annotations

sub get_transcripts_from_GFF{
    my ($self,$file) = @_;
    
    open ( IN, "<$file") || $self->throw("could not open input file $file");
    my @transcripts;
    
    my %trans;
    
    while(<IN>){
	chomp;
	my @e = split;
	next unless ( $e[3] && $e[4] ); 
	my $exon = $self->exon_from_gff( $_ );
	
	if ( $exon && $exon->start && $exon->end ){
	    push( @{$trans{$exon->group_tag}}, $exon );
	}
    }
    foreach my $tag ( keys %trans ){
	my $transcript = ClusterMerge::Transcript->new();
	$transcript->dbID( $tag );
	my $strand;
	foreach my $exon ( @{ $trans{$tag} } ){
	    $transcript->add_Exon($exon);
	}
	push( @transcripts, $transcript );
    }
    close(IN);
    return @transcripts;
}

############################################################

sub exon_from_gff {
    my ($self,$gffstring) = @_;
    
    #print STDERR "gff_string: $gffstring\n";
    
    #chr6 	human_cdna	exon	1030942	1031591	100	+	0	6604
    #chr6 	human_cdna	exon	1034227	1034285	100	+	0	6604
    #chr6 	human_cdna	exon	1035280	1035339	100	+	0	6604
    my ($seqname, 
	$source, 
	$primary, 
	$start, 
	$end, 
	$score, 
	$strand, 
	$frame, 
	@group) 
	= split(/\s+/, $gffstring);
    
    my $exon = ClusterMerge::Exon->new();
    $exon->seqname($seqname);
    $exon->source_tag($source);
    $exon->start($start);
    $exon->end($end);
    $exon->primary_tag($primary);

    #$frame = "." unless( $frame =~ /^\d+$/);
    
    if ( $frame =~ /^\d+$/){
	$exon->phase($frame);
    }
    else{
	$exon->phase(0);
    }
    $exon->frame($exon->phase);
    
    #my $phase = ( 3 - $frame )%3;
    #$exon->phase($phase);
    #$exon->end_phase( ( $exon->phase + $exon->length)%3 );
    
    if ($score eq '.'){
	$exon->score(0);
    }
    elsif ( defined($score) ){
	$exon->score( $score );
    }
    else{
	$exon->score(0);
    }
    
    if ( $strand eq '-' ) { $exon->strand(-1); }
    elsif ( $strand eq '+' ) { $exon->strand(1); }
    elsif ( $strand eq '.' ) { $exon->strand(0); }
    
    ############################################################
    # warning: it parses only the first element of the group
    $exon->transcript_tag($group[0]);
    
    #print "group-tag = $group[0]\n";
    return $exon;
}

############################################################

sub exon_from_gtf {
    my ($self,$gtfstring) = @_;
    
    #print STDERR "gff_string: $gffstring\n";
    
    #1 cdna CDS 475556 475670 . + 0 gene_id "853"; transcript_id "1048"; exon_number "14"; 
    #1 cdna CDS 480393 480532 . + 0 gene_id "853"; transcript_id "1048"; exon_number "15"; 
    #1 cdna CDS 482243 482345 . + 0 gene_id "853"; transcript_id "1048"; exon_number "16"; 
    my ($seqname, 
	$source, 
	$primary, 
	$start, 
	$end, 
	$score, 
	$strand, 
	$frame, 
	$gene_field,
	$gene_tag,
	$transcript_field,
	$transcript_tag,
	$exon_field,
	$exon_tag)
	= split(/\s+/, $gtfstring);
    
    my $exon = ClusterMerge::Exon->new();
    $exon->seqname($seqname);
    $exon->source_tag($source);
    $exon->start($start);
    $exon->end($end);
    $exon->primary_tag($primary);

    #$frame = "." unless( $frame =~ /^\d+$/);
    
    if ( $frame =~ /^\d+$/){
	$exon->phase($frame);
    }
    else{
	$exon->phase(0);
    }
    $exon->frame($exon->phase);
    
    #my $phase = ( 3 - $frame )%3;
    #$exon->phase($phase);
    #$exon->end_phase( ( $exon->phase + $exon->length)%3 );
    
    if ($score eq '.'){
	$exon->score(0);
    }
    elsif ( defined($score) ){
	$exon->score( $score );
    }
    else{
	$exon->score(0);
    }
    
    if ( $strand eq '-' ) { $exon->strand(-1); }
    elsif ( $strand eq '+' ) { $exon->strand(1); }
    elsif ( $strand eq '.' ) { $exon->strand(0); }
    
    $gene_tag =~/\"(.*)\"\;/;
    $exon->gene_tag($1);

    $transcript_tag =~/\"(.*)\"\;/;
    $exon->transcript_tag($1);
    
    $exon_tag=~/\"(.*)\"\;/;
    $exon->exon_tag($1);

    return $exon;
}

############################################################

1;

