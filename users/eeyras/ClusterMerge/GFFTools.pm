package ClusterMerge::GFFTools;

use vars qw(@ISA);
use strict;

use ClusterMerge::Root;
use ClusterMerge::Transcript;
use ClusterMerge::Exon;

@ISA = qw(ClusterMerge::Root);


############################################################
# this method gets all the genes from a GTF file
# the GTF format contains the genes and transcript
# structures, hence transcripts are put into ClusterMerge::Transcript objects
# and genes are created as ClusterMerge::TranscriptCluster objects, i.e. they are
# a set of transcripts. You get the trasncripts from the gene like this:
#   if $gene is a ClusterMerge::TranscriptCluster object
#   my @transcripts = @{$gene->get_Transcripts};
#
#   the gene identifier is stored as dbID:
#   my $gene_id = $gene->dbID;
# see the ClusterMerge::TranscriptCluster for more details

sub get_genes_from_GFF{
    my ($self,$file) = @_;

    my @transcripts = $self->get_transcripts_from_GFF($file);
    
    my @genes = ClusterMerge::TranscriptUtils->(\@transcripts,"mix");
    # these genes are ClusterMerge::TranscriptCluster objects
    return @genes;
}

sub get_transcripts_from_GFF{
    my ($self,$file) = @_;
    
    my $verbose = 0;

    my $coding_label = "coding_exon";
    my $exon_label   = "exon";
  
    open ( IN, "<$file") || $self->throw("could not open input file $file");
    my @transcripts;
    
    my %tag2transcript;
    my %tag2codingtranscript;
    my %transcript2codingtranscript;
    
    while(<IN>){
	my @e = split;
	next unless ( $e[3] && $e[4] ); 
	
	if ( $e[2] eq $exon_label ){
	    my $exon = $self->exon_from_gff( $_ );
	    
	    # create a transcript
	    unless ( $tag2transcript{$exon->transcript_tag} ){
		my $trans_id = $exon->transcript_tag;
		$tag2transcript{$trans_id} =  ClusterMerge::Transcript->new();
		$tag2transcript{$trans_id}->dbID( $trans_id );
		$tag2transcript{$trans_id}->type( $exon->source_tag );
		push( @transcripts, $tag2transcript{$trans_id} );
	    }
	    $tag2transcript{$exon->transcript_tag}->add_Exon($exon);
	    
	}
	elsif( $e[2] eq $coding_label ){
	    my $coding_exon = $self->exon_from_gff( $_);
	    	    
	    print "seen a coding transcript: $_\n" if $verbose;
	    # create a coding transcript
	    unless ( $tag2codingtranscript{$coding_exon->group_tag} ){
		my $trans_id = $coding_exon->transcript_tag;
		$tag2codingtranscript{$trans_id} =  ClusterMerge::Transcript->new();
		$tag2codingtranscript{$trans_id}->dbID( $trans_id );
		$tag2codingtranscript{$trans_id}->type( $coding_exon->source_tag );
	    }
	    unless ( $transcript2codingtranscript{$tag2transcript{$coding_exon->transcript_tag}} ){
		$transcript2codingtranscript{$tag2transcript{$coding_exon->transcript_tag}} = $tag2codingtranscript{$coding_exon->transcript_tag};
		$tag2transcript{$coding_exon->group_tag}->coding_transcript($tag2codingtranscript{$coding_exon->transcript_tag});
	    }
	    $tag2codingtranscript{$coding_exon->transcript_tag}->add_Exon($coding_exon);
	}
    }
    close(IN);
    if ($verbose){
	foreach my $t (@transcripts){
	    my $cod = "NO";
	    if ($t->coding_transcript){
		$cod = "YES";
	    }
	    print "found transcript: ".$t->dbID." coding $cod\n";
	}
    }

    return @transcripts;
}



############################################################
# read annotations

sub get_simple_transcripts_from_GFF{
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

1;

