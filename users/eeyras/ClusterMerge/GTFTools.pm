package ClusterMerge::GTFTools;

use vars qw(@ISA);
use strict;

use ClusterMerge::Root;
use ClusterMerge::Transcript;
use ClusterMerge::TranscriptCluster;
use ClusterMerge::Exon;

@ISA = qw(ClusterMerge::Root);



############################################################
# example:
#Y       EnsEMBL start_codon     8899746 8899748 .       +       .       gene_id "ENSG00000182061.1"; transcript_id "ENST00000316670.2"; exon_id "ENSE00001236096.2";
#Y       EnsEMBL exon    8899920 8900077 .       +       .       gene_id "ENSG00000182061.1"; transcript_id "ENST00000316670.2"; exon_id "ENSE00000981599.2";
#Y       EnsEMBL CDS     8899920 8900077 .       +       2       gene_id "ENSG00000182061.1"; transcript_id "ENST00000316670.2"; exon_id "ENSE00000981599.2";
#Y       EnsEMBL exon    8915926 8916025 .       -       .       gene_id "ENSG00000183652.1"; transcript_id "ENST00000332464.1"; exon_id "ENSE00001298905.1";
#Y       EnsEMBL CDS     8915926 8916025 .       -       0       gene_id "ENSG00000183652.1"; transcript_id "ENST00000332464.1"; exon_id "ENSE00001298905.1";
#Y       EnsEMBL exon    8915694 8915839 .       -       .       gene_id "ENSG00000183652.1"; transcript_id "ENST00000332464.1"; exon_id "ENSE00001315537.1";
#Y       EnsEMBL CDS     8915694 8915839 .       -       2       gene_id "ENSG00000183652.1"; transcript_id "ENST00000332464.1"; exon_id "ENSE00001315537.1";
#Y       EnsEMBL exon    8914091 8914189 .       -       .       gene_id "ENSG00000183652.1"; transcript_id "ENST00000332464.1"; exon_id "ENSE00001327631.2";
#Y       EnsEMBL CDS     8914094 8914189 .       -       0       gene_id "ENSG00000183652.1"; transcript_id "ENST00000332464.1"; exon_id "ENSE00001327631.2";
#Y       EnsEMBL stop_codon      8914091 8914093 .       -       .       gene_id "ENSG00000183652.1"; transcript_id "ENST00000332464.1"; exon_id "ENSE00001327631.2";

sub get_transcripts_from_GTF{
    my ($self,$file) = @_;
    
    open ( IN, "<$file") || $self->throw("could not open input file $file");
    my @transcripts;
    
    my %tag2transcript;
    my %tag2codingtranscript;
    my %transcript2codingtranscript;
    
    while(<IN>){
	my @e = split;
	next unless ( $e[3] && $e[4] ); 
	
	if ( $e[2] =~/exon/ ){
	    my $exon = $self->exon_from_gtf( $_ );

	    # create a transcript
	    unless ( $tag2transcript{$exon->transcript_tag} ){
		my $trans_id = $exon->transcript_tag;
		$tag2transcript{$trans_id} =  ClusterMerge::Transcript->new();
		$tag2transcript{$trans_id}->dbID( $$trans_id );
		push( @transcripts, $tag2transcript{$trans_id} );
	    }
	    $tag2transcript{$exon->transcript_tag}->add_Exon($exon);
	    
	}
	elsif( $e[2] =~/CDS/ ){
	    my $coding_exon = $self->exon_from_gtf( $_);

	    # create a coding transcript
	    unless ( $tag2codingtranscript{$coding_exon->group_tag} ){
		my $trans_id = $coding_exon->transcript_tag;
		$tag2codingtranscript{$trans_id} =  ClusterMerge::Transcript->new();
		$tag2codingtranscript{$trans_id}->dbID( $coding_exon->transcript_tag );
	    }
	    unless ( $transcript2codingtranscript{$tag2transcript{$coding_exon->transcript_tag}} ){
		$transcript2codingtranscript{$tag2transcript{$coding_exon->transcript_tag}} = $tag2codingtranscript{$coding_exon->transcript_tag};
		$tag2transcript{$coding_exon->group_tag}->coding_transcript($tag2codingtranscript{$coding_exon->transcript_tag});
	    }
	    $tag2codingtranscript{$coding_exon->transcript_tag}->add_Exon($coding_exon);
	}
    }
    close(IN);
    return @transcripts;
}

############################################################
# this method is like the "get_genes_from_GTF" method
# except that it retrieves one specific gene_id
#
# this is handy when we want to run different queries on different
# machines.
#
# The idea of the method is very simple: it does a grep on the GTF file
# with the requested gene_id, and then it reads the gene/transcript/exon
# structure as usual piping the output of the grep.
# We use pipe rather than writing to a file and calling the 
# "get_genes_from_GTF" method as it seems to be faster.

sub get_gene_from_GTF{
    my ($self,$file,$gene_id) = @_;

    my $verbose = 0;

    my @transcripts;
    my @genes;

    my %tag2gene;
    my %tag2transcript;
    my %tag2codingtranscript;
    my %transcript2codingtranscript;

    open (IN, "egrep $gene_id $file | ");
    while(<IN>){
	my @e = split;
	next unless ( $e[3] && $e[4] ); 
	
	if ( $e[2] =~/exon/ ){
	    my $exon = $self->exon_from_gtf( $_ );

	        # create a transcript
	    unless ( $tag2transcript{$exon->transcript_tag} ){
		my $trans_id = $exon->transcript_tag;
		$tag2transcript{$trans_id} =  ClusterMerge::Transcript->new();
		$tag2transcript{$trans_id}->dbID( $trans_id );
		push( @transcripts, $tag2transcript{$trans_id} );
		
		# create a gene
		push (@{$tag2gene{$exon->gene_tag}}, $tag2transcript{$trans_id});
	    }
	    $tag2transcript{$exon->transcript_tag}->add_Exon($exon);
	}
	elsif( $e[2] =~/CDS/ ){
	    my $coding_exon = $self->exon_from_gtf( $_);

	        #print "GTF: detected coding exon\n";
	        # create a coding transcript
	    unless ( $tag2codingtranscript{$coding_exon->group_tag} ){
		my $trans_id = $coding_exon->transcript_tag;
		$tag2codingtranscript{$trans_id} =  ClusterMerge::Transcript->new();
		$tag2codingtranscript{$trans_id}->dbID( $trans_id );
		#print "GTF: created coding transcript\n";
	    }
	    unless ( $transcript2codingtranscript{$tag2transcript{$coding_exon->transcript_tag}} ){
		$transcript2codingtranscript{$tag2transcript{$coding_exon->transcript_tag}} = $tag2codingtranscript{$coding_exon->transcript_tag};
		$tag2transcript{$coding_exon->group_tag}->coding_transcript($tag2codingtranscript{$coding_exon->transcript_tag});
		#print "GTF: attached coding transcript\n";
	    }
	    $tag2codingtranscript{$coding_exon->transcript_tag}->add_Exon($coding_exon);
	}
    }
    close(IN);
    
    ############################################################
    # create genes (we do this here and not in the fly
    # as everytime we create a TranscriptCluster we reset the boundaries 
    # and here we cannot do it without exons
    my %gene_hash;
    foreach my $id ( keys %tag2gene ){
	my @t = @{$tag2gene{$id}};
	delete $tag2gene{$id};
	my $gene = ClusterMerge::TranscriptCluster->new();
	$gene->dbID( $id );
	$gene->put_Transcripts(@t);
	$gene_hash{$id} = $gene;
	push( @genes, $gene );
    }
    if ($verbose){
	print "found genes:\n";
	foreach my $g (@genes){
	    $g->print_Cluster;
	}
    }

    return (\@genes, \%gene_hash); # ClusterMerge::TranscriptCluster objects
}

############################################################
# this method gets all the genes from a GTF file
# teh GTF format contains the genes and transcript
# structures, hence transcripts are put into ClusterMerge::Transcript objects
# and genes are created as ClusterMerge::TranscriptCluster objects, i.e. they are
# a set of transcripts. You get the trasncripts from the gene like this:
#   if $gene is a ClusterMerge::TranscriptCluster object
#   my @transcripts = @{$gene->get_Transcripts};
#
#   the gene identifier is stored as dbID:
#   my $gene_id = $gene->dbID;
# see the ClusterMerge::TranscriptCluster for more details

sub get_genes_from_GTF{
    my ($self,$file) = @_;

    my $verbose = 0;
    
    open ( IN, "<$file") || $self->throw("could not open input file $file");
    my @transcripts;
    my @genes;

    my %tag2gene;
    my %tag2transcript;
    my %tag2codingtranscript;
    my %transcript2codingtranscript;

    while(<IN>){
	my @e = split;
	next unless ( $e[3] && $e[4] ); 
	
	if ( $e[2] =~/exon/ ){
	    my $exon = $self->exon_from_gtf( $_ );

	    # create a transcript
	    unless ( $tag2transcript{$exon->transcript_tag} ){
		my $trans_id = $exon->transcript_tag;
		$tag2transcript{$trans_id} =  ClusterMerge::Transcript->new();
		$tag2transcript{$trans_id}->dbID( $trans_id );
		push( @transcripts, $tag2transcript{$trans_id} );
		
		# create a gene
		push (@{$tag2gene{$exon->gene_tag}}, $tag2transcript{$trans_id});
	    }
	    $tag2transcript{$exon->transcript_tag}->add_Exon($exon);
	}
	elsif( $e[2] =~/CDS/ ){
	    my $coding_exon = $self->exon_from_gtf( $_);

	    #print "GTF: detected coding exon\n";
	    # create a coding transcript
	    unless ( $tag2codingtranscript{$coding_exon->group_tag} ){
		my $trans_id = $coding_exon->transcript_tag;
		$tag2codingtranscript{$trans_id} =  ClusterMerge::Transcript->new();
		$tag2codingtranscript{$trans_id}->dbID( $trans_id );
		#print "GTF: created coding transcript\n";
	    }
	    unless ( $transcript2codingtranscript{$tag2transcript{$coding_exon->transcript_tag}} ){
		$transcript2codingtranscript{$tag2transcript{$coding_exon->transcript_tag}} = $tag2codingtranscript{$coding_exon->transcript_tag};
		$tag2transcript{$coding_exon->group_tag}->coding_transcript($tag2codingtranscript{$coding_exon->transcript_tag});
		#print "GTF: attached coding transcript\n";
	    }
	    $tag2codingtranscript{$coding_exon->transcript_tag}->add_Exon($coding_exon);
	}
    }
    close(IN);
    
    ############################################################
    # create genes (we do this here and not in the fly
    # as everytime we create a TranscriptCluster we reset the boundaries 
    # and here we cannot do it without exons
    my %gene_hash;
    foreach my $id ( keys %tag2gene ){
	my @t = @{$tag2gene{$id}};
	delete $tag2gene{$id};
	my $gene = ClusterMerge::TranscriptCluster->new();
	$gene->dbID( $id );
	$gene->put_Transcripts(@t);
	$gene_hash{$id} = $gene;
	push( @genes, $gene );
    }
    if ($verbose){
	print "found genes:\n";
	foreach my $g (@genes){
	    $g->print_Cluster;
	}
    }

    return (\@genes, \%gene_hash); # ClusterMerge::TranscriptCluster objects
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

