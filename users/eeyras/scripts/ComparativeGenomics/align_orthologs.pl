#!/usr/bin/perl

use strict;


# INTRODUCIR TRANSCRIPT_ID (HUMANO O RATON, TIENE QUE DAR LO MISMO)

#print "enter the gene id or the transcript id: \n";
#my $id = <STDIN>;
#chomp ($id);

#print "enter the organism: \n";
#my $especie = <STDIN>;
#chomp ($especie);

#print "enter whether your sequence is a gene or a transcript: \n";
#my $tog = <STDIN>;
#chomp($tog);

# BUSCAR EL GENE_ID CORRESPONDIENTE AL T_ID E IMPRIMIR G_ID (HASH DE G_ID Y T_ID)

if (!open(TRANSCRIPTSGENEID,"< transcripts_to_genes.gff")) {
    print "t_to_g.pl: no se puede abrir el fichero transcripts_to_genes.gff\n";
    exit;
}

my %transcripts;              # hash donde registraremos la matriz transcript-gene
my %genes;                    # hash donde registraremos la matriz gene-transcript
my @ids;                      # vector donde registraremos los ids (de transcript o gene)

#my $tog = $ARGV[0];          # variable introducida indicadora de transcript o gen
#my $id = $ARGV[1];           # variable introducida indicadora de id
    
while (<TRANSCRIPTSGENEID>){
    chomp;
    @ids = split " ", $_;  
    
    if ($ids[0] =~ m/(ENST|ENSMUST)/ && $ids[1] =~ m/(ENSG|ENSMUSG)/ ) {
	push ( @{$transcripts{$ids[0]}}, $ids[1] );
	push ( @{$genes{$ids[1]}}, $ids[0] );
    }
}

#my $gene_id;

#if ( uc($tog) eq 'TRANSCRIPT' ) {
#    my $string = join " ", @{$transcripts{$id}};
#    #print "---> $id es el transcript correspondiente al gen ",$string,"\n";
#    $gene_id = $string;
#}

#if ( uc($tog) eq 'GENE' ){
#    my $string = join " ", @{$genes{$id}};
#    #print "---> $id es el gen correspondiente al transcript ",$string,"\n";
#    $gene_id = $id;
#}

close TRANSCRIPTSGENEID;
    
# BUSCAR EL ORTOLOGO DE G_ID E IMPRIMIR OG_ID
    
if (!open(ORTOLOGOS,"< ortologos.fa")) {
    print "ortologos.pl: no se puede abrir el fichero ortologos.fa\n";
    exit;
}
    
my %Hortologos;              # hash donde registraremos la matriz de ortologos human-mouse
my %Mortologos;              # hash donde registraremos la matriz de ortologos mouse-human
my @genes;                   # vector donde registraremos los gene_id

while (<ORTOLOGOS>){
    chomp;
    @genes = split "\t", $_;  
    
    if ($genes[0] =~ m/ENS/ && $genes[1] =~ m/ENSMUS/ ) {
	#print "genes ortologos: $genes[0] <--> $genes[1]\n";
	push ( @{$Hortologos{$genes[0]}}, $genes[1] );
	push ( @{$Mortologos{$genes[1]}}, $genes[0] );
    }
}



if (!open(HUMANGENES,"< human_uniq.fa")) {
    print "no se puede abrir el fichero inicial\n";
    exit;
}

while (<HUMANGENES>) {
    chomp;
    my $especie = "HUMAN";     # variable introducida indicadora de especie
    my $id = $_;     # variable introducida indicadora del gen
    my $tog = "GENE";
    print "$id \n";
    
    if (!$genes{$id} || !$Hortologos{$id}) {
	print "no existe transcrito u ortologo de $id \n";
	next;
    }

    my $gene_id2;
    my @orto_genes;
    
    if ( uc($especie) eq 'HUMAN' ) {
	$gene_id2 = join " ", @{$Hortologos{$id}};
	@orto_genes = @{$Hortologos{$id}};
	#print "OUTPUT $id \t",$gene_id2,"\n";
    }
    
    if ( uc($especie) eq 'MOUSE' ){
	$gene_id2 = join " ", @{$Mortologos{$id}};
	@orto_genes = @{$Mortologos{$id}};
	#print "OUTPUT $id - ",$gene_id2,"\n";
    }
    
    my $orto_index=0;
    
  ORTO_GENE:
    while ( $orto_index < scalar @orto_genes ){
	my $gene_id2 = $orto_genes[$orto_index];
	my @htrans;
	my @mtrans;
	print "OUTPUT GENE $id - ",$gene_id2,"\n";

	if (!$genes{$gene_id2}) {
	    print "no existe transcrito para el ortologo $gene_id2\n";
	    $orto_index = $orto_index + 1;
	    next;
	}

	if ( uc($especie) eq 'HUMAN' ) {
	    @htrans = @{$genes{$id}};
	    @mtrans = @{$genes{$gene_id2}};
	}
	
	if ( uc($especie) eq 'MOUSE' ){
	    @htrans = @{$genes{$gene_id2}};
	    @mtrans = @{$genes{$id}};
	}

	close(ORTOLOGOS);

	close(TRANSCRIPTSGENEID);

	# ALINEAR (CLUSTALW) T_ID CON OT1_ID 
	
	my $i = 0;
	my $j = 0;
	
	while ($i < scalar(@htrans)){
	    while ($j < scalar(@mtrans)){
		my $hpid = $htrans[$i];
		my $mpid = $mtrans[$j];
				
		my $hgff = get_gff($hpid, "human_transcripts.gff");
		my $mgff = get_gff($mpid, "mouse_transcripts.gff");
		
		
		if ( scalar( @$hgff ) != 1 && scalar( @$mgff ) != 1 ){
		    
		    my $hseq = get_seq($hpid,"/disc8/home/marbla04/human_peptides.fa");
		    my $mseq = get_seq($mpid,"/disc8/home/marbla04/mouse_peptides.fa");
		    
		    open(FA,">pep.fa");
		    print FA ">$hpid\n";
		    print FA @$hseq;
		    print FA ">$mpid\n";
		    print FA @$mseq;
		    close(FA);
		    
		    system("/disc8/bin/clustalw1.8/clustalw pep.fa");
		    
		    # clustalw produce pep.aln
		     
		    # APLICAR EXSTRAL.PL (SI NO SALE, VOLVER AL PASO ANTERIOR Y HACERLO CON OT2_ID...)
		    
		    open(HGFF,"> hum.gff");
		    print HGFF @$hgff;
		    close(HGFF);
		    
		    open(MGFF,"> mus.gff");
		    print MGFF @$mgff;
		    close(MGFF);
			    
		    open(EX, "/disc8/home/marbla04/exstral.pl  hum.gff mus.gff pep.aln | "); 
		    
		    while(<EX>){
			print "-->".$_;	
			chomp;
			# # ENST00000005593.1  ENSMUST00000016463.1    3   3   3   0.98 1.00  1.00
			next if (!/\#/);
			my @result = split '\s+', $_;
			#print "[1] = $result[1]\n";
			#print "[2] = $result[2]\n";
			#print "[6] = $result[6]\n";
			#print "[7] = $result[7]\n";
			#print "[8] = $result[8]\n";
			if ($result[8] == 1){
			    print "OUTPUT TRAN $result[1] - $result[2] \t GENE $id - ",$gene_id2,"\n";
			}
		    }
		}
		else{
		    if ( scalar( @$hgff ) == 1 ){
			print "$hpid tiene un solo exon \n"}
		    if ( scalar( @$mgff ) == 1 ){
			print "$mpid tiene un solo exon \n"}
		}
		$j = $j + 1;
	    }
	    $i = $i + 1;
	}
	$orto_index = $orto_index + 1;
    }  # fin del bucle ORTO_GENE 
} #fin del megabucle
    
sub get_gff{
    my $t_id = $_[0];
    my $file = $_[1];
    my $gff;
    my @lines;
    
    open(GFF,"<$file"); 
    
    while(<GFF>){
	
	if ( /$t_id/ && /CDS/ ){
	    my @entries = split;
	    my $linea = $t_id."\tensembl\tCDS\t".
		$entries[3]."\t".$entries[4]."\t.\t".
		$entries[6]."\t".$entries[7]."\t".$t_id."\n";
	    push(@lines, $linea);
	}
    }
    
    close(GFF);
    
    $gff = \@lines;
    return $gff;
}



sub get_seq{
    my $t_id = $_[0];
    my $file = $_[1];
    my $sequence;
    my @lines;
    
    #print "reading $file\n";
    #print "looking for $t_id\n";
    
    open(PEP,"<$file"); 
    
  FILE:
    while(<PEP>){
	
	if ( /\>/ && /$t_id/ ){
	    #print $_;
	    
	    while(<PEP>){
		
		if ( !/\>/ ){
		    #print $_;
		    push( @lines, $_ );
		}
		else{
		    last FILE;
		}
	    }
	}
    }
    $sequence = \@lines;
    
    close(PEP);
    
    return $sequence;
}
