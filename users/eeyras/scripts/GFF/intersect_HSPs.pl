#!/usr/local/bin/perl -w

use strict;
use ClusterMerge::Exon;
use ClusterMerge::ExonUtils;
# this program reads in one shot all the annotations in gff format contained in the
# file given as a first argument, and intersects them with the annotations given in another GFF format file
# as a second parameter. the intersection is given by the standard output STDOUT

if (scalar(@ARGV) < 2) {
    print "$0 hsps1.gff hsps2.gff\n";
    print "Description: script to obtain the intersections of exons/HSPs purely on coordinate and strand level\n";
    print "             It reports the HSPs from hsps1.gff that intersects with HSPs in hsps2.gff\n";
    exit(1);
}

my $pred_file = $ARGV[0];
my $ann_file = $ARGV[1]; 

############################################################

my %forward_hsps;
my %reverse_hsps;

my %tag;
my $ann_type;
open (ANN,"<$ann_file") or die("cannot open $ann_file");
while (<ANN>){
    #chr1	SGP_v1.0	Internal	32133	32262	31.70	-	0	chr1_1
    
    chomp;
    my @e = split;      
    next unless ($e[3] && $e[4]);
    
    my $hsp = exon_from_gff ($_,\%tag);

    if( $e[6] eq '+') {
	push ( @{$forward_hsps{$e[0]}}, $hsp );
    } 
    elsif ( $e[6] eq '-' ){
	push ( @{$reverse_hsps{$e[0]}}, $hsp );
    }
    else{
	push ( @{$forward_hsps{$e[0]}}, $hsp );
	push ( @{$reverse_hsps{$e[0]}}, $hsp );
    }
    unless ($ann_type){
	$ann_type = $e[1];
    }
}
close(ANN);

############################################################

my $pred_type;
open (PRED,"<$pred_file") or die("cannot open $pred_file");
while (<PRED>){
    chomp;
    my @e = split;      
    next unless ($e[3] && $e[4]);
    my $hsp = exon_from_gff ($_,\%tag);
    if( $e[6] eq '+') {
	push ( @{$reverse_hsps{$e[0]}}, $hsp );
    } 
    elsif ( $e[6] eq '-' ){
	push ( @{$forward_hsps{$e[0]}}, $hsp );
    }
    else{
	push ( @{$reverse_hsps{$e[0]}}, $hsp );
    	push ( @{$forward_hsps{$e[0]}}, $hsp );
    }
    unless ( $pred_type ){
	$pred_type = $e[1];
    }
}
close(PRED);

############################################################

#print "anntype = $ann_type predtype = $pred_type\n";

foreach my $tid ( keys %forward_hsps ){
    
    my ($clusters,$exon2cluster)= ClusterMerge::ExonUtils->_cluster_Exons(@{$forward_hsps{$tid}} );
    
    foreach my $c ( @$clusters ){
	my @ann  = grep { $_->source_tag eq $ann_type }  @{$c->get_Exons};
	my @pred = grep { $_->source_tag eq $pred_type } @{$c->get_Exons};
	
	if ( @ann && @pred ){
	    foreach my $hsp (@pred){
		print gff_string($hsp,$tag{$hsp})."\n";
	    }
	}
    }
    delete $forward_hsps{$tid};
}

foreach my $tid ( keys %reverse_hsps ){
    
    my ($clusters,$exon2cluster)= ClusterMerge::ExonUtils->_cluster_Exons(@{$reverse_hsps{$tid}} );
    
    foreach my $c ( @$clusters ){
	my @ann  = grep { $_->source_tag eq $ann_type }  @{$c->get_Exons};
	my @pred = grep { $_->source_tag eq $pred_type } @{$c->get_Exons};
	
	if ( @ann && @pred ){
	    foreach my $hsp (@pred){
		print gff_string($hsp,$tag{$hsp})."\n";
	    }
	}
    }
    delete $reverse_hsps{$tid};
}

############################################################

sub gff_string{
  my ($exon, $trans_id) = @_;
  my ($str,$source,$primary_tag,$score,$frame,$name,$strand);
  
  $score      = $exon->score()     || ".";
  $source     = $exon->source_tag  || "merged"; 
  $primary_tag= $exon->primary_tag || "exon";
  $frame      = $exon->phase       || ".";
  $strand = $exon->strand();
  if(! $strand) {
      $strand = ".";
  } elsif( $strand == 1 ) {
      $strand = '+';
  } elsif ( $exon->strand == -1 ) {
      $strand = '-';
  }
  $name        = $exon->seqname();
  
  $str = join("\t",
	      $name,
	      $source,
	      $primary_tag,
	      $exon->start(),
	      $exon->end(),
	      $score,
	      $strand,
	      $frame);
  
  $str .= "\t$trans_id";
  return $str;
}


############################################################

############################################################

sub exon_from_gff {
    my ($gffstring,$tag) = @_;
    
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

    $frame = 0 unless( $frame =~ /^\d+$/);
    $exon->phase($frame);

    #my $phase = ( 3 - $frame )%3;
    #$exon->phase($phase);
    #$exon->end_phase( ( $exon->phase + $exon->length)%3 );
    
    if ($score eq '.'){
	$exon->score(0);
    }
    elsif ( $score ){
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
    $tag{ $exon } = $group[0];
    return $exon;
}
