#!/usr/local/bin/perl -w

use strict;
use Bio::SeqIO;	
use Bio::Seq;
use Getopt::Long;

my $chunks;
my $file;
my $length;
my $outfile;

&GetOptions( 		
			'file:s'   => \$file,
			'length:s'  => \$length,
			'outfile:s' => \$outfile,
			'chunks:s' => \$chunks,	   
			);	

unless ( $file && $length && $outfile){
    print STDERR "script to split one fasta sequence into multifasta\n";
    print STDERR "Usage: $0 -file <file.fa> -length <length of chunks> -outfile <output file>\n";
    exit (0);
}


my $max_length = get_length( $file );

my $start = 1;
my $end   = $start + $length - 1;

#my @path = split '\/', $file;
#my $outfile = $path[-1]."_split";
#print "outfile = $outfile\n";

open (OUT, ">$outfile") or die ("cannot open file $outfile for writing");
my $out = Bio::SeqIO->new(-format =>'Fasta',
			  -fh     => \*OUT,
			  );

#print "start = $start - end = $end\n";
while( $start < $max_length ){
    
    if ( $end > $max_length ){
	$end = $max_length;
    }
    #print "INSIDE start = $start - end = $end\n";
    my $seq = get_subseq( $file, $start, $end );
        
    my $seq_object = Bio::Seq->new();
    $seq_object->seq($seq);
    $seq_object->display_id("$file"."_"."$start"."_"."$end");
    
    $out->write_seq($seq_object);

    
    $start += $length - 200;
    $end   = $start + $length - 1;
}
    
close(OUT);

############################################################

sub get_subseq{
    my ( $file, $start, $end ) = @_;
    
    my $command = "chr_subseq $file $start $end |";
    
    open( SEQ, $command ) || die("Error running command $command");
    my $seq = <SEQ>;
    chomp $seq;
    close( SEQ );
    return $seq;
}

############################################################

sub get_length{
    my ( $file ) = @_;
    
    my $command = "fastalength $file |";
    open( LEN, $command ) || die("Error running command $command");
    my $length = 0;
    while (<LEN>){
	chomp;
	
	my @entries  = split '\s+',$_;
	#print "entries[0] = $entries[0]\n";
	#print "entries[1] = $entries[1]\n";
	$length = $entries[0];
    }
    #print " Length = $length\n";
    return $length;
}
