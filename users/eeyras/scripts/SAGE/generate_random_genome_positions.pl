#!/usr/local/bin/perl -w

use strict;

# this program generates a number of random locations in the genome
# in GFF format. The number and size of the locations
# is given as input.

if (scalar(@ARGV) < 2) {
    print STDERR "$0 location-number location-size\n";
    exit(1);
}

my $locations = $ARGV[0]; 
my $size      = $ARGV[1];
my $chr_length_file = $ARGV[2];

my $seq_dir = "/seq/genomes/G.gallus/golden_path_200402/chromFa/";

my ($chrname,$start,$end,$strand,$id);
my @chrs;
my %chr_length;
my %chr_file;
my $longest = 0;

if ( $chr_length_file ){
    open (IN, "<$chr_length_file") or die("cannot open file $chr_length_file");
    while(<IN>){
	chomp;
	my ($chr_name,$chr_length) = split;
	$chr_length{$chr_name} = $chr_length;
	$longest = $chr_length{$chr_name} if ($chr_length{$chr_name} > $longest);
	$chr_file{$chr_name}   = $seq_dir."/".$chr_name.".fa";
	push (@chrs,$chr_name);
    }
    close(IN);
}
else{
    open ( LS, "ls -al $seq_dir | egrep chr | egrep -v M | ") or die("cannot open ls -al $seq_dir");
    
    while(<LS>){
	chomp;
	my @e = split;
	my $chr_name = $e[-1];
	$chr_name =~ s/\.fa//;
	#print "chr_name = $chr_name\n";
	push (@chrs, $chr_name);
	$chr_file{$chr_name}   = $seq_dir."/".$chr_name.".fa";
	#print "chr_file = $chr_file{$chr_name}\n";
	$chr_length{$chr_name} = get_chr_length( $chr_file{$chr_name} );
	#print "chr_length = $chr_length{$chr_name}\n";
	$longest = $chr_length{$chr_name} if ($chr_length{$chr_name} > $longest);
	
	print "$chr_name\t$chr_length{$chr_name}\n";
    }
    close(LS);
}
my $counter = 0;

############################################################
# we use the 'rejection sampling algorithm'

 LOCATION:
    while ( $counter < $locations ){
	
	############################################################
	# choose a random number  to select one chromosome
	my $rand = random_number(1,scalar(@chrs));
	my $chrname = $chrs[$rand-1];
	
	############################################################
	# select randomly the strand
	my $strand_rand  = random_number(1,2);
	my $strand;
	$strand = -1 if $strand_rand == 1;
	$strand = 1 if $strand_rand == 2;
	
	############################################################
	# select a start point in the chromosome using the length
	# of the longest, so that the probability is weighted by the
	# relative lengths
	my $start = random_number( 1, $longest - $size + 1);
	#print "$start is the random start between 1 and ".($longest - $size + 1)."\n";
	
	my $end   = $start + $size - 1;
	
	#print "random choice: $chrname.$start-$end:$strand\t";
	#print "rejected\n" if ( ( $start > $chr_length{$chrname} - $size +1) ||  has_Ns($chrname,$start,$end) );


	############################################################
	# reject if you fall off the chromosome:
	next LOCATION if ( $start > $chr_length{$chrname} - $size +1);
	
	############################################################
	# reject if the sequence picked has got any N's
	next LOCATION if ( has_Ns($chrname,$start,$end) );


	#print "accepted\n";
	############################################################
	# otherwise accept
	$counter++;
	my $string = join "\t", ($chrname,"random","tag",$start,$end,".",$strand,".",$counter);
	print $string."\n";
    }



############################################################

sub random_number{
    my ($min,$max) = @_;
    
    my $random = rand();
    
    # translate into random numbers between 1 and 10 (both included)
    my $x = int (($max - $min + 1)*$random)+1;

    return $x;
}
    
############################################################

sub get_chr_length{
    my ( $chr_file ) = @_;
    
    my $command = "fastalength $chr_file |";
    open( LEN, $command ) || die("Error running command $command");
    my $length = 0;
    while (<LEN>){
	chomp;
	my @entries  = split '\s+',$_;
	$length = $entries[0];
	
    }
    return $length;
}
 
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

sub has_Ns{
    my ($chrname,$start,$end) = @_;
    my $file = $chr_file{$chrname};
    my $seq = get_subseq($file, $start, $end );
    
    my @Ns = $seq =~/n/ig;
    
    #print "Ns = @Ns\n";
    return 1 if ( @Ns);
    return 0;
}
 
############################################################   
