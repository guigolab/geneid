#!/usr/local/bin/perl -w

use strict;
use Getopt::Long;

my $dir;

&GetOptions( 
			 'dir:s'   => \$dir,
			 );

unless( $dir ){
	print STDERR "Usage: $0 -dir <dir_name>\n";
	exit(0);
}

############################################################
open ( LS, "ls $dir | ") or die("cannot do LS on $dir");

while (<LS>){
    chomp;
	my $chr = $_;
	my $chr_file = $dir."/".$chr;
	my $length = get_chr_length( $chr_file );
	
	$chr =~ s/\.fa//;
	print STDOUT "$chr\t$length\n";
}
close(LS);


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
