#!/usr/local/bin/perl -w

use strict;
use Getopt::Long;


my $in_dir;
my $out_dir;

&GetOptions( 
	     'in_dir:s'  => \$in_dir,
	     'out_dir:s' => \$out_dir,
	     );
unless( $in_dir && $out_dir ){
    print STDERR "$0 -in_dir -out_dir\n";
    exit(0);
}


open ( LS_QUERY, "ls $in_dir | ") or die("cannot do LS on $in_dir");
my @files;
while (<LS_QUERY>){
    chomp;
    push ( @files, $_ );
}

my $count_machine = 0;
my $count_files = 0;
my $command = "cat ";

foreach my $file (@files){
    my $this_file = $in_dir."/".$file;
    $command .= " $this_file ";
    $count_files++;
    print STDERR "file:$count_files machine:$count_machine\n";
    if ( $count_machine <12 ){
	if ( $count_files == 4 ){
	    $count_machine++;
	    $command .= " > $out_dir/catted"."_".$count_machine."_".$count_files;
	    system("$command");
	    $command = "cat ";
	    $count_files = 0;
	}
    }
    if ( $count_machine >= 12 ){
	if ( $count_files == 3 ){
	    $count_machine++;
	    $command .= " > $out_dir/catted"."_".$count_machine."_".$count_files;
	    system("$command");
	    $command = "cat ";
	    $count_files = 0;
	}
    }
}

		   
