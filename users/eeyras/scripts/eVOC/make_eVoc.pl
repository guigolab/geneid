#!/usr/bin/perl

# Before parsing the files, it is good to 
# check any hidden tab-like character
# which could hinder the parsing:

#                     NUL                  \0
#                     <alert character>    \a
#                     <backspace>          \b
#                     <form-feed>          \f
#                     <newline>            \n
#                     <carriage return>    \r
#                     <tab>                \t
#                     <vertical tab>       \v

# with 
# od -c filename
# one can check for this

use strict;

############################################################
# files with the annotations

my @in_files = ('anatomicalsystem_annotations.tsv',
		'celltype_annotations.tsv',
		'developmentalstage_annotations.tsv',
		'pathology_annotations.tsv'); 

############################################################
# file with the map of clones to ESTs

my $in_clone2ests ="../EST2CloneLibrary.tsv";

############################################################
# output files

my @out_files = ('EST_anatomicalsystem.txt',
		 'EST_celltype.txt',
		 'EST_developmentalstage.txt',
		 'EST_pathology.txt'); 

my $out_clone2ests = "Clone2ESTs.txt";

############################################################
# read the ESTs for each clone 
# put it also into a tab delimited file:

my %clone2est;
open (IN,"<$in_clone2ests");
open (OUT,">$out_clone2ests");
while(<IN>){
    chomp;
    my @list = split '\[';
    
    shift @list;
    
    my $clone = shift @list;
    $clone =~ s/\]//;
    my ($c,$clone_lib) = split ":", $clone;

    my $id = shift @list;
    $id =~ /dbEST\s+ID\:(\d+)/;
    my $dbEST_id = $1;

    foreach my $e ( @list ){
	
	$e =~ /\:([a-zA-Z0-9]+)\]/;
	my $est_id = $1;
	
	push( @{$clone2est{$clone_lib}}, $est_id );

	print OUT $clone_lib."\t";
	print OUT $dbEST_id."\t";
	print OUT $est_id."\n";
    }
}
#####################################
close(OUT);
close(IN);


############################################################
# write each annotation file as 
# 'est_id'  'clone_name'  'clone_id'  'colon-separated-annotation'
for(my $k=0; $k<4; $k++ ){
  
    open (IN,"<$in_files[$k]");
    open (OUT,">$out_files[$k]");
    my $clone_lib;
    my $dbEST_id;
    my @annotations = ();
    
  READ_LINE:
    while(<IN>){
	
	############################################################
	# eliminate the \n and possible \r characters
	# the modifier 'o' precompiles the regex for efficiency
	$_ =~ s/\r?\n\r?$//o; # chomp;
	my @list = split '\t', $_;
	
      COLUMN:
	for(my $i=0; $i< @list; $i++){
	    
	    if ( $list[$i] =~/Clone Library\[Name:(.*?)\]\[dbEST ID:(\d+)\]/ ){
		
		$clone_lib= $1;
		$dbEST_id = $2;		
		foreach my $est_id ( @{$clone2est{$clone_lib}} ){
		    print OUT $est_id."\t";
		    print OUT $clone_lib."\t".$dbEST_id."\t";
		    for(my $j=0; $j<$i; $j++){
			print OUT "$annotations[$j]:";
		    }
		    print OUT "\n";
		}
		next READ_LINE;
	    }
	    
	    ############################################################
	    # update the paths (this keeps the common nodes)
	    if ($list[$i]){
		$annotations[$i] = $list[$i];
	    }
	    
	}
    }
    close(OUT);
    close(IN);
    
}


