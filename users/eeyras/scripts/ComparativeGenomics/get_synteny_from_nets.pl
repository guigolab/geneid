#!/usr/local/bin/perl -w

use strict;

#bin    level   tName   tStart  tEnd    strand  qName   qStart  qEnd    chainId ali     score   qOver   qFar    qDup    type    tN      qN      tR      qR      tNewR   qNewR   tOldR   qOldR   tTrf     qTrf
#0       1       chr7    92403592        149996044       +       chr6    3063416 48830107        14      16487834        378348525       -1      -1      620584  top     25000   1556412 2634070517791422 13076051        13889658        9416204 1515610 746802  1002934
#1465    2       chr7    115400414       115406540       +       chr6    16773070        16777788        0       0       0       -1      -1      -1      gap     0       0       1899    2167    0        1295    1144    0       0       95
#1

#cut -f3,4,5,6,7,8,9 table.net.test | more
##tName  tStart  tEnd    strand  qName   qStart  qEnd
#chr7    92403592        149996044       +       chr6    3063416 48830107
#chr7    115400414       115406540       +       chr6    16773070        16777788
#chr7    115406740       115406751       +       chr6    16777990        16777990


my %chr_features;
while (<>){
    chomp;
    my @e = split;
    next unless $e[2] =~/chr/;
    
    # ($target_chr, $target_start, $target_end, $query_chr, $query_start, $query_end, $query_strand) 
    
    my $f = [ $e[2],$e[3],$e[4],$e[6],$e[7],$e[8],$e[5] ];
    push (@{$chr_features{$e[2]}{$e[6]} }, $f);
}    

############################################################
# chain features into longer ones:

my $length = 1000;

my @slices;
foreach my $target_chr ( keys %chr_features ){
    
    my @target_slices;
    foreach my $query_chr ( keys %{$chr_features{$target_chr}} ){
	
	my @chr_features = sort { $a->[1] <=> $b->[1] } @{ $chr_features{$target_chr}{$query_chr} };
    
	my ($target_chr, $target_start, $target_end, $query_chr, $query_start, $query_end, $query_strand) = @{$chr_features[0]};
	foreach (my $i=1; $i<scalar(@chr_features); $i++ ){
	    
	    if ( $chr_features[$i]->[1] - $target_end - 1 > 10000
		 ||
		 !($chr_features[$i]->[3] eq $query_chr)
		 ||
		 !($chr_features[$i]->[6] eq $query_strand)
		 ){
		
		push (@target_slices, [$target_chr, $target_start, $target_end, $query_chr, $query_start, $query_end, $query_strand]);
		my $s = join "\t", ($target_chr, $target_start, $target_end, $query_chr, $query_start, $query_end, $query_strand);
		print $s."\n";
		($target_chr, $target_start, $target_end, $query_chr, $query_start, $query_end, $query_strand) = @{$chr_features[$i]};
	    }
	    if ( $chr_features[$i]->[2] > $target_end ){
		$target_end = $chr_features[$i]->[2];
	    }
	    if ( $chr_features[$i]->[4] < $query_start ){
		$query_start = $chr_features[$i]->[4];
	    }
	    if ( $chr_features[$i]->[5] > $query_end ){
		$query_end = $chr_features[$i]->[5];
	    }
	    
	}
	push (@target_slices, [$target_chr, $target_start, $target_end, $query_chr, $query_start, $query_end, $query_strand]);
	my $s = join "\t", ($target_chr, $target_start, $target_end, $query_chr, $query_start, $query_end, $query_strand);
	print $s."\n";
    }

    # take the longest non-overlapping regions:
    my @slices = sort { length($b) <=> length($a) } @target_slices;
    my @chosen_slices;
    push (@chosen_slices, shift @slices);
  SLICE:
    foreach my $slice ( @slices ){
	
	my $found = 0;
	foreach my $chosen (@chosen_slices){
	    
	    if ( !($slice->[2]< $chosen->[1]
		   ||
		   $slice->[1] > $chosen->[2]
		   )){
		my $found = 1;
		next SLICE;
	    }
	}
	if ($found == 0){
	    push (@chosen_slices,$slice);
	}
    }

    foreach my $slice (@chosen_slices){
	my $s = join "\t", @$slice;
	print "CHOSEN$s"."\n";
    }
}
 

sub length{
    my $slice = shift;
    return $slice->[2] - $slice->[1] + 1;
}
