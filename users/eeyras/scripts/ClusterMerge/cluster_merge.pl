#!/usr/bin/perl  -w

############################################################
#
# Release 0.2.1
#
# script to run the ClusterMerge Algorithm (Eduardo Eyras)
#
# written by Eduardo Eyras (eeyras@imim.es)
#
##########################################################################
#                                                                        #
#  This program is free software; you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation; either version 2 of the License, or     #
#  (at your option) any later version.                                   #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program; if not, write to the Free Software           #
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.             #
##########################################################################


######################################################################
#
#  If you use this program in your analyses, please cite:
#
#  Eyras E, Caccamo M, Curwen V, Clamp M.
#  ESTGenes: alternative splicing from ESTs in Ensembl.
#  Genome Res. 2004 May;14(5):976-87. 
#
######################################################################

use strict;
use Getopt::Long;
use ClusterMerge::ClusterMerge;
use ClusterMerge::Transcript;
use ClusterMerge::TranscriptUtils;
use ClusterMerge::Exon;
use ClusterMerge::GFFTools;

my $gff_format = 1;
my $input;
my $output;

############################################################
# deafult options

my $gff               = 1;
my $gtf               = 0;
my $comp_level        = 3;
my $splice_mismatch   = 10;
my $intron_mismatch   = 0;
my $min_order         = 1;
my $internal_splice_overlap = 10;
my $sets;
my $help;

############################################################

&GetOptions( 
	    #'gff_format:s'        => \$gff_format,
	     'input:s'                   => \$input,
	     'output:s'                  => \$output,
	     'comp_level:n'              => \$comp_level,
	     'splice_mismatch:n'         => \$splice_mismatch,
	     'intron_mismatch:n'         => \$intron_mismatch,
	     'internal_splice_overlap:n' => \$internal_splice_overlap,
	     'sets'                      => \$sets,
	     'help'                      => \$help,
	     );

if ( $help || !$input || !$output){
  &usage;
  exit(0);
}

  
############################################################
# READ INPUT FILE
#
open ( IN, "<$input") || die("could not open input file $input");

my @transcripts;
my %tag2transcript;

#my %forward_predictions;
#my %reverse_predictions;

while(<IN>){
    chomp;
    my @e = split;
    next unless ( $e[3] && $e[4] );
    my $exon;
    if ($gff){
	$exon = ClusterMerge::GFFTools->exon_from_gff ($_);
    }
    elsif($gtf){
	$exon = ClusterMerge::GFFTools->exon_from_gtf ($_);
    }
    $exon->source_tag("prediction");
    my $strand = $exon->strand;
    unless ( $tag2transcript{$exon->group_tag} ){
	my $trans_id = $exon->group_tag;
	$tag2transcript{$trans_id} =  ClusterMerge::Transcript->new();
	$tag2transcript{$trans_id}->dbID( $exon->transcript_tag );
	$tag2transcript{$trans_id}->type("prediction");
	#if ($strand == 1){
	#    push ( @{$forward_predictions{$exon->seqname}}, $tag2transcript{$exon->transcript_tag} );
	#}
	#else{
	#    push ( @{$reverse_predictions{$exon->seqname}}, $tag2transcript{$exon->transcript_tag} );
	#}
	push ( @transcripts, $tag2transcript{$exon->transcript_tag} );
    }
    $tag2transcript{$exon->group_tag}->add_Exon($exon);
}


############################################################

print STDERR "running ClusterMerge with parameters:\n";
print STDERR "comp_level              = $comp_level\n";
print STDERR "splice_mismatch         = $splice_mismatch\n";
print STDERR "intron_mismatch         = $intron_mismatch\n"; 
print STDERR "min_order               = $min_order\n";
print STDERR "internal_splice_overlap = $internal_splice_overlap\n";

############################################################

my $cluster_merge = 
    ClusterMerge::ClusterMerge->new(
				    -transcripts                   => \@transcripts,
				    -comparison_level              => $comp_level,
				    -splice_mismatch               => $splice_mismatch,
				    -intron_mismatch               => $intron_mismatch,
				    -minimum_order                 => $min_order,
				    -internal_splice_overlap       => $internal_splice_overlap,
				    );



#my $t1 = time;
$cluster_merge->run;
#my $t2 = time;
#print STDERR "TIME TO RUN CLUSTERMERGE: ".($t2-$t1)." *********************\n";

open ( OUT, ">$output") || die("could not open input file $output");

############################################################
# can retrieve the non-redundant sets:
if ( $sets ){
  
  # list of listrefs, each one cointaining a list of transcript objects
  my @sets = $cluster_merge->sub_clusters;
  
  my $count = 0;
  foreach my $set ( @sets ){
      $count++;
      print OUT "set $count:\n";
      foreach my $transcript ( @$set ){
	  foreach my $exon ( @{$transcript->get_all_Exons} ){
	      print OUT &gff_string($exon,$transcript)."\n";
	  }
      }
  } 
}

############################################################
# or the merged transcripts:
else{
    my @merged_transcripts = $cluster_merge->output;
    foreach my $transcript ( @merged_transcripts ){
	foreach my $exon ( @{$transcript->get_all_Exons} ){
	    print OUT &gff_string($exon,$transcript)."\n";
	}
    }
}

close (OUT);


############################################################

sub gff_string{
  my ($exon, $transcript) = @_;
  my ($str,$source,$primary_tag,$score,$frame,$name,$strand);
  
  if( $exon->can('score') ) {
    $score = $exon->score();
  }
  $score = '.' unless defined $score;
  
  if( $exon->can('frame') ) {
    $frame = $exon->frame();
  }
  $frame = '.' unless defined $frame;
  
  $strand = $exon->strand();
  if(! $strand) {
    $strand = ".";
  } elsif( $strand == 1 ) {
    $strand = '+';
  } elsif ( $exon->strand == -1 ) {
    $strand = '-';
  }
  
  $name        = $exon->seqname();
  $source      = "merged2";
  $primary_tag = "exon";
  
  $str = join("\t",
	      $name,
	      $source,
	      $primary_tag,
	      $exon->start(),
	      $exon->end(),
	      $score,
	      $strand,
	      $frame);
  
  my $tag_str = $transcript->type || $transcript->dbID;
  $str .= "\t$tag_str";
  
  my $transcript_label = '';
  my $exon_label = '';
  foreach my $e ( @{$exon->exon_evidence} ){
      $transcript_label .= $e->transcript_tag.":";
      if ( $e->dbID || $e->stable_id ){
	  $exon_label .= ($e->dbID || $e->stable_id).":";
      }
  }
  $str .= "\t$transcript_label\t$exon_label";
      
  return $str;
}


############################################################

sub usage {
    print STDERR <<EOF
    script to run the ClusterMerge algorithm
    Release 0.2.1
    Usage: cluster_merge.pl options
    Where options are:

    -input           : name of the input file in gff format

    -output          : name of the output file

    -comp_level      : comparison level (default = 3 )
                       1 --> strict: exact exon matching. 
                       Does not use any other parameteres passed in. Example:

                       #####-----#####-----#####
                       #####-----#####-----#####

                       2 --> allow edge exon mismatches. 
                       Uses the parameter 'exon_match' and 'internal_splice_overlap' if defined. 
                       Example:
    
                         ###-----#####-----#######
                       #####-----#####-----#####-----#####

                       3 ---> allow internal mismatches. 
                       Uses the parameters 'exon_match', 'splice_mismatch' and 'internal_splice_overlap' 
                       if defined. Example:

                       #####---########----######
                       #####-----######----####------#####

                       4 ---> allow intron mismatches. 
                       This one can use all the parameters if they have been defined.

                       ################----#######
                       #####-----######----####------#####
  
                       5 ---> loose match. It allows intron mismatches if so desired. There is no limitation on
                       the number of mismatches at the splice-sites. Examples:

                       #################----#######           and      #######------####----#######  
                       #####-----######----####------#####               #####-----######----####------#####
                    
                       would be merged as redundant.

     -splice_mismatch: maximum number of non-opverlapping nucleotides allowed in splice sites 
                       ( not used at comp_level = 5 )

     -intron_mismatch: maximum number of non-opverlapping nucleotides allowed in introns (default = 0 )

     -internal_splice_overlap: (default = 0 ) 
                       number of base pairs (N) we allow an external exon overlap
                       an intron in another transcript:
                                       |--N--|
                       ######-------##########
                      #######-------####-------------#######
    
     -exon_match     : TRUE if we want both transcripts to match 1-to-1 all their exons

     -min_order      : minimum number of transcripts required to be in a cluster to create a merged transcript 
                       ( default = 1 )

     -sets           : outputs the non redundant lists instead of the merged transcript 
                       (switched off by default)

     -gff            : input is in GFF format
     
     -gtf            : input is in GTF format

     -help           : outputs this help
EOF
}

