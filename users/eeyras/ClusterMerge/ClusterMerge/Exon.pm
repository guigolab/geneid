# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code

=pod 

=head1 NAME

ClusterMerge::Exon

=head1 SYNOPSIS

    $ex = ClusterMerge::Exon->new();

    $ex->start(10);
    $ex->end(100);

Examples of creating an exon

    # start = 1208, end = 1506, forward strand
    $ex = Exon->new(1208,1506,1) 
    
    Start and end coordinates are always stored with start < end. If they are 
    input in the reverse order they will be swapped over.  The value for the 
    strand will be kept as its input value;

    Strand values:  + or  1 = forward strand
                    - or -1 = reverse strand
                    . or  0 = unknown strand

    $ex->phase(0);         # Sets the phase of the exon
    $ex->end_phase(1);      # sets the end_phase of the exon

    Phase values  are 0,1,2


=head1 DESCRIPTION

Exon object.  

=head1 CONTACT

eae@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal 
methods are usually preceded with a_

=cut


# Let the code begin...


package ClusterMerge::Exon;
use vars qw(@ISA);
use strict;
use ClusterMerge::Root;


@ISA = qw(ClusterMerge::Root);


=head2 new

=cut

sub new {
    my($caller,@args) = @_;
    
    my $self = {};

    if(ref $caller) {
	bless $self, ref $caller;
    } 
    else {
	bless $self, $caller;
    }
    
    $self->{'_gsf_tag_hash'} = {};
    $self->{'_gsf_sub_array'} = [];
    $self->{'_parse_h'} = {};
    $self->{'_is_splittable'} = 0;
    
    my ($start,$end,$strand,$frame,$score,$analysis,$seqname, $source_tag,
	$primary_tag, $percent_id, $p_value, $phase, $end_phase) =
	    
	    $self->_rearrange([qw(START
				  END
				  STRAND
				  FRAME
				  SCORE
				  ANALYSIS
				  SEQNAME
				  SOURCE_tag
				  PRIMARY_TAG
				  PERCENT_ID
				  P_VALUE
				  PHASE
				  END_PHASE
				  )],@args);
    
    #  $gff_string && $self->_from_gff_string($gff_string);
    
    if ( defined $analysis  && $analysis ne "")   { $self->analysis($analysis)};
    if ( defined ($start) && $start ne "" )       { $self->start($start)};
    if ( defined ($end )  && $end   ne "" )       { $self->end($end)}
    if ( defined $strand  && $strand ne "")       { $self->strand($strand)}
    if ( defined $frame  && $frame ne "")         { $self->frame($frame)}
    if ( defined $score  && $score ne "")         { $self->score($score)}
    if ( defined $seqname && $seqname ne "")      { $self->seqname($seqname)};
    if ( defined $percent_id && $percent_id ne ""){ $self->percent_id($percent_id)};
    if ( defined $p_value && $p_value ne "")      { $self->p_value($p_value)};
    if ( defined $phase && $phase ne "")          { $self->phase($phase)};
    if ( defined $end_phase && $end_phase ne "")  { $self->end_phase($end_phase)};
    
    return $self;
}

=head2 end_phase

  Arg [1]    : (optional) int $end_phase
  Example    : $end_phase = $feat->end_phase;
  Description: Gets/Sets the end phase of the exon.
               end_phase = number of bases from the last incomplete codon of 
               this exon.
               Usually, end_phase = (phase + exon_length)%3
               but end_phase could be -1 if the exon is half-coding and its 3 
               prime end is UTR.
  Returntype : int
  Exceptions : warning if end_phase is called without an argument and the
               value is not set.
  Caller     : general

=cut

sub end_phase {
  my ($self,$endphase) = @_;
  if ( defined($endphase) ){
    $self->{_end_phase} = $endphase;
  }
  if ( !defined( $self->{_end_phase} ) ){
    $self->throw("No end phase set in Exon. You must set it explicitly. $!" .
	      "Caller: ".caller);
  }
  return $self->{_end_phase};
}

=pod

=head2 phase

  my $phase = $exon->phase;
  $exon->phase(2);

Get or set the phase of the Exon, which tells the
translation machinery, which makes a peptide from
the DNA, where to start.

The Ensembl phase convention can be thought of as
"the number of bases of the first codon which are
on the previous exon".  It is therefore 0, 1 or 2
(or -1 if the exon is non-coding).  In ascii art,
with alternate codons represented by B<###> and
B<+++>:

       Previous Exon   Intron   This Exon
    ...-------------            -------------...

    5'                    Phase                3'
    ...#+++###+++###          0 +++###+++###+...
    ...+++###+++###+          1 ++###+++###++...
    ...++###+++###++          2 +###+++###+++...

Here is another explanation from Ewan:

Phase means the place where the intron lands
inside the codon - 0 between  codons, 1 between
the 1st and second base, 2 between the second and
3rd  base. Exons therefore have a start phase and
a end phase, but introns have just one phase.

=cut

sub phase {
  my ($self,$value) = @_;
  
  if (defined($value)) {
    # Value must be 0,1,2, or -1 for non-coding
    if ($value =~ /^(-1|0|1|2)$/) {
      #print STDERR "Setting phase to $value\n";
      $self->{'phase'} = $value;
    } else {
      $self->throw("Bad value ($value) for exon phase. Should only be" .
		   " -1,0,1,2\n");
    }
  }
  return $self->{'phase'};
}



=head2 frame

=cut

sub frame {
    my ($self,$value) = @_;
    
    if( defined $value ) {
	$self->{_frame} = $value;
    }
    return $self->{_frame};
    
    
    # frame is mod 3 of the translation point
    if( $self->phase == -1 ) {
	return '.'; # gff convention for no frame info
    }
    if( $self->phase == 0 ) {
	return $self->start%3;
    }
    
    if( $self->phase == 1 ) {
	return ($self->start+2)%3;
    }
    
    if( $self->phase == 2 ) {
	return ($self->start+1)%3;
    }
    
    $self->throw("bad phase in exon ".$self->phase);
    
}



=head2 type

  Arg [1]    : (optional) $value
  Example    : Gets/Sets th etype of this exon
  Description: Returns the type of the exon (Init, Intr, Term)
  Returntype : string
  Exceptions : none

=cut

sub type {
  my ($self,$value) = @_;
  
  if (defined($value)) {
    $self->{'type'} = $value;
  }
  return $self->{'type'};
}


=head2 stable_id

 Title   : stable_id
 Usage   : $obj->stable_id
 Function: 
 Returns : value of stable_id
 Args    : 


=cut

sub stable_id{

    my ($self,$value) = @_;
    

    if( defined $value ) {
      $self->{'_stable_id'} = $value;
    }

    return $self->{'_stable_id'};

}

sub dbID{

    my ($self,$value) = @_;
    

    if( defined $value ) {
      $self->{_dbID} = $value;
    }
    return $self->{_dbID};

}


# Inherited methods
# but you do have all the SeqFeature documentation: reproduced here
# for convenience...

=pod

=head1 Methods inherited from SeqFeature

=head2 start

 Title   : start
 Usage   : $start = $feat->start
 Function: Returns the start coordinate of the feature
 Returns : integer
 Args    : none

=cut

sub start{
    my ($self,$value) = @_;

    if (defined($value)) {
        if ($value !~ /^\-?\d+/ ) {
        $self->throw("$value is not a valid start");
    }
	
    $self->{'_gsf_start'} = $value
   }

    return $self->{'_gsf_start'};

}

=head2 end

 Title   : end
 Usage   : $end = $feat->end
 Function: Returns the end coordinate of the feature
 Returns : integer
 Args    : none

=cut

sub end{
    my ($self,$value) = @_;

    if (defined($value)) {
        if( $value !~ /^\-?\d+/ ) {
            $self->throw("[$value] is not a valid end");
        }
	
        $self->{'_gsf_end'} = $value;
    }

   return $self->{'_gsf_end'};
}


=head2 strand

 Title   : strand
 Usage   : $strand = $feat->strand()
           $feat->strand($strand)
 Function: get/set on strand information, being 1,-1 or 0
 Returns : -1,1 or 0
 Args    : none


=cut

sub strand {
    my ($self,$value) = @_;

    if (defined($value)) {
        if( $value eq '+' ) { $value = 1; }
        if( $value eq '-' ) { $value = -1; }
        if( $value eq '.' ) { $value = 0; }

        if( $value != -1 && $value != 1 && $value != 0 ) {
            $self->throw("$value is not a valid strand info");
        }
        $self->{'_gsf_strand'} = $value;
    }

    return $self->{'_gsf_strand'};
}



=head2 length

 Title   : length
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :

=cut 

sub length{
   my ($self,@args) = @_;

   return $self->end - $self->start +1;
}

############################################################

=head2 gff_string

 Title   : gff_string
 Usage   : $str = $feat->gff_string
 Function: provides the feature information in GFF
           version 2 format.
 Returns : A string
 Args    : None

=cut

sub gffstring {
   my ($self) = @_;

   my $str;

   my $strand = "+";
   
   if ((defined $self->strand)&&($self->strand == -1)) {
     $strand = "-";
   }
   
   $str .= (defined $self->seqname)     ?   $self->seqname."\t"      :  "unknown\t";
   $str .= (defined $self->source_tag)  ?   $self->source_tag."\t"   :  "merged\t";
   $str .= (defined $self->primary_tag) ?   $self->primary_tag."\t"  :  "exon\t";
   $str .= (defined $self->start)       ?   $self->start."\t"        :  ".\t";
   $str .= (defined $self->end)         ?   $self->end."\t"          :  ".\t";
   $str .= (defined $self->score)       ?   $self->score."\t"        :  ".\t";
   $str .= (defined $self->strand)      ?   $strand."\t"             :  ".\t";
   $str .= (defined $self->phase)       ?   $self->phase."\t"        :  ".\t";
   eval{
     $str .= (defined $self->end_phase) ?   $self->end_phase."\t"        :  ".\t";
   };
   $str .= (defined $self->group_tag)   ?   $self->group_tag."\t"        :  ".\t";
   return $str;
}

############################################################

=head2 group_tag

 Title   : group_tag
 Usage   : $tag = $exon->group_tag()
           $exon->group_tag('chr1_123');
 Function: Store/Returns the transcript/group tag for an exon
 Returns : a string
 Args    : none


=cut

sub group_tag{
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{_group_tag} = $arg;
    }
    
    return $self->{_group_tag};
}

############################################################

=head2 source_tag

 Title   : source_tag
 Usage   : $tag = $feat->source_tag()
           $feat->source_tag('genscan');
 Function: Returns the source tag for a feature,
           eg, 'genscan'
 Returns : a string
 Args    : none


=cut

sub source_tag{
    my ($self,$arg) = @_;

    if (defined($arg)) {
        $self->{_source_tag} = $arg;
    }
    
    return $self->{_source_tag};
}



=head2 primary_tag

 Title   : primary_tag
 Usage   : $tag = $feat->primary_tag()
           $feat->primary_tag('exon')
 Function: get/set on the primary tag for a feature,
           eg 'exon'
 Returns : a string
 Args    : none


=cut

sub primary_tag{
  my ($self,$arg) = @_;
  
  if (defined($arg)) {
      $self->{_primary_tag} = $arg;
  }
  return $self->{_primary_tag};
}

=head2 score

 Title   : score
 Usage   : $score = $feat->score()
           $feat->score($score)
 Function: get/set on score information
 Returns : float
 Args    : none if get, the new value if set


=cut

sub score {
    my ($self,$value) = @_;

    if(defined ($value) ) {
	if( $value !~ /^[+-]?\d+\.?\d*(e-\d+)?/ ) {
	    $self->warn("'$value' is not a valid score - setting it to zero");
	    $value = 0;
	}
	$self->{'_gsf_score'} = $value;
    }
    
    return $self->{'_gsf_score'};
}

=head2 has_tag
    
 Title   : has_tag
 Usage   : $value = $self->has_tag('some_tag')
 Function: Returns the value of the tag (undef if
           none)
 Returns :
 Args    :


=cut

sub has_tag{
   my ($self,$tag) = (shift, shift);

   return exists $self->{'_gsf_tag_hash'}->{$tag};
}

=head2 add_tag_value

 Title   : add_tag_value
 Usage   : $self->add_tag_value('note',"this is a note");
 Returns : nothing
 Args    : tag (string) and value (any scalar)


=cut

sub add_tag_value{
   my ($self,$tag,$value) = @_;

   if( !defined $self->{'_gsf_tag_hash'}->{$tag} ) {
       $self->{'_gsf_tag_hash'}->{$tag} = [];
   }

   push(@{$self->{'_gsf_tag_hash'}->{$tag}},$value);
}

=head2 each_tag_value

=cut

sub each_tag_value {
   my ($self,$tag) = @_;
   if( ! exists $self->{'_gsf_tag_hash'}->{$tag} ) {
       $self->throw("asking for tag value that does not exist $tag");
   }

   return @{$self->{'_gsf_tag_hash'}->{$tag}};
}


=head2 all_tags

 Title   : all_tags
 Usage   : @tags = $feat->all_tags()
 Function: gives all tags for this feature
 Returns : an array of strings
 Args    : none


=cut

sub all_tags{
   my ($self,@args) = @_;

   return keys %{$self->{'_gsf_tag_hash'}};
}


sub transcript_id{
    my ($self,$id) = @_;
    if ( $id ){
	$self->{_transcript_id} = $id;
    }
    return $self->{_transcript_id};
}




=head2 seqname

  Arg [1]    : string $seqname
  Example    : $seqname = $self->seqname();
  Description: Obtains the seqname of this features sequence.  This is set
               automatically when a sequence with a name is attached, or may
               be set manually.
  Returntype : string
  Exceptions : none
  Caller     : general, attach_seq

=cut

sub seqname{
   my ($self,$seqname) = @_;

   if(defined $seqname) {
       $self->{_seqname} = $seqname;
   } 
   return $self->{_seqname};
}


=head1 Range methods

=head2 overlaps

  Title   : overlaps
  Usage   : if($feat->overlaps($r)) { do stuff }
            if($feat->overlaps(200)) { do stuff }
  Function: tests if $feat overlaps $r
  Args    : a RangeI to test for overlap with, or a point
  Returns : true if the Range overlaps with the feature, false otherwise

=cut 

sub overlaps{
    my ( $self, $range ) = @_;
    if ( $self->strand == $range->strand &&
	 !( $self->start > $range->end || $self->end < $range->start ) ){
	return 1;
    }
    return 0;
}



=head2 contains

  Title   : contains
  Usage   : if($feat->contains($r) { do stuff }
  Function: tests whether $feat totally contains $r
  Args    : a RangeI to test for being contained
  Returns : true if the argument is totaly contained within this range


=head2 equals

  Title   : equals
  Usage   : if($feat->equals($r))
  Function: test whether $feat has the same start, end, strand as $r
  Args    : a RangeI to test for equality
  Returns : true if they are describing the same range


=head1 Geometrical methods

These methods do things to the geometry of ranges, and return
triplets (start, stop, strand) from which new ranges could be built.

=cut

=head2 intersection

  Title   : intersection
  Usage   : ($start, $stop, $strand) = $feat->intersection($r)
  Function: gives the range that is contained by both ranges
  Args    : a RangeI to compare this one to
  Returns : nothing if they don''t overlap, or 
            a new exon based on the range that they do overlap
=cut

sub intersection{
    my ($exon1, $exon2 ) = @_;

    my ($start1,$end1) = ( $exon1->start, $exon1->end );
    my ($start2,$end2) = ( $exon2->start, $exon2->end );
    
    my $end   = $exon1->min( $end1,  $end2 );
    my $start = $exon1->max( $start1, $start2);

    return ($start,$end);
}

sub min{
    my ($self,$min, @values) = @_;
    foreach my $v (@values){
	if ($v < $min){
	    $min = $v;
	}
    }
    return $min;
}

sub max{
    my ($self,$max, @values) = @_;
    foreach my $v (@values){
	if ($v > $max){
	    $max = $v;
	}
    }
    return $max;
}



=head2 union

  Title   : union
  Usage   : ($start, $stop, $strand) = $feat->union($r);
          : ($start, $stop, $strand) = Bio::RangeI->union(@ranges);
  Function: finds the minimal range that contains all of the ranges
  Args    : a range or list of ranges to find the union of
  Returns : the range containing all of the ranges

=cut



1;





