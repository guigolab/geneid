package Geneid::Isocore;
use strict;
# use vars qw(@ISA);

############################################################
# 
# Author: T. Alioto
# Date: 060404
# POD documentation - main docs before the code

=head1 NAME

Geneid::Isocore - Handler for Isocore objects in GeneID parameter files

=head1 SYNOPSIS

    use Geneid::Param;
    use Geneid::Isocore;

    my ($profile,$length,$newoff) = getKmatrix($true_seqs,$false_seqs,$order,$offset);
    my $param = Geneid::Param->new(); #make a new parameter object
    $param->readParam($pin); #read in a template parameter file

    for (my $i = 0;$i < $param->numIsocores ; $i++){
	    if (!defined @{$param->isocores}[$i]->set_profile($label,$length,$newoff,$cutoff,$order,$profile)){die "error in setting profile\n";} #set a new profile where $profile is the profile name.
    }  

    $param->writeParam; # write out the parameter file


=head1 DESCRIPTION

Geneid::Isocore is a handler module for the isocores in a Geneid 
parameter file.  

It is preferable to make a new isocore from an old one, preferably 
through the Param module. Read a parameter file, change the necessary
profiles and parameters using the get, copy, set methods for profiles 
and the get/set access methods for all the other parameters. 

=head1 OBJECT METHODS

See below for more detailed summaries. 


=head1 AUTHOR - Tyler Alioto

Email: talioto@imim.es

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

my %profilenames = (
		    Start_profile => 1,
		    U12_Branch_point_profile => 1,
		    Branch_point_profile => 1,
		    Poly_Pyrimidine_Tract_profile => 1,
		    U12gtag_Acceptor_profile => 1,
		    U12atac_Acceptor_profile => 1,
		    Acceptor_profile => 1,
		    U12gtag_Donor_profile => 1,
		    U12atac_Donor_profile => 1,
		    U2gcag_Donor_profile => 1,
		    U2gta_Donor_profile => 1,
		    U2gtg_Donor_profile => 1,
		    U2gty_Donor_profile => 1,
		    Donor_profile => 1,
		    Stop_profile => 1,
		    PolyA_Signal_profile => 1,
		   );
sub new 
  {
      my ($class) = @_; 
      my $self = {
		  _boundaries_of_isochore => [0,100],
		  _Absolute_cutoff_exons => [qw(-150 -150 -150 -150 0)],
		  _Coding_cutoff_oligos => [qw(-100 -150 -150 -150)],
		  _Site_factor => [0.6,0.6,0.6,0.6,0.6],
		  _Exon_factor => [0.4,0.4,0.4,0.4],
		  _HSP_factor => [1.0,1.0,1.0,1.0,1.0],
		  _Evidence_Factor => 1,
		  _Evidence_Weight => 0,
		  _RSS_Markov_Score => undef,
		  _RSS_Donor_Score_Cutoff => undef,
		  _RSS_Acceptor_Score_Cutoff => undef,
		  _U12_Splice_Score_Threshold => undef,
		  _U12_Exon_Score_Threshold => undef,
		  _U12_Exon_weight => undef,
		  _Exon_weights => [-4,-4,-4,-4,0],
		  _Start_profile => undef,
		  _U12_Branch_point_profile => undef,
		  _Poly_Pyrimidine_Tract_profile => undef,
		  _U12gtag_Acceptor_profile => undef,
		  _U12atac_Acceptor_profile => undef,
		  _Acceptor_profile => undef,
		  _U12gtag_Donor_profile => undef,
		  _U12atac_Donor_profile => undef,
		  _Donor_profile => undef,
		  _PolyA_Signal_profile => undef,
		  _Stop_profile => {label=>'Stop_Profile',length=>4,offset=>0,cutoff=>-9999,order=>0,afactor=>"",bfactor=>"",acc_context=>"",dist=>"",optdist=>"",penalty_factor=>"",profile=>[
																							       [1,'A',0],
																							       [1,'C',0],
																							       [1,'G',0],
																							       [1,'T',0],
																							       [2,'A',-9999],
																							       [2,'C',-9999],
																							       [2,'G',-9999],
																							       [2,'T',0.000],
																							       [3,'A',0.000],
																							       [3,'C',-9999],
																							       [3,'G',0.000],
																							       [3,'T',-9999],
																							       [4,'A',0.000],
																							       [4,'C',-9999],
																							       [4,'G',0.000],
																							       [4,'T',-9999]
																							       ]},
		  _Markov_order => undef,
		  _Markov_Initial_probability_matrix => undef,
		  _Markov_Transition_probability_matrix => undef,
		  _maximum_number_of_donors_per_acceptor_site => 5
		 };
      bless $self, $class;
  }


#constructor
sub boundaries_of_isochore {
    my ( $self, $val ) = @_;
    $self->{_boundaries_of_isochore} = $val if defined($val);
    return $self->{_boundaries_of_isochore};
}
sub Absolute_cutoff_exons {
    my ( $self, $val ) = @_;
    $self->{_Absolute_cutoff_exons} = $val if defined($val);
    return $self->{_Absolute_cutoff_exons};
}
sub Coding_cutoff_oligos {
    my ( $self, $val ) = @_;
    $self->{_Coding_cutoff_oligos} = $val if defined($val);
    return $self->{_Coding_cutoff_oligos};
}
sub Site_factor {
    my ( $self, $val ) = @_;
    $self->{_Site_factor} = $val if defined($val);
    return $self->{_Site_factor};
}
sub Exon_factor {
    my ( $self, $val ) = @_;
    $self->{_Exon_factor} = $val if defined($val);
    return $self->{_Exon_factor};
}
sub HSP_factor {
    my ( $self, $val ) = @_;
    $self->{_HSP_factor} = $val if defined($val);
    return $self->{_HSP_factor};
}
sub Evidence_Factor {
    my ( $self, $val ) = @_;
    $self->{_Evidence_Factor} = $val if defined($val);
    return $self->{_Evidence_Factor};
}
sub Evidence_Weight {
    my ( $self, $val ) = @_;
    $self->{_Evidence_Weight} = $val if defined($val);
    return $self->{_Evidence_Weight};
}
sub RSS_Markov_Score {
    my ( $self, $val ) = @_;
    $self->{_RSS_Markov_Score} = $val if defined($val);
    return $self->{_RSS_Markov_Score};
}
sub RSS_Donor_Score_Cutoff {
    my ( $self, $val ) = @_;
    $self->{_RSS_Donor_Score_Cutoff} = $val if defined($val);
    return $self->{_RSS_Donor_Score_Cutoff};
}
sub RSS_Acceptor_Score_Cutoff {
    my ( $self, $val ) = @_;
    $self->{_RSS_Acceptor_Score_Cutoff} = $val if defined($val);
    return $self->{_RSS_Acceptor_Score_Cutoff};
}
sub U12_Splice_Score_Threshold {
    my ( $self, $val ) = @_;
    $self->{_U12_Splice_Score_Threshold} = $val if defined($val);
    return $self->{_U12_Splice_Score_Threshold};
}
sub U12_Exon_Score_Threshold {
    my ( $self, $val ) = @_;
    $self->{_U12_Exon_Score_Threshold} = $val if defined($val);
    return $self->{_U12_Exon_Score_Threshold};
}
sub U12_Exon_weight {
    my ( $self, $val ) = @_;
    $self->{_U12_Exon_weight} = $val if defined($val);
    return $self->{_U12_Exon_weight};
}

# sub U12gtag_Exon_weights {
#     my ( $self, $val ) = @_;
#     $self->{_U12gtag_Exon_weights} = $val if defined($val);
#     return $self->{_U12gtag_Exon_weights};
# }
# sub U12atac_Exon_weights {
#     my ( $self, $val ) = @_;
#     $self->{_U12atac_Exon_weights} = $val if defined($val);
#     return $self->{_U12atac_Exon_weights};
# }

=head2 Exon_weights

 Title   : Exon_weights
 Usage   : $isocore->Exon_weights(\@ew);
 Function:
 Example : $isocore->Exon_weights([split $line_with_exon_weights]); # set the exon weights, anon. array with 4 values
           $isocore->Exon_weights;   # get the exon weights
 Returns : an array reference
 Args    :

=cut

sub Exon_weights {
    my ( $self, $val ) = @_;
    $self->{_Exon_weights} = $val if defined($val);
    return $self->{_Exon_weights};
}
sub set_profile {
    my ( $self,$profile_name,$length,$offset,$cutoff,$order,$afactor,$bfactor,$acc_con,$dist,$optdist,$penalty,$profile ) = @_;
    return undef unless exists $profilenames{$profile_name};
    $self->{"_".$profile_name} = {label=>$profile_name,length=>$length,offset=>$offset,cutoff=>$cutoff,order=>$order,afactor=>$afactor,bfactor=>$bfactor,acc_context=>$acc_con,dist=>$dist,optdist=>$optdist,penalty_factor=>$penalty,profile=>$profile} if defined($profile_name);
    return $self->{"_".$profile_name};
}
sub copy_profile {
    my ( $self,$source ) = @_;
    $self->{"_".$source->{label}} = $source if defined($source);
    return $self->{"_".$source->{label}};
}
sub get_profile {
    my ( $self,$profile_name) = @_;
    return undef unless exists $profilenames{$profile_name};
    return $self->{"_".$profile_name};
}
sub Markov_order {
    my ( $self, $val ) = @_;
    $self->{_Markov_order} = $val if defined($val);
    return $self->{_Markov_order};
}
sub Markov_Initial_probability_matrix {
    my ( $self, $val ) = @_;
    $self->{_Markov_Initial_probability_matrix} = $val if defined($val);
    return $self->{_Markov_Initial_probability_matrix};
}
sub Markov_Transition_probability_matrix {
    my ( $self, $val ) = @_;
    $self->{_Markov_Transition_probability_matrix} = $val if defined($val);
    return $self->{_Markov_Transition_probability_matrix};
}
sub maximum_number_of_donors_per_acceptor_site {
    my ( $self, $val ) = @_;
    $self->{_maximum_number_of_donors_per_acceptor_site} = $val if defined($val);
    return $self->{_maximum_number_of_donors_per_acceptor_site};
}
1;
