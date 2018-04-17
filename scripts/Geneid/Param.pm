package Geneid::Param;
use strict;
use vars qw(@ISA);
use Geneid::Isocore;
use Geneid::GeneModel;
#@ISA = qw(Bio::Root::Root);

############################################################
# 
# Author: T. Alioto
# Date: 060404

############################################################
					##### Strings #####
my $head = "# geneid parameter file";
my $label_comments = "# Comment lines must start with '#'\n";
my $label_noscore = "# Non-homology penalty\n";
my $label_numiso = "# Number of isochores\n";
my $label_isohead = "# PARAMETERS FROM ISOCHORE 1\n";
my $label_gc = "# \%GC\n";
my $label_exon_score_thresh = "# Exons score: cutoffs\n";
my $label_exon_score_factors = "# Exon score: factors\n";
my $label_profiles = "# Site prediction: Position Weight Arrays\n# Length, Offset, Cutoff and Order (Markov model)\n";
my $label_prof_props = "#len offset cutoff order a b acc_context min_dist opt_dist pen_scale\n";
my $label_profile_array = "# Transition probabilities at every position\n";
my $label_coding_potential = "# Exon prediction: Markov model\n";
my $label_initial_prob = "# Initial probabilities at every codon position\n";
my $label_transition_prob = "# Transition probabilities at every codon position\n";
my $label_donors_per_acceptor = "# Donors per acceptor to build exons\n";
my $label_gene_model = "# GENE MODEL: Rules about gene assembling (GenAmic)\n";
my $blank_line = "\n";

my $header_bkgd_sub_flank_len = "BKGD_SUBTRACT_FLANK_LENGTH\n";

############################################################

sub new
  {
      my ($class,$species) = @_;
      if ($species) {
	  $head = "# geneid parameter file: $species";
      }
      my $self = {
		  _head => "$head\n",
		  _noScore => 0,
		  _bkgdSubFlankLen => 500,
		  _numIsocores => 1,
		  _isocores => undef,
		  _geneModel => Geneid::GeneModel->new()
		 };
      bless $self, $class;
  }

sub readParam
  {
      my ( $self, $paramfile ) = @_;
      open (P, "<$paramfile") or die "Couldn't open parameter file $paramfile: $!\n";
      my @isocores;
      my $NOSCORE  = 0;
      my $BKGD_SUB_FLANK_LEN = 0;
      my $gm = "";
      my $line;
      #$line = <P>;
      #$head = $line;
      #$self->{_head}=$head;
      while ($line = <P>) {
	  next if $line =~/^#/;
	  chomp $line;
	  if ($line =~ /NO_SCORE/) {
	      $NOSCORE = <P>;
	      $NOSCORE =~ s/\s//g;
	  } elsif ($line =~ /$header_bkgd_sub_flank_len/) {
	      $BKGD_SUB_FLANK_LEN = <P>;
	      $BKGD_SUB_FLANK_LEN =~ s/\s//g;
	  } elsif($line =~ /number_of_isochores/) {
	      my $numiso = <P>;
	      $numiso =~ s/\s//g;	
	      $self->{_numIsocores} = $numiso;
	      for (my $iso = 1; $iso <= $numiso; $iso++) {
		  my $maxdonors = 0;
		  my @ary;
		  my $isocore = Geneid::Isocore->new();
		  while (<P>) {
		      next if m/^#/;
		      last if $maxdonors;
		      if (/boundaries_of_isochore/i) {
			  $_ = <P>;
			  chomp;
			  #print STDERR $_;
			  #@ary = split;
			  $isocore->boundaries_of_isochore([split]);
		      } elsif (/cutt?off_exons/i) {
			  $_ = <P>;
			  chomp;
			  $isocore->Absolute_cutoff_exons([split]);
		      } elsif (/cutoff_oligos/i) {
			  $_ = <P>;
			  chomp;
			  $isocore->Coding_cutoff_oligos([split]);
		      } elsif (/Site_factor/i) {
			  $_ = <P>;
			  chomp;
			  $isocore->Site_factor([split]);
		      } elsif (/Exon_factor|Weight_of_oligos_in_exons/i) {
			  $_ = <P>;
			  chomp;
			  $isocore->Exon_factor([split]);
		      } elsif (/HSP_factor/i) {
			  $_ = <P>;
			  chomp;
			  $isocore->HSP_factor([split]);
		      } elsif (/RSS_Markov_Score/i) {
			  $_ = <P>;
			  chomp;
			  s/\s//g;
			  $isocore->RSS_Markov_Score($_);
		      } elsif (/Evidence_Factor/i) {
			  $_ = <P>;
			  chomp;
			  s/\s//g;
			  $isocore->Evidence_Factor($_);
		      } elsif (/Evidence_Weight/i) {
			  $_ = <P>;
			  chomp;
			  s/\s//g;
			  $isocore->Evidence_Weight($_);
		      } elsif (/RSS_Donor_Score_Cutoff/i) {
			  $_ = <P>;
			  chomp;
			  s/\s//g;
			  $isocore->RSS_Donor_Score_Cutoff($_);
		      } elsif (/RSS_Acceptor_Score_Cutoff/i) {
			  $_ = <P>;
			  chomp;
			  s/\s//g;
			  $isocore->RSS_Acceptor_Score_Cutoff($_);
		      } elsif (/U12_Splice_Score_Threshold/i) {
			  $_ = <P>;
			  chomp;
			  s/\s//g;
			  $isocore->U12_Splice_Score_Threshold($_);
		      } elsif (/U12_Exon_Score_Threshold/i) {
			  $_ = <P>;
			  chomp;
			  s/\s//g;
			  $isocore->U12_Exon_Score_Threshold($_);
		      } elsif (/U12_Exon_weight/i) {
			  $_ = <P>;
			  chomp;
			  s/\s//g;
			  $isocore->U12_Exon_weight($_);
			  # 					} elsif  (/U12atac_Exon_weights/i){
			  # 						$_ = <P>;
			  # 						chomp;
			  # 						$isocore->U12atac_Exon_weights([split]);
		      } elsif (/(Exon_weights|Exon_weigths)/i) {
			  $_ = <P>;
			  chomp;
			  $isocore->Exon_weights([split]);
		      } elsif (/(\w+)_profile/i) {
			  my $profile_name = $1."_profile";
			  while ($_ = <P>) {
			      last if $_ !~ m/^#/;
			  }
			  chomp;
			  my ($length,$offset,$cutoff,$order,$afactor,$bfactor,$acc_con,$dist,$optdist,$penalty) = split; # if $profile_name ne "Stop_profile";
			  $afactor = 0 if !$afactor;
			  $bfactor = 1 if !$bfactor;
			  $acc_con = "" if !$acc_con;
			  $dist = "" if !$dist;
			  $optdist = "" if !$optdist;
			  $penalty = "" if !$penalty;
			  my @profile;
			  while (<P>) {
			      next if m/^#/;
			      last if m/^\s/;
			      last if m/^[^\d]/;
			      chomp;
			      my @e = split;
			      push @profile, \@e;
			  }
			  $isocore->set_profile($profile_name,$length,$offset,$cutoff,$order,$afactor,$bfactor,$acc_con,$dist,$optdist,$penalty,\@profile);
		      } elsif (/Markov_oligo_logs_file|Markov_order/i) {
			  $_ = <P>;
			  chomp;
			  s/\s//g;
			  $isocore->Markov_order($_);
		      } elsif (/Markov_Initial_probability_matrix/i) {
			  my @profile2;
			  while (<P>) {
			      last if  m/^\s/;
			      last if m/^[^ACGTacgt]/;
			      next if m/^#/;
			      chomp;
			      my @f = split;
			      push @profile2, \@f;
			  }
			  $isocore->Markov_Initial_probability_matrix(\@profile2);
		      } elsif (/Markov_Transition_probability_matrix/i) {
			  my @profile3;
			  while (<P>) {
			      last if m/^\s/;
			      last if m/^[^ACGTacgt]/;
			      next if m/^#/;
			      chomp;
			      my @g = split;
			      push @profile3, \@g;
			  }
			  $isocore->Markov_Transition_probability_matrix(\@profile3);
		      } elsif (/maximum_number_of_donors_per_acceptor_site/i) {
			  $_ = <P>;
			  chomp;
			  s/\s//g;
			  $isocore->maximum_number_of_donors_per_acceptor_site($_);
			  $maxdonors = $_;
		      } 
		  }
		  push @isocores, $isocore;
				#print STDERR "Read isochore\n";
	      }
	  } elsif ($line =~ /General_Gene_Model/i) {
	      while (my $gmline = <P>) {
		  $gm .= $gmline;
	      }	
	      if ($gm !~ /\n$/) {
		  $gm.="\n";
	      }
	  } 
      }
      close P;	
      $self->{_isocores} = \@isocores;
      $self->{_noScore} = $NOSCORE;
      $self->{_bkgdSubFlankLen} = $BKGD_SUB_FLANK_LEN;
      $self->{_geneModel}->string_to_model($gm);
  }


sub writeParam {
    my ( $self, $filename ) = @_;
    open OUT, ">-";
    if ($filename) {
	open OUT, ">$filename";
    }
    print OUT $self->{_head},$label_comments,$blank_line;
    print OUT $label_noscore,"NO_SCORE\n",$self->{_noScore},$blank_line,$blank_line;
    print OUT "$header_bkgd_sub_flank_len",$self->{_bkgdSubFlankLen},$blank_line,$blank_line;
    my $numisocores = 0;
    if ($numisocores = $self->{_numIsocores}) {
	print OUT $label_numiso,"number_of_isochores\n",$numisocores,$blank_line,$blank_line;
    } else {
	die "undefined: number_of_isochores\n";
    }
	 
    for (my $i = 0;$i<$numisocores;$i++) {
	if (!defined  $self->{_isocores}->[$i]) {
	    my $iso_string = $i + 1;
	    print STDERR "Isochore $iso_string not defined: skipping...\n";
	    next;
	}
	my $val = 0;
	$label_isohead =~ s/\d+/$i+1/e;
	print OUT $label_isohead,$blank_line;
	if ($val = $self->{_isocores}->[$i]->boundaries_of_isochore) {
	    print OUT $label_gc,"boundaries_of_isochore\n",join(" ",@$val),$blank_line,$blank_line;
	} else {
	    die "undefined: boundaries_of_isochore";
	}

	if ($val = $self->{_isocores}->[$i]->Absolute_cutoff_exons) {
	    print OUT $label_exon_score_thresh,"Absolute_cutoff_exons\n",join(" ",@$val),$blank_line,$blank_line;
	} else {
	    die "undefined: Absolute_cutoff_exons";
	}
	if ($val = $self->{_isocores}->[$i]->Coding_cutoff_oligos) {
	    print OUT "Coding_cutoff_oligos\n",join(" ",@$val),$blank_line,$blank_line;
	} else {
	    die "undefined: Coding_cutoff_oligos";
	}

	if ($val = $self->{_isocores}->[$i]->Site_factor) {
	    print OUT $label_exon_score_factors,"Site_factor\n",join(" ",@$val),$blank_line,$blank_line;
	} elsif ($val = $self->{_isocores}->[$i]->Exon_factor) {
	    print OUT $label_exon_score_factors,"Site_factor\n",join(" ",map { 1 - $_ } @$val),$blank_line,$blank_line; 
	} else {
	    die "undefined: Site_factor";
	}
	if ($val = $self->{_isocores}->[$i]->Exon_factor) {
	    print OUT "Exon_factor\n",join(" ",@$val),$blank_line,$blank_line;
	} else {
	    die "undefined: Exon_factor";
	}
	if ($val = $self->{_isocores}->[$i]->HSP_factor) {
	    print OUT "HSP_factor\n",join(" ",@$val),$blank_line,$blank_line;
	} else {
	    print OUT "HSP_factor\n","1.0 1.0 1.0 1.0 1.0",$blank_line,$blank_line;
	}
	if ($val = $self->{_isocores}->[$i]->Evidence_Factor) {
	    print OUT "Evidence_Factor\n",$val,$blank_line,$blank_line;
	} else {
	    print OUT "Evidence_Factor\n","1.0",$blank_line,$blank_line;
	}
	if ($val = $self->{_isocores}->[$i]->Evidence_Weight) {
	    print OUT "Evidence_Weight\n",$val,$blank_line,$blank_line;
	} else {
	    print OUT "Evidence_Weight\n","0",$blank_line,$blank_line;
	}
	if ($val = $self->{_isocores}->[$i]->RSS_Markov_Score) {
	    print OUT "RSS_Markov_Score\n",$val,$blank_line,$blank_line;
	} else {
	    print OUT "RSS_Markov_Score\n","5.0",$blank_line,$blank_line;
	}
	if ($val = $self->{_isocores}->[$i]->RSS_Donor_Score_Cutoff) {
	    print OUT "RSS_Donor_Score_Cutoff\n",$val,$blank_line,$blank_line;
	} else {
	    print OUT "RSS_Donor_Score_Cutoff\n","3.5",$blank_line,$blank_line;
	}
	if ($val = $self->{_isocores}->[$i]->RSS_Acceptor_Score_Cutoff) {
	    print OUT "RSS_Acceptor_Score_Cutoff\n",$val,$blank_line,$blank_line;
	} else {
	    print OUT "RSS_Acceptor_Score_Cutoff\n","3.5",$blank_line,$blank_line;
	}

	if ($val = $self->{_isocores}->[$i]->U12_Splice_Score_Threshold) {
	    print OUT "U12_Splice_Score_Threshold\n",$val,$blank_line,$blank_line;
	}
	if ($val = $self->{_isocores}->[$i]->U12_Exon_Score_Threshold) {
	    print OUT "U12_Exon_Score_Threshold\n",$val,$blank_line,$blank_line;
	}
	if ($val = $self->{_isocores}->[$i]->U12_Exon_weight) {
	    print OUT "U12_Exon_weight\n",$val,$blank_line,$blank_line;
	}

	if ($val = $self->{_isocores}->[$i]->Exon_weights) {
	    print OUT "Exon_weights\n",join(" ",@$val),$blank_line,$blank_line;
	} else {
	    die "undefined: Exon_weights";
	}

	if ($val = $self->{_isocores}->[$i]->get_profile("Start_profile")) {
	    print OUT $label_profiles,"Start_profile\n",join(" ",($val->{length},$val->{offset},$val->{cutoff},$val->{order})),$blank_line,$label_profile_array;
	    foreach my $e (@{$val->{profile}}) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	} else {
	    die "undefined: Start_profile";
	}
	if ($val = $self->{_isocores}->[$i]->get_profile("Branch_point_profile")) {
	    print OUT "Branch_point_profile\n$label_prof_props",join(" ",($val->{length},$val->{offset},$val->{cutoff},$val->{order},$val->{afactor},$val->{bfactor},$val->{acc_context},$val->{dist},$val->{optdist},$val->{penalty_factor})),$blank_line,$label_profile_array;
	    foreach my $e (@{$val->{profile}}) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	}
	if ($val = $self->{_isocores}->[$i]->get_profile("U12_Branch_point_profile")) {
	    print OUT "U12_Branch_point_profile\n$label_prof_props",join(" ",($val->{length},$val->{offset},$val->{cutoff},$val->{order},$val->{afactor},$val->{bfactor},$val->{acc_context},$val->{dist},$val->{optdist},$val->{penalty_factor})),$blank_line,$label_profile_array;
	    foreach my $e (@{$val->{profile}}) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	}
	if ($val = $self->{_isocores}->[$i]->get_profile("Poly_Pyrimidine_Tract_profile")) {
	    print OUT "Poly_Pyrimidine_Tract_profile\n$label_prof_props",join(" ",($val->{length},$val->{offset},$val->{cutoff},$val->{order},$val->{afactor},$val->{bfactor},$val->{acc_context},$val->{dist},$val->{optdist},$val->{penalty_factor})),$blank_line,$label_profile_array;
	    foreach my $e (@{$val->{profile}}) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	} 
	if ($val = $self->{_isocores}->[$i]->get_profile("U12gtag_Acceptor_profile")) {
	    print OUT "U12gtag_Acceptor_profile\n",join(" ",($val->{length},$val->{offset},$val->{cutoff},$val->{order},$val->{afactor},$val->{bfactor})),$blank_line,$label_profile_array;
	    foreach my $e (@{$val->{profile}}) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	} 
	if ($val = $self->{_isocores}->[$i]->get_profile("U12atac_Acceptor_profile")) {
	    print OUT "U12atac_Acceptor_profile\n",join(" ",($val->{length},$val->{offset},$val->{cutoff},$val->{order},$val->{afactor},$val->{bfactor})),$blank_line,$label_profile_array;
	    foreach my $e (@{$val->{profile}}) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	} 
	if ($val = $self->{_isocores}->[$i]->get_profile("Acceptor_profile")) {
	    print OUT "Acceptor_profile\n",join(" ",($val->{length},$val->{offset},$val->{cutoff},$val->{order},$val->{afactor},$val->{bfactor})),$blank_line,$label_profile_array;
	    foreach my $e (@{$val->{profile}}) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	} else {
	    die "undefined: Acceptor_profile";
	}
	if ($val = $self->{_isocores}->[$i]->get_profile("U12gtag_Donor_profile")) {
	    print OUT "U12gtag_Donor_profile\n",join(" ",($val->{length},$val->{offset},$val->{cutoff},$val->{order},$val->{afactor},$val->{bfactor})),$blank_line,$label_profile_array;
	    foreach my $e (@{$val->{profile}}) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	}
	if ($val = $self->{_isocores}->[$i]->get_profile("U12atac_Donor_profile")) {
	    print OUT "U12atac_Donor_profile\n",join(" ",($val->{length},$val->{offset},$val->{cutoff},$val->{order},$val->{afactor},$val->{bfactor})),$blank_line,$label_profile_array;
	    foreach my $e (@{$val->{profile}}) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	}
	if ($val = $self->{_isocores}->[$i]->get_profile("U2gcag_Donor_profile")) {
	    print OUT "U2gcag_Donor_profile\n",join(" ",($val->{length},$val->{offset},$val->{cutoff},$val->{order},$val->{afactor},$val->{bfactor})),$blank_line,$label_profile_array;
	    foreach my $e (@{$val->{profile}}) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	}
	if ($val = $self->{_isocores}->[$i]->get_profile("U2gta_Donor_profile")) {
	    print OUT "U2gta_Donor_profile\n",join(" ",($val->{length},$val->{offset},$val->{cutoff},$val->{order},$val->{afactor},$val->{bfactor})),$blank_line,$label_profile_array;
	    foreach my $e (@{$val->{profile}}) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	}
	if ($val = $self->{_isocores}->[$i]->get_profile("U2gtg_Donor_profile")) {
	    print OUT "U2gtg_Donor_profile\n",join(" ",($val->{length},$val->{offset},$val->{cutoff},$val->{order},$val->{afactor},$val->{bfactor})),$blank_line,$label_profile_array;
	    foreach my $e (@{$val->{profile}}) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	}
	if ($val = $self->{_isocores}->[$i]->get_profile("U2gty_Donor_profile")) {
	    print OUT "U2gty_Donor_profile\n",join(" ",($val->{length},$val->{offset},$val->{cutoff},$val->{order},$val->{afactor},$val->{bfactor})),$blank_line,$label_profile_array;
	    foreach my $e (@{$val->{profile}}) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	}
	if ($val = $self->{_isocores}->[$i]->get_profile("Donor_profile")) {
	    print OUT "Donor_profile\n",join(" ",($val->{length},$val->{offset},$val->{cutoff},$val->{order},$val->{afactor},$val->{bfactor})),$blank_line,$label_profile_array;
	    foreach my $e (@{$val->{profile}}) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	} else {
	    die "undefined: Donor_profile";
	}
	if ($val = $self->{_isocores}->[$i]->get_profile("Stop_profile")) {
	    print OUT "Stop_profile\n",join(" ",($val->{length},$val->{offset},$val->{cutoff},$val->{order})),$blank_line,$label_profile_array;
	    foreach my $e (@{$val->{profile}}) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	} else {
	    die "undefined: Stop_profile";
	}
	if ($val = $self->{_isocores}->[$i]->get_profile("PolyA_Signal_profile")) {
	    print OUT "PolyA_Signal_profile\n",join(" ",($val->{length},$val->{offset},$val->{cutoff},$val->{order})),$blank_line,$label_profile_array;
	    foreach my $e (@{$val->{profile}}) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	} 
	if ($val = $self->{_isocores}->[$i]->Markov_order) {
	    print OUT $label_coding_potential,$label_initial_prob,"Markov_order\n",$val,$blank_line;
	} else {
	    die "undefined: Markov_order\n";
	}
	if ($val = $self->{_isocores}->[$i]->{_Markov_Initial_probability_matrix}) {
	    print OUT "Markov_Initial_probability_matrix\n";
	    foreach my $e (@$val) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	} else {
	    die "undefined:Markov_Initial_probability_matrix\n";
	}
	if ($val = $self->{_isocores}->[$i]->Markov_Transition_probability_matrix) {
	    print OUT $label_transition_prob,"Markov_Transition_probability_matrix\n";
	    foreach my $e (@$val) {
		print OUT join(" ",@$e),$blank_line;
	    }
	    print OUT $blank_line;
	} else {
	    die "undefined: Markov_Transition_probability_matrix\n";
	}

	if ($val = $self->{_isocores}->[$i]->maximum_number_of_donors_per_acceptor_site) {
	    print OUT $label_donors_per_acceptor,"maximum_number_of_donors_per_acceptor_site\n",$val,$blank_line,$blank_line;
	} else {
	    die "undefined: $label_donors_per_acceptor";
	}
    }
    if (defined $self->{_geneModel}) {
	print OUT $label_gene_model,"General_Gene_Model\n",$self->{_geneModel}->modelString;
    } else {
	die "undefined: $label_gene_model\n";
    }
    close OUT;
}

sub noScore {
    my ( $self, $noscore ) = @_;
    $self->{_noScore} = $noscore if defined($noscore);
    return $self->{_noScore};
}
sub bgkdSubtractFlankLength {
    my ( $self, $bsfl ) = @_;
    $self->{_bkgdSubFlankLen} = $bsfl if defined($bsfl);
    return $self->{_bkgdSubFlankLen};
}
sub geneModel {
    my ( $self, $gm ) = @_;
    $self->{_geneModel} = $gm if defined($gm);
    return $self->{_geneModel};
}
sub isocores {
    my ( $self, $iso ) = @_; #iso is a ref to array of GeneID::Isocore objects
    $self->{_isocores} = $iso if defined($iso);
    return $self->{_isocores};
}
sub numIsocores {
    my ( $self, $numiso ) = @_; 
    $self->{_numIsocores} = $numiso if defined($numiso);
    return $self->{_numIsocores};
}

1;
