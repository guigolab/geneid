package Geneid::GeneModel;
use strict;
# use vars qw(@ISA);

############################################################
# 
# Author: T. Alioto
# Date: 090422
# POD documentation - main docs before the code

=head1 NAME

Geneid::GeneModel - Handler for GeneModel objects in GeneID parameter files

=head1 SYNOPSIS

    use Geneid::Param;
    use Geneid::Isocore;
    use Geneid::GeneModel;

    my ($profile,$length,$newoff) = getKmatrix($true_seqs,$false_seqs,$order,$offset);
    my $param = Geneid::Param->new(); #make a new parameter object
    $param->readParam($pin); #read in a template parameter file

    for (my $i = 0;$i < $param->numIsocores ; $i++){
	    if (!defined @{$param->isocores}[$i]->set_profile($label,$length,$newoff,$cutoff,$order,$profile)){die "error in setting profile\n";} #set a new profile where $profile is the profile name.
    }  
    $param->geneModel = Geneid::GeneModel->new();




=head1 DESCRIPTION

Geneid::GeneModel is a handler module for the Gene Model in a Geneid 
parameter file.  


=head1 OBJECT METHODS

See below for more detailed summaries. 


=head1 AUTHOR - Tyler Alioto

Email: tyler.alioto@crg.es

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


sub new
  {
      my ($class, $m) = @_;
      my $self = {
		  _model_string => undef,
		  _model => {
			     IntragenicRules => undef,
			     ExternalRules   => undef,
			     IntergenicRules => undef,
			     BeginEndRules   => undef
			     }
		 };
      
      bless $self, $class;
      if (defined $m){
	  string_to_model($m);
      }
      return $self;
  }

sub eraseModel{
    my ( $self, $val ) = @_;
    $self->{_model}->{IntragenicRules} = undef;
    $self->{_model}->{ExternalRules} = undef;
    $self->{_model}->{IntergenicRules} = undef;
    $self->{_model}->{BeginEndRules} = undef;
}
sub addIntragenicRule{
    my ($self, $uc, $de, $minD, $maxD, $block) = @_;
    push @{$self->{_model}->{IntragenicRules}},{uc=>$uc, de=>$de, md=>$minD, Md=>$maxD, block=>$block};
}
sub intronRange{
    my ($self, $min, $max) = @_;
    if (defined $self->{_model}->{IntragenicRules}) {
	foreach my $rule (@{$self->{_model}->{IntragenicRules}}) {
	    $rule->{md}=$min;
	    $rule->{Md}=$max;
	}
    }
}
sub externalRange{
    my ($self, $min, $max) = @_;
    if (defined $self->{_model}->{ExternalRules}) {
	foreach my $rule (@{$self->{_model}->{ExternalRules}}) {
	    $rule->{md}=$min;
	    $rule->{Md}=$max;
	}
    }
}
sub intergenicRange{
    my ($self, $min, $max) = @_;
    if (defined $self->{_model}->{IntergenicRules}) {
	foreach my $rule (@{$self->{_model}->{IntergenicRules}}) {
	    $rule->{md}=$min;
	    $rule->{Md}=$max;
	}
    }
}

sub addExternalRule{
    my ($self, $uc, $de, $minD, $maxD, $block) = @_;
    push @{$self->{_model}->{ExternalRules}},{uc=>$uc, de=>$de, md=>$minD, Md=>$maxD, block=>$block};

}

sub addIntergenicRule{
    my ($self, $uc, $de, $minD, $maxD, $block) = @_;
    push @{$self->{_model}->{IntergenicRules}},{uc=>$uc, de=>$de, md=>$minD, Md=>$maxD, block=>$block};

}

sub addBeginEndRule{
    my ($self, $uc, $de, $minD, $maxD, $block) = @_;
    push @{$self->{_model}->{BeginEndRules}},{uc=>$uc, de=>$de, md=>$minD, Md=>$maxD, block=>$block};

}
sub modelString{
    my ( $self, $val ) = @_;
    if (defined($val)){
	$self->{_model_string} = $val;
    }else{
	$self->model_to_string();
    }
    return $self->{_model_string};
}
sub string_to_model{
    my ( $self, $ms ) = @_;
    my $type = 'intra';
    for (split /^/, $ms) {
	#print STDERR "$_\n\n";
	next if $_!~/\w/;
	chomp;
	if (m/^#.*INTRA/i){
	    $type = 'intra';
	    #print STDERR "here\n";
	}elsif(m/^#.*External/i){
	    $type = 'ext';
	}elsif(m/^#.*INTER/i){
	    $type = 'inter';
	}elsif(m/^#.*END/i){
	    $type = 'begin_end';
	}else{
	    #print STDERR "$type\n";
	    if ($type eq 'intra'){
		my $block = 0;
		my @f = split;
		my @uc = split ":",$f[0];
		my @de = split ":",$f[1];
		my ($minD,$maxD) = split ":",$f[2];
		$block = 1 if exists $f[3];
		$self->addIntragenicRule([@uc],[@de],$minD,$maxD,$block);
	    }elsif($type eq 'ext'){
		my $block = 0;
		my @f = split;
		my @uc = split ":",$f[0];
		my @de = split ":",$f[1];
		my ($minD,$maxD) = split ":",$f[2];
		$block = 1 if exists $f[3];
		$self->addExternalRule([@uc],[@de],$minD,$maxD,$block);
	    }elsif($type eq 'inter'){
		my $block = 0;
		my @f = split;
		my @uc = split ":",$f[0];
		my @de = split ":",$f[1];
		my ($minD,$maxD) = split ":",$f[2];
		$block = 1 if exists $f[3];
		$self->addIntergenicRule([@uc],[@de],$minD,$maxD,$block);
	    }elsif($type eq 'begin_end'){
		my $block = 0;
		my @f = split;
		my @uc = split ":",$f[0];
		my @de = split ":",$f[1];
		my ($minD,$maxD) = split ":",$f[2];
		$block = 1 if exists $f[3];
		$self->addBeginEndRule([@uc],[@de],$minD,$maxD,$block);
	    }
	}
    }
}
sub model_to_string{
    my ( $self, $val ) = @_;
    my $ms = "";
    if (defined $self->{_model}->{IntragenicRules}) {
	#print STDERR "print INTRAgenic\n";
	$ms .= "# INTRAgenic connections\n";
	foreach my $rule (@{$self->{_model}->{IntragenicRules}}) {
	    $ms .= sprintf "%-40s   ",join(':',@{$rule->{uc}});
	    $ms .= sprintf "%-40s   ",join(':',@{$rule->{de}});
	    $ms .= sprintf "%-20s   ",join(':',($rule->{md},$rule->{Md}));
	    $ms .= ' block' if $rule->{block};
	    $ms .= "\n";
	}
    }
    if (defined $self->{_model}->{ExternalRules}) {
	$ms .= "# External features\n";
	foreach my $rule (@{$self->{_model}->{ExternalRules}}) {
	    $ms .= sprintf "%-40s   ",join(':',@{$rule->{uc}});
	    $ms .= sprintf "%-40s   ",join(':',@{$rule->{de}});
	    $ms .= sprintf "%-20s   ",join(':',($rule->{md},$rule->{Md}));
	    $ms .= ' block' if $rule->{block};
	    $ms .= "\n";
	}
    }
    if (defined $self->{_model}->{IntergenicRules}) {
	$ms .= "# INTERgenic connections\n";
	foreach my $rule (@{$self->{_model}->{IntergenicRules}}) {
	    $ms .= sprintf "%-40s   ",join(':',@{$rule->{uc}});
	    $ms .= sprintf "%-40s   ",join(':',@{$rule->{de}});
	    $ms .= sprintf "%-20s   ",join(':',($rule->{md},$rule->{Md}));
	    $ms .= ' block' if $rule->{block};
	    $ms .= "\n";
	}
    }
    if (defined $self->{_model}->{BeginEndRules}) {
	$ms .= "# BEGINNING and END of prediction\n";
	foreach my $rule (@{$self->{_model}->{BeginEndRules}}) {
	    $ms .= sprintf "%-40s   ",join(':',@{$rule->{uc}});
	    $ms .= sprintf "%-40s   ",join(':',@{$rule->{de}});
	    $ms .= sprintf "%-20s   ",join(':',($rule->{md},$rule->{Md}));
	    $ms .= ' block' if $rule->{block};
	    $ms .= "\n";
	}
    }
    $self->{_model_string} = $ms;
}
sub useDefault{
    my ($self) = @_;
    my $block = 0;

    my $minD = 20;
    my $maxD = 40000;
    my @uc = qw(First+ Internal+);
    my @de = qw(Internal+ Terminal+);
    $self->addIntragenicRule([@uc],[@de],$minD,$maxD,$block);
    @uc = qw(Internal- Terminal-);
    @de = qw(First- Internal-);
    $self->addIntragenicRule([@uc],[@de],$minD,$maxD,$block);

    $minD = 50;
    $maxD = 4000;
    @uc = qw(Promoter+);
    @de = qw(First+ Single+);
    $self->addExternalRule([@uc],[@de],$minD,$maxD,$block);
    @uc = qw(First- Single-);
    @de = qw(Promoter-);
    $self->addExternalRule([@uc],[@de],$minD,$maxD,$block);
    @uc = qw(Terminal+ Single+);
    @de = qw(aataaa+);
    $self->addExternalRule([@uc],[@de],$minD,$maxD,$block);
    @uc = qw(aataaa-);
    @de = qw(Terminal- Single-);
    $self->addExternalRule([@uc],[@de],$minD,$maxD,$block);

    $minD = 500;
    $maxD = 'Infinity';
    @uc = qw(aataaa+ Terminal+ Single+);
    @de = qw(Single+ First+ Promoter+);
    $self->addIntergenicRule([@uc],[@de],$minD,$maxD,$block);
    @uc = qw(aataaa+ Terminal+ Single+);
    @de = qw(aataaa- Terminal- Single-);
    $self->addIntergenicRule([@uc],[@de],$minD,$maxD,$block);
    @uc = qw(Single- First- Promoter-);
    @de = qw(aataaa- Terminal- Single-);
    $self->addIntergenicRule([@uc],[@de],$minD,$maxD,$block);
    @uc = qw(Single- First- Promoter-);
    @de = qw(Single+ First+ Promoter+);
    $self->addIntergenicRule([@uc],[@de],$minD,$maxD,$block);

    $minD = 0;
    $maxD = 'Infinity';
    @uc = qw(Begin+);
    @de = qw(First+ Internal+ Terminal+ Single+);
    $self->addBeginEndRule([@uc],[@de],$minD,$maxD,$block);
    @uc = qw(Begin-);
    @de = qw(First- Internal- Terminal- Single-);
    $self->addBeginEndRule([@uc],[@de],$minD,$maxD,$block);
    @uc = qw(First+ Internal+ Terminal+ Single+);
    @de = qw(End+);
    $self->addBeginEndRule([@uc],[@de],$minD,$maxD,$block);
    @uc = qw(First- Internal- Terminal- Single-);
    @de = qw(End-);
    $self->addBeginEndRule([@uc],[@de],$minD,$maxD,$block);

}
1;
