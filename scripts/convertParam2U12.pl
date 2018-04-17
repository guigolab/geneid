#!/usr/bin/env perl
use strict;
use Getopt::Long;
use Geneid::Param;
use Geneid::Isocore;

my @U12atts = qw(RSS_Markov_Score RSS_Donor_Score_Cutoff RSS_Acceptor_Score_Cutoff U12_Splice_Score_Threshold U12_Exon_Score_Threshold U12_Exon_weight);
my @U12profiles = qw(U12_Branch_point_profile U12gtag_Acceptor_profile U12atac_Acceptor_profile U12gtag_Donor_profile U12atac_Donor_profile);# Poly_Pyrimidine_Tract_profile
my $u12file = 0;
#my $parampath = "~/projects/u12/geneid_training/param/";

my $pin = "-";

GetOptions(
		'u12param:s'		=> \$u12file,
		'pin|paramin:s'	=> \$pin,
	   );
my $usage = "Usage: $0 -pin <paramfile_in> -u12param <u12_paramfile>\n";
print $usage and exit unless ($u12file);
my $param = Geneid::Param->new();
my $u12param = Geneid::Param->new();
$param->readParam($pin);
$u12param->readParam($u12file);
my $u12iso = @{$u12param->isocores}[0];

for (my $i = 0;$i < $param->numIsocores ; $i++){
	#my $targetiso = @{$param->isocores}[$i];
	foreach my $attribute (@U12atts){
		@{$param->isocores}[$i]->{"_".$attribute} = $u12iso->{"_".$attribute};
	}
	foreach my $u12profile (@U12profiles){
		@{$param->isocores}[$i]->copy_profile($u12iso->get_profile($u12profile)); 
	}
}   
#$param->geneModel($u12param->geneModel);
$param->writeParam;
