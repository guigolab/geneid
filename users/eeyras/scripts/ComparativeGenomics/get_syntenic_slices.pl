#!/usr/local/bin/perl -w

use strict;  
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

# compara_db
my $compara_dbname = 'ensembl_compara_25_1';
my $compara_dbhost = 'ensembldb.sanger.ac.uk';
my $user = "anonymous";
my $compara_config =  '';

my $compara_db = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
							      -user      => $user,
							      -dbname    => $compara_dbname,
							      -host      => $compara_dbhost,
							      -conf_file => $compara_config,
							      );
my $adaptor = $compara_db->get_DnaAlignFeatureAdaptor;


my $chr_name = "X";
my $chr_start = 1;
my $chr_end   = 10100000;

my $target_species = "Mus musculus";
my $focus_species  = "Homo sapiens";

my $target_assembly_type = "NCBIM33";
my $focus_assembly_type  = "NCBI34";



my @features = @{$adaptor->fetch_all_by_species_region($focus_species,
						       $focus_assembly_type,
						       $target_species,
						       $target_assembly_type,
						       $chr_name,
						       $chr_start, 
						       $chr_end,
						       'WGA'
						       )};

print scalar(@features)." features\n";
foreach my $f (@features){
    print $f->gff_string."\n";
}
