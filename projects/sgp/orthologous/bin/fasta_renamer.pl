#!/usr/bin/perl -w
# This is perl, version 5.005_03 built for i386-linux
#line 1960 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
# $Id: fasta_renamer.pl,v 1.1 2001-10-01 14:09:27 jabril Exp $
#line 1786 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
#
use strict;
#line 1320 "/projects/datasets/orthologous/_docs/ORTHOLOGOUS_dataset.nw"
#
# fasta_renamer.pl infile descfile new_seq_id > outfile
#
#     Replacing sequence name 
#     for single sequence fasta files
#     and reformating sequence to 80 cols.
#
use lib qw( /usr/lib/perl5/site_perl/5.005/ /usr/lib/perl5/5.00503/ ) ;
use Bio::Seq;
use Bio::SeqIO;

my ($infile,$descfile,$newid) = @ARGV;

my $seqin  = Bio::SeqIO->new(-format => 'FASTA', -file => "$infile");
my $seqout = Bio::SeqIO->new(-format => 'FASTA', -fh => \*STDOUT);

open(DESC,"> $descfile");
while (my $sequence = $seqin->next_seq()) {
    my ($sid,$len,$seq,$desc);
    print STDERR "### READING FASTA... $newid\n";
    $sid  = $sequence->display_id();
    $desc = $sequence->desc();
    $len  = $sequence->length();
    $seq  = $sequence->seq();
    $sid =~ s/\s+/\_/og;
    $seq =~ tr/[a-z]/[A-Z]/;
    print STDERR "### WRITING FASTA... $newid\n";
    print DESC "$newid $sid $len $desc\n";
    $sequence->display_id($newid);
    $sequence->desc('');
    $sequence->seq($seq);
    $seqout->write_seq($sequence);
}; # while 
close(DESC);
exit(0);
