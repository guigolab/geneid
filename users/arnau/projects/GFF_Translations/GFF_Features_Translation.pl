#!/usr/local/bin/perl -w

# To do
# * Parse the gene systematic identifier
# * GFF Parsing in a different package

# Issue warnings about suspicious programming.
use warnings 'all';

# Must declare and initialize all variables
use strict;

# be prepare for command-line options/arguments
use Getopt::Std;

use Data::Dumper;

# Bioperl
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::Location::Simple;
use Bio::Location::Split;

sub help {
return <<"END_HELP";
Description: Translate Features in GFF format
Usage:

    promoter_extraction.pl [-h] -s {Sequence in FASTA format} -f {Features in GFF format} -t {translation table}
	-h help
	-s Sequence in FASTA format
	-f Features in GFF format
	-t Translation table (Standard => 1, Bacterial => 11) - default is the Standard Genetic Code
	
Examples using some combinations:
	GFF_Features_Translation.pl -s AC005155.fa -f AC005155.geneid.gff -t 1

END_HELP

}

BEGIN {

    use vars qw /%opts/;
    
    # these are switches taking an argument (a value)
    my $switches = 'hf:s:t:';
    
    # Get the switches
    getopts($switches, \%opts);
    
    # If the user does not write nothing, skip to help
    if (defined($opts{h}) || !defined($opts{f}) || !defined($opts{s})){
	print STDERR help;
	exit 0;
    }
    
}

my $_debug = 1;

my $seqfile;
my $featurefile;
my $translation_table = 1;

defined $opts{f} and $featurefile       = $opts{f};
defined $opts{s} and $seqfile           = $opts{s};
defined $opts{t} and $translation_table = $opts{t};

# Check that the files exist !!!

if (not -f $featurefile) {
    print STDERR "Error, can't find feature file, $featurefile\n";
    exit 0;
}
if (not -f $seqfile) {
    print STDERR "Error, can't find sequence file, $seqfile\n";
    exit 0;
}

if ((not ($translation_table =~ /^\d+$/)) || ($translation_table > 22)) {
   print STDERR "the translation code, $translation_table, is not recognized !\n";
   print STDERR "must be a number between 1 and 21\n";
   exit 0;
}

# The sequence object factory

my $sin  = Bio::SeqIO->new(
			   -file   => $seqfile,
			   -format => 'fasta'
			   );

# The peptide sequence factory

my $sout = Bio::SeqIO->new(
			   -format => 'fasta'
			   );

if ($_debug) {
    print STDERR "processing the GFF file...\n";
}
my $featuresPerSeqId = parse_gff ($featurefile);

if ($_debug) {
    print STDERR "processing done.\n";
    print STDERR "processing the sequence fasta file\n";
}

while ( my $seqobj = $sin->next_seq() ) {
    if ($_debug) {
	print STDERR "processing sequence ",$seqobj->id," first 10 bases ",$seqobj->subseq(1,10),"\n";
    }
    
    my $features = $featuresPerSeqId->{$seqobj->id};
    
    if ((defined $features) && (@$features > 0)) {
	foreach my $f (@$features) {
	    # The gene name
	    my $geneName = $f->display_name;

	    $f->attach_seq ($seqobj);
	    my $splicedseq = $f->spliced_seq();
	    my $pepobj     = $splicedseq->translate(undef, undef, undef, $translation_table);
	    if (defined $pepobj) {
		
		# Changing the id to match the gene name

		if (defined $geneName) {
		    if ($_debug) {
			print STDERR "changing peptide identifier to gene name, $geneName\n";
		    }

		    $pepobj->id ($geneName);

		}

		if ($_debug) {
		    print STDERR "writing out translation for " . $pepobj->id . "\n";
		}

		$sout->write_seq($pepobj);
	    }
	    else {
		print STDERR "Error, can't translate feature!\n";
	    }
	}
	
    }
    else {
	if ($_debug) {
	    print STDERR "No features for sequence, " . $seqobj->id . "\n";
	} 
    }
}

if ($_debug) {
    print STDERR "processing done.\n";
}

##
# End
##


sub parse_gff {
    my ($featurefile) = @_;
    my $features;
    
    open GENEID, "$featurefile" or die "can't open gff file, $featurefile!\n";

    my @features = ();
    my $seqName = undef;
    
    while (<GENEID>) {
	next if /^##/;
	my $line = $_;
	chomp $line;
	
	if ($_debug) {
	    print STDERR "new line : $line\n";
	}
	
	# extracting the number of exons for the current gene
	
	if ($line =~ /\# Gene [^\s]+\s+\(\w+\)\. (\d+) exons\..+/) {
	    
	    my $nb_exons = $1;
	    my $i = 1;
	    
	    if ($_debug) {
		print STDERR "new gene of $nb_exons exons found\n";
	    }
	    
	    # CDS locations
	    my $splitlocation = new Bio::Location::Split();
	    my $codon_start   = undef;
	    my $geneName      = undef;

	    # build the cds feature by joining the exon together
	    
	    while ($i <= $nb_exons) {
		$line = <GENEID>;
		chomp $line;

		if ($_debug) {
		    print STDERR "parsing line, $line...\n";
		}

		# remove the white spaces just keep the tabulations
		$line =~ s/ //g;
		
		my $featureName = extractFeatureName ($line);
	
		if ($_debug) {
		    print STDERR "featureName, $featureName\n";
		}
	
		# Make sure we are parsing an exon feature
		if ((defined $featureName) && ($featureName =~ /exon|single|first|internal|terminal/i)) {
		    
		    if ($line =~ /^([^\t]+)\t+(\S+)\t+(\w+)\t+(\d+)\t+(\d+)\t+([^\t]*)\t+(\+|-)\t+(\d*)\t+(.+)/) {
			if ($_debug) {
			    print STDERR "feature line : $line\n";
			}
			
			my $seqName_tmp = $1;
			my $start       = $4;
			my $end         = $5;
			my $strand      = 1;
			my $strand_tmp  = $7;
			if ($strand_tmp =~ /-/) {
			    $strand = -1;
			}
			my $frame       = $8;
			$geneName = $9;
			if (defined $geneName) {
			    chomp $geneName;
			}
			
			# parsing problem when the strand is "-", can not get the frame in that case !!!
			# so substitute - by + to make the parsing work
			if (not defined $frame) {
			    $line =~ s/-/\+/g;
			    $line =~ /(\w+)\t+(\S+)\t+(\w+)\t+(\d+)\t+(\d+)\t+(\S+)\t+(\+|-)\t+(\d+)\t+(.+)/;
			    $frame = $8;
			    
			    $geneName = $9;
			    chomp $geneName;
			}
			# specify the frame for the first exon, which might not start with a MET
			if (($strand==1 && $i==1) || ($strand==-1 && $i==$nb_exons)) {
			    if ($frame =~ /0/) {
				$codon_start = 1;
			    }
			    elsif ($frame =~ /1/) {
				$codon_start = 2;
			    }
			    elsif ($frame =~ /2/) {
				$codon_start = 3;
			    }
			}
			    
			if ($_debug) {
			    print STDERR "found gene name, $geneName\n";
			}
			    
			if (not defined $seqName) {
			    if ($_debug) {
				print STDERR "Assigning the sequence name, $seqName\n";
			    }
			    $seqName = $seqName_tmp;
			}
			elsif ($seqName ne $seqName_tmp) {
			    
			    print STDERR "Adding " . @features . " features to sequence, $seqName\n";
			    
			    $features->{$seqName} = [ @features ];
			    @features = ();
			    $seqName  = $seqName_tmp;
			}
			
			my $sublocation = new Bio::Location::Simple (
								     -start  => $start,
								     -end    => $end,
								     -strand => $strand
								     );
			$splitlocation->add_sub_Location ($sublocation);
			
			$i++;
		    }
		    else {
			print STDERR "can't parsed current GFF feature line!\n";
		    }
		}
		else {
		    print STDERR "not an exon line\n";
		}
	    } # While
	    
	    # build the feature object

	    if ($_debug) {
		print STDERR "instanciating feature object...\n";
	    }

	    my $feature = new Bio::SeqFeature::Generic (
							-location     => $splitlocation,
							-primary      => 'CDS',
							-display_name => "$geneName",
							);
	    $feature->add_tag_value ('codon_start', $codon_start);

	    if ($_debug) {
		print STDERR "instanciation done.\n";
	    }
	    
	    push (@features, $feature);
	}
    }

    if ($_debug) {
	print STDERR "Adding " . @features . " features to sequence, $seqName\n";
    }
    
    $features->{$seqName} = [ @features ];
    
    close GENEID;
    
    return $features;
    
}

sub extractFeatureName {
    my ($gff_line) = @_;

    if ($gff_line =~ /^[^\t]+\t+[^\t]+\t+([^\t]+)\t+.+/) {
	my $featureName = $1;
	return $featureName;
    }
    
    return undef;

}
