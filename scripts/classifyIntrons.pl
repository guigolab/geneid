#!/usr/bin/env perl

#author: Tyler Alioto (tyler.alioto@gmail.com)
#former script name: exon2intronGFFflankseqScoreBDF2.pl

#make sure the required modules are in your PERL5LIB path
#if you don't have one set modify and uncomment the following line to include the location of the Bioperl modules and the SeqOp module distributed with geneid
#use lib "$ENV{HOME}/myperlmods/";

use Bio::Tools::GFF;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Seq;
use Bio::FeatureIO;
#use Bio::FeatureIO::gff;
use Bio::SeqFeature::Annotated;
use Bio::SeqFeature::Gene::Transcript;
use Bio::SeqFeature::Gene::Exon;
use Data::Dumper;
use strict;
use Bio::DB::Fasta;
use SeqOp;
use File::Temp qw/ tempfile tempdir /;
my $tmpdir = $ENV{TMPDIR};
my $SCORE = 0;
my $GROUP = '';
my $DONOR = '';
my $FRAME = ".";
my $donor = 0;
my $STRAND = "+";
my $SEQNAME = "";
#my @gff = ();
my $u12isdonor = 0;

my $version = 3; #### IMPORTANT: must be gff3 due to bug in Bio::Tools::GFF
my $oversion = 3;
my $userscore = undef;
my $gmapdb = 'hg18';
my $lcsid = 0;
my $flank_length = 100;		# tx seq
my $donss_length = 15;		# genomic seq
my $accss_length = 53;		# genomic seq
my $ftype = "exon";
my $scoreIntrons = 1;
my $dn = 1;
my $suf = ".fa";
my $introngff;
my $fafile = 0;
my $param = "~/param/human.101007.scoring.param";
my $exec ="geneid";
my $in = "-";
my $species = 0;
my $fa = "";
GetOptions(
	   'i|input:s'   => \$in,
	   'v:s'       => \$version,
	   'ov:s'      =>$oversion,
	   'us:s'      => \$userscore,
	   #'gmapdb|db|d:s' => \$gmapdb,
	   's|fa:s' => \$fa,
	   'Chr2chr'   => \$lcsid,
	   'f|flank:s' => \$flank_length,
	   't|type:s'  => \$ftype,
	   'score!'     => \$scoreIntrons,
	   'p|param:s' => \$param
	  );

#die "IMPORTANT: must be gff3 due to bug in Bio::Tools::GFF" if $version != 3;
my $db = Bio::DB::Fasta->new($fa);
# get chr lengths
my %chrlen;
my @ids = $db->ids;
foreach my $id (@ids){
    print STDERR "$id\n";
	$chrlen{$id}=$db->length($id);
}
print STDERR "Reading GFF input...\n";
my $gffio = Bio::Tools::GFF->new(-file => $in, -gff_version => $version); 
my $gffout = Bio::Tools::GFF->new(-fh => \*STDOUT, -gff_version => $oversion); 
my $gfferror = Bio::Tools::GFF->new(-fh => \*STDERR, -gff_version => $oversion); # loop over the input stream
my %transcripts;
while (my $feature = $gffio->next_feature()) {
    next if ($feature->primary_tag !~ /$ftype/);
    if ($lcsid && $feature->seq_id =~ /^Chr/) {
	my $sid = $feature->seq_id; $sid=~s/^Chr/chr/; $feature->seq_id($sid);
    }
    my %att;
    #print STDERR $feature->gff_string;
    my @fields = split /\t/,$feature->gff_string;
    #print STDERR $fields[8],"\n";
    foreach my $key ( $feature->get_all_tags ) {
	my @values = $feature->get_tag_values($key);
	$att{$key}=\@values;
    }
    if (exists $att{'Parent'}) {
	$att{transcript_id}=$att{Parent};
    } elsif (exists $att{transcript_id}) {
    } elsif (exists $att{group}) {
	$att{transcript_id}=$att{group};
    } elsif (exists $att{GenePrediction}) {
	$att{transcript_id}=$att{GenePrediction};
    } elsif (exists $att{Match}) {
	$att{transcript_id}=$att{Match};
    } elsif (exists $att{proteinId}) {
	$att{transcript_id}=$att{proteinId};
    } else {
	my $grp = $fields[8];
	$grp =~ s/\s+/_/g;
	$att{transcript_id}=[$grp];
	#print STDERR "$grp\n";
    }
    foreach my $p (@{$att{transcript_id}}) {
	my $newfeature = Bio::SeqFeature::Gene::Exon->new(-start=>$feature->start,
							  -end=>$feature->end,
							  -seq_id=>$feature->seq_id,
							  -source_tag=>$feature->source_tag,
							  -primary_tag=>$feature->primary_tag,
							  -strand=>$feature->strand,
							  -frame=>$feature->frame
							 );
	my $original_p = $p;
	$p = $feature->seq_id ."_". $feature->strand() ."_$p";
	foreach my $key ( $feature->get_all_tags ) {
	    my @values = $feature->get_tag_values($key);
	    $att{$key}=\@values;
	    foreach my $val (@values) {
		$newfeature->add_tag_value($key,$val);
	    }
	}
	if (exists $transcripts{$p}) {
	    $transcripts{$p}->add_exon($newfeature);
	} else {
	    my $transcript = Bio::SeqFeature::Gene::Transcript->new();
	    $transcript->add_exon($newfeature);
	    $transcript->display_name($original_p);
	    $transcripts{$p}=$transcript;
	}
    }
}
$gffio->close;

print STDERR scalar (keys %transcripts), " transcripts found\n";

foreach my $t (keys %transcripts) {
    my @lastGRPgff;
    my $exon_number = 0;
    my $DONOR = 0;
    my $donor = 0;
    my $SCORE = 0;
    my $GROUP = '';
    my $FRAME = ".";
    my $STRAND = "+";
    my $SEQNAME = "";
    my $tseq = get_tseq($transcripts{$t},$db);
    my $prevlen = 0;
    foreach my $feature (sort {$a->start <=> $b->start} $transcripts{$t}->exons) {
	$exon_number++;
	my %att;
	foreach my $key ( $feature->get_all_tags) {
	    my @values = $feature->get_tag_values($key);
	    $att{$key}=\@values;
	}
	if (exists $att{'Parent'}) {
	    $att{transcript_id}=$att{Parent};
	} elsif (exists $att{transcript_id}) {
	} elsif (exists $att{group}) {
	    $att{transcript_id}=$att{group};
	} elsif (exists $att{GenePrediction}) {
	    $att{transcript_id}=$att{GenePrediction};
	} elsif (exists $att{Match}) {
	    $att{transcript_id}=$att{Match};
	} else {
	}
	if ($exon_number == 1) {
	    $GROUP = $att{transcript_id}->[0];
	    $STRAND = $feature->strand;
	    $SEQNAME = $feature->seq_id;

	    $DONOR=$feature->end +1;
	    $SCORE = $feature->score if defined $feature->score;
	} else {
	    my $newintron = Bio::SeqFeature::Generic->new(-start=>$DONOR,
							  -end=>$feature->start -1,
							  -seq_id=>$feature->seq_id,
							  -source_tag=>$feature->source_tag,
							  -primary_tag=>"Intron",
							  -strand=>$feature->strand,
							  -score=>$feature->score,
							  -frame=>$feature->frame,
							 );
	    #print STDERR $feature->frame,"\n";
	    foreach my $annot (keys %att) {
		foreach my $val (@{$att{$annot}}) {
		    if ($annot eq "Target") {
			$annot="target";
			my ($n,$st,$en,$tstrand)= split /\s+/,$val;
			if ((($feature->strand > 0)&&($tstrand eq "+"))||(($feature->strand < 0)&&($tstrand eq "-"))) {
			    $en = $st; $st = $st-1;
			} else {
			    $st = $en; $en = $st+1;
			}
			$val = "$n $st $en $tstrand";
		    }
		    $newintron->add_tag_value($annot,$val);
		}
	    }
	    $newintron->score(($SCORE + $feature->score)/2) if defined $feature->score;
	    $SCORE = $newintron->score if defined $newintron->score;
	    if ($feature->frame ne ".") {
		if ($newintron->strand == 1) {
		    $newintron->frame((3 - $feature->frame)%3 );
		} else {
		    $newintron->frame((($feature->end - $feature->start + 1) - $feature->frame)%3 );
		}
	    }
	    $newintron->score($userscore) if ! defined $newintron->score  && defined $userscore;
	    if (($newintron->end < $newintron->start)||( ($newintron->end - $newintron->start + 1) > 5000000 )||( ($newintron->end - $newintron->start + 1) < 30 )) {
		print STDERR ('Warning: abnormal intron size (',($newintron->end - $newintron->start + 1),")\n");
		$newintron->add_tag_value('bad_length',$newintron->end - $newintron->start + 1);
		$gfferror->write_feature($newintron);
	    } 
		#add flanking seqs to tags
		my $tx_len = length($tseq);
		my $pos = $newintron->strand > 0 ? $prevlen : $tx_len - $prevlen ;
		my $don_flank = "";
		my $acc_flank = "";
		my $padlength = 0;
		if ($pos - $flank_length < 0) {
		    $padlength = $flank_length - $pos;
		    $don_flank = 'N' x $padlength;
		    $don_flank .= substr($tseq,0, $flank_length - $padlength);
		} else {
		    $don_flank = substr($tseq,$pos - $flank_length, $flank_length);
		}
		if ($pos + $flank_length > $tx_len) {
		    $padlength = $flank_length + $pos - $tx_len;
		    $acc_flank = substr($tseq,$pos, $flank_length - $padlength);
		    $acc_flank .= 'N' x $padlength;
		} else {
		    $acc_flank = substr($tseq,$pos, $flank_length);
		}
		my $acc_seq = "";
		my $don_seq = "";
		if ($newintron->strand > 0) {
		    $don_seq = SeqOp::get_seq_BioDBFasta($db, $newintron->seq_id,$newintron->start, $newintron->start + 14, "+");
		    $acc_seq = SeqOp::get_seq_BioDBFasta($db, $newintron->seq_id,$newintron->end - 52, $newintron->end, "+");
		} else {
		    $acc_seq = SeqOp::get_seq_BioDBFasta($db, $newintron->seq_id,$newintron->start, $newintron->start + 52, "-");
		    $don_seq = SeqOp::get_seq_BioDBFasta($db, $newintron->seq_id,$newintron->end - 14, $newintron->end, "-");
		}
		$newintron->add_tag_value('acc_flank',$acc_flank);
		$newintron->add_tag_value('don_flank',$don_flank);
		$newintron->add_tag_value('acceptor_seq',$acc_seq);
		$newintron->add_tag_value('donor_seq',$don_seq);
		scoreIntron($newintron);
		push @lastGRPgff, $newintron;
	    
	    $DONOR=$feature->end + 1;
	}
	$prevlen = $prevlen + $feature->end - $feature->start +1;
    }
    if (scalar @lastGRPgff) {

	if ($lastGRPgff[0]->strand eq "1") {
	    my $counter = 1;
	    foreach my $int (@lastGRPgff) {
		$int->remove_tag('ID') if $int->has_tag('ID');
		$int->add_tag_value('ID',$transcripts{$t}->display_name."_intron_$counter");
		$int->add_tag_value('inum',$counter);
		$gffout->write_feature($int);
		$counter++;
	    }
	} else {
	    my $counter = scalar @lastGRPgff;
	    foreach my $int (@lastGRPgff) {
		$int->remove_tag('ID') if $int->has_tag('ID');
		$int->add_tag_value('ID',$transcripts{$t}->display_name."_intron_$counter");
		$int->add_tag_value('inum',$counter);
		$gffout->write_feature($int);
		$counter--;
	    }
	}
	@lastGRPgff = ();
    }
}
$gffout->close;
$gfferror->close;

sub get_tseq{
    my $tx = shift;
    my $fadb = shift;
    my @exons = $tx->exons();
    my $strand = $exons[0]->strand;
    my $seq = "";
    if ($strand > 0) {
	foreach my $feature (sort {$a->start <=> $b->start} $tx->exons) {
	    my $s  = SeqOp::get_seq_BioDBFasta($fadb, $feature->seq_id,$feature->start, $feature->end, "+");
	    $seq .= $s;
	}
    } else {
	foreach my $feature (sort {$b->start <=> $a->start} $tx->exons) {
	    my $s  = SeqOp::get_seq_BioDBFasta($fadb, $feature->seq_id,$feature->start, $feature->end, "-");
	    $seq .= $s;
	}
    }
    return $seq;
}
sub get_iseq{
    my $tx = shift;
    my $fadb = shift;
    my @exons = $tx->exons();
    my $strand = $exons[0]->strand;
    my $seq = "";
    if ($strand > 0) {
	foreach my $feature (sort {$a->start <=> $b->start} $tx->exons) {
	    my $s  = SeqOp::get_seq_BioDBFasta($fadb, $feature->seq_id,$feature->start, $feature->end, "+");
	    $seq .= $s;
	}
    } else {
	foreach my $feature (sort {$b->start <=> $a->start} $tx->exons) {
	    my $s  = SeqOp::get_seq_BioDBFasta($fadb, $feature->seq_id,$feature->start, $feature->end, "-");
	    $seq .= $s;
	}
    }
    return $seq;
}

sub scoreIntron{
    my $intron = shift;
    my %DonGTAGScore;
    my %AccGTAGScore;
    my %DonATACScore;
    my %AccATACScore;
    my %DonU2Score;
    my %AccU2Score;
    my %DonSeq;
    my %AccSeq;
    my %names;
    my %AccName;
    my %AccGFF;
    my %FullName;
    my %scoreatt;
    my $seqname;
    my $adn = "";
    my $ddn = "";
    my $counter = 1;
    my $DonGTAGScore;
    my $AccGTAGScore;
    my $DonATACScore;
    my $AccATACScore;
    my $DonU2Score;
    my $AccU2Score;
    my $DonSeq;
    my $AccSeq;
    my $names;
    my $AccName;
    my $AccGFF;
    my $FullName;
    my $fasta ="";
    $seqname = $counter++;
    $adn = "";
    $ddn = "";
    $seqname = '';

    # GFF data structure
    my $seqid = 0;
    my $source = 1;
    my $feature  = 2;
    my $start    = 3;
    my $end      = 4;
    my $score    = 5;
    my $strand   = 6;
    my $frame = 7;
    my $group = 8;
    my %original_att;
    foreach my $key ( $intron->get_all_tags ) {
	my @values = $intron->get_tag_values($key);
	$original_att{$key}=$values[0];
    }
    my $faname = $intron->seq_id;
    $scoreatt{U2_acceptor_score} = -100;
    $scoreatt{U12_acceptor_score} = -100;
    $scoreatt{atype} = "unknown";
    $scoreatt{U2_donor_score} = -100;
    $scoreatt{U12_donor_score} = -100;
    $scoreatt{dtype} = "unknown";
    my $as_subtype = "unknown";
    my $ds_subtype = "unknown";
    my $u2ds_subtype = "unknown";
    my $u2as_subtype = "unknown";
    my $u12ds_subtype = "unknown";
    my $u12as_subtype = "unknown";
    my $u2bp_score = 0.00;
    my $u12bp_score = 0.00;
    my $ppt_score = 0.00;
    my $u2bp_pos = 0;
    my $u12bp_pos = 0;
    my $ppt_pos = 0;
    $DonSeq="";
    $AccSeq="";

    my $subseq  = substr($original_att{don_flank},-10,10).substr($original_att{donor_seq},0,15); #get_genome($gmapdb, $faname,$j, $k, $gff[6]);
    $scoreatt{donor} = substr($subseq, 10, 2);
    #my $frag = SeqOp::chr_subseq($fasta, $j,$k,"+",$int);
    my $dfrag_file = new File::Temp( DIR => $ENV{TMPDIR},SUFFIX => '.fa');
    my $afrag_file = new File::Temp( DIR => $ENV{TMPDIR},SUFFIX => '.fa');
    my $dgid = new File::Temp( DIR => $ENV{TMPDIR},SUFFIX => '.gff3');
    my $agid = new File::Temp( DIR => $ENV{TMPDIR},SUFFIX => '.gff3');
    open OUT, ">$dfrag_file";
    print OUT ">$faname\n$subseq\n";
    close OUT;
    `$exec -UW3do -P $param $dfrag_file > $dgid`;
    open (DTMP, "<$dgid");


    my @geneidGFF = ();
    while (<DTMP>) {
	next if  m/^#/;
	chomp;
	@geneidGFF = ();
	@geneidGFF = split;
	$geneidGFF[$group] =~ s/;$//;
	my %att = split(/[=;]/,$geneidGFF[$group]);
	if ( $geneidGFF[$start] eq "10") {
	    if ($att{type} eq "U2") {
		if ($geneidGFF[$score] > $scoreatt{U2_donor_score}) {
		    $scoreatt{U2_donor_score} = $geneidGFF[$score];
		    $u2ds_subtype = $att{subtype} if exists $att{subtype};
		    #$scoreatt{dseq} = $att{seq}  if exists $att{seq};
		}
	    } else {
		if ($geneidGFF[$score] > $scoreatt{U12_donor_score}) {
		    $scoreatt{U12_donor_score} = $geneidGFF[$score];
		    $u12ds_subtype = $att{subtype} if exists $att{subtype};
		    #$scoreatt{dseq} = $att{seq}  if exists $att{seq};
		}
	    }
	}
    }
    close DTMP;
    if (($scoreatt{U12_donor_score} > 0)&&($scoreatt{U12_donor_score} > $scoreatt{U2_donor_score})) {
	$scoreatt{dtype} = "U12";
	$ds_subtype = $u12ds_subtype;
    } elsif ($scoreatt{U2_donor_score} > -100) {
	$scoreatt{dtype} = "U2";
	$ds_subtype = $u2ds_subtype;
    } else {
	$scoreatt{dtype} = "unknown";
	$ds_subtype = "unknown";
    }

    #now acceptor on plus

    my $accsubseq  = substr($original_att{acceptor_seq},-53,53).substr($original_att{acc_flank},0,10);
    $scoreatt{acceptor} = substr($accsubseq, 51, 2);
    open OUT, ">$afrag_file";
    print OUT ">$faname\n$accsubseq\n";
    close OUT;
    `$exec -UW3ao -P $param $afrag_file > $agid`;
    open (TMP, "<$agid");

    @geneidGFF = ();
    while (<TMP>) {
	next if  m/^#/;
	chomp;
	@geneidGFF = ();
	@geneidGFF = split;
	$geneidGFF[$group] =~ s/;$//;
	my %att = split(/[=;]/,$geneidGFF[$group]);
	if ( $geneidGFF[$start] eq "53") {
	    if ($att{type} eq "U2") {
		if ($geneidGFF[$score] > $scoreatt{U2_acceptor_score}) {
		    $scoreatt{U2_acceptor_score} = $geneidGFF[$score];
		    $u2as_subtype = $att{subtype} if exists $att{subtype};
		    $scoreatt{U2_BP_score} = $att{bp_score} if exists $att{bp_score};
		    $scoreatt{ppt_score} = $att{ppt_score} if exists $att{ppt_score};
		    $scoreatt{U2_BP_pos} = $att{bp_pos} if exists $att{bp_pos};
		    $scoreatt{ppt_pos} = $att{ppt_pos} if exists $att{ppt_pos};
		    #$scoreatt{aseq} = $att{seq}  if exists $att{seq};
		}
	    } else {
		if ($geneidGFF[$score] > $scoreatt{U12_acceptor_score}) {
		    $scoreatt{U12_acceptor_score} = $geneidGFF[$score];
		    $u12as_subtype = $att{subtype} if exists $att{subtype};
		    $scoreatt{U12_BP_score} = $att{bp_score} if exists $att{bp_score};
		    $scoreatt{U12_BP_pos} = $att{bp_pos} if exists $att{bp_pos};
		    #$scoreatt{aseq} = $att{seq}  if exists $att{seq};
		}
	    }
	}
    }
    close TMP;
    if (($scoreatt{U12_acceptor_score} > 0)&&($scoreatt{U12_acceptor_score} > $scoreatt{U2_acceptor_score})) {
	$scoreatt{atype} = "U12";
	$as_subtype = $u12as_subtype;
    } elsif ($scoreatt{U2_acceptor_score} > -100) {
	$scoreatt{atype} = "U2";
	$as_subtype = $u2as_subtype;
    } else {
	$scoreatt{atype} = "unknown";
	$as_subtype = "unknown";
    }

    ### determine intron type
    if (($scoreatt{U2_donor_score} + $scoreatt{U2_acceptor_score})>($scoreatt{U12_donor_score} + $scoreatt{U12_acceptor_score})) {
	$scoreatt{introntype}="U2";
    } elsif (2*($scoreatt{U12_donor_score} + $scoreatt{U12_acceptor_score})>= 4 && ($scoreatt{U12_donor_score} > -1) && ($scoreatt{U12_acceptor_score} > -1)) {
	$scoreatt{introntype}="U12";
    } elsif ($scoreatt{dtype} ne "unknown" && $scoreatt{atype} ne "unknown") {
	$scoreatt{introntype}="U2";
    } else {
	$scoreatt{introntype}="non-canonical";
    }
    foreach my $nk (sort keys %scoreatt) {
	$intron->add_tag_value($nk,$scoreatt{$nk});
    }
    $intron->add_tag_value('loc',$intron->seq_id.":".$intron->start."-".$intron->end)
}
