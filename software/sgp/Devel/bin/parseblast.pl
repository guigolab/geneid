#!/usr/bin/perl 
#
# $Id: parseblast.pl,v 1.2 2000-07-18 16:51:38 jabril Exp $
#

use strict;
use Getopt::Long;

Getopt::Long::Configure("bundling");

my $PROGRAM = " parseblast.pl ";
my $VERSION = ' $Id: parseblast.pl,v 1.2 2000-07-18 16:51:38 jabril Exp $ ';

my ($hsp_flg, $gff_flg, $fullgff_flg, $aplot_flg, 
	$comment_flg, $split_flg, $help_flg, $err_flg,
    $default_flg, $aln_flg, $bit_flg, $ids_flg
	) = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
my ($prt, $main, $seqflg, $hsp, $fragment, $param
	) = (0, 0, 0, 0, 0, 0);
my ($prog_params, $program, $version, $seqname);
my ($scQ, $scS);   # reSCale_Query/Subject : 1 to re-scalate lengths.
my ($query_name, $db_name, $score, $descr,
	$ori, $end, $seq, $txt, $tt, $pt);
my (@seqlist, %prgseq, %dbase, %query, %cnt,
	%desc, %sco, %hsp_start, %hsp_end, %hsp_seq);

###################################################
## Getting command-line options.
###################################################

GetOptions( "H|hsp"            => \$hsp_flg     ,
			"G|gff"            => \$gff_flg     ,
			"F|fullgff"        => \$fullgff_flg ,
			"A|aplot"          => \$aplot_flg   ,
			"b|bit-score"      => \$bit_flg     ,
			"i|identity-score" => \$ids_flg     ,
			"c|comments"       => \$comment_flg ,
			"v|verbose"        => \$err_flg     ,
			"h|help|\?"        => \$help_flg     );
#			"a|align-score "   => \$aln_flg     ,
#			"s|split-output"   => \$split_flg   ,
#   -a, --align-score    : set <score> field to Alignment Score.
#	-s, --split-output : output each sequence match in a separate file in the current directory.

{   # first choose disables any other command-line option.
	$hsp_flg     && ($gff_flg = 0, $fullgff_flg = 0, $aplot_flg = 0, last);
	$fullgff_flg && ($hsp_flg = 0, $gff_flg     = 0, $aplot_flg = 0, last);
	$gff_flg     && ($hsp_flg = 0, $fullgff_flg = 0, $aplot_flg = 0, last);
	$aplot_flg   && ($hsp_flg = 0, $fullgff_flg = 0, $gff_flg   = 0, last);
	$default_flg = 1;
};

{   # first choose disables any other command-line option.  
    $bit_flg && ( $aln_flg = 0, $ids_flg = 0, last);
    $ids_flg && ( $aln_flg = 0, $bit_flg = 0, last);
	$aln_flg = 1;
};

($help_flg) && &prt_help;

###################################################
## Subroutines
###################################################

#
# Print help to STDERR.

sub prt_help {
	open(HELP, "| more") ;
	print HELP <<EndOfHelp;
PROGRAM:
          $PROGRAM
          $VERSION

USAGE:	parseblast.pl [options] <blast.results>

COMMAND-LINE OPTIONS:

	-H, --hsp            : prints output in HSP format.
	-G, --gff            : prints output in GFF format.
	-F, --fullgff        : prints output in GFF "alignment" format.
	-A, --aplot          : prints output in APLOT "GFF" format.
	-b, --bit-score      : set <score> field to Bit Score (default Alignment Score).
	-i, --identity-score : set <score> field to Identity Score (default Alignment).
	-c, --comments       : include parameters from blast program as comments.
	-v, --verbose        : warnings sent to Standard Error.
	-h, --help           : shows this help pages.

OUTPUT FORMATS:

	"S_" stands for "Subject_Sequence" and "Q_" for "Query_Sequence".
	<Program> name is taken from input blast file.
	<Strands> are calculated from <start> and <end> positions on original blast file.
	<Frame> is obtained from the blast file if is present else is set to ".".
    <SCORE> is set to Alignment Score by default, you can change it with "-b" and "-i".

[GFF]
 <Q_Name> <Program> hsp <Q_Start> <Q_End> <SCORE> <Q_Strand> <Q_Frame> <S_Name>

[FULL GFF] (Alignment)
 <Q_Name> <Program> hsp <Q_Start> <Q_End> <SCORE> <Q_Strand> <Q_Frame> ...
   ... Target "<S_Name>" <S_Start> <S_End> E_value <E_Value> Strand <S_Strand> Frame <S_Frame>

[APLOT] (Alignment)
 <Q_Name>:<S_Name> <Program> hsp <Q_Start>:<S_Start> <Q_End>:<S_End> <SCORE> ...
   ... <Q_Strand>:<S_Strand> <Q_Frame>:<S_Frame> <BitScore>:<HSP_Number> \# E_value <E_Value>

[HSP]
 <Program> <DataBase> : ...
   ... <IdentityMatches> <Min_Length> <IdentityScore> <AlignmentScore> <BitScore> <E_Value> <P_Sum> : ...
   ... <Q_Name> <Q_Start> <Q_End> <Q_Strand> <Q_Frame> : ...
   ... <S_Name> <S_Start> <S_End> <S_Strand> <S_Frame> : <S_FullDescription>

[DEFAULT]
 SCORE(<HSP_Number>): <AlignmentScore>
 BITSC(<HSP_Number>): <BitScore>
 EXPEC(<HSP_Number>): <E_Value> Psum(<P_Sum>)
 IDENT(<HSP_Number>): <IdentityMatches>/<Min_Length> : <IdentityScore> \%
 T_GAP(<HSP_Number>): <TotalGaps(BothSeqs)>
 FRAME(<HSP_Number>): <Q_Frame>/<S_Frame>
 MXLEN(<HSP_Number>): <Max_Length>
 QUERY(<HSP_Number>): length <Q_Length> : gaps <Q_TotalGaps> : ...
   ... <Q_Start> <Q_End> : <Q_Strand> : <Q_Frame> : <Q_FullSequence>
 SBJCT(<HSP_Number>): length <S_Length> : gaps <S_TotalGaps> : ...
   ... <S_Start> <S_End> : <S_Strand> : <S_Frame> : <S_FullSequence>

AUTHOR:  $PROGRAM (C) 2000 - Josep F. Abril

EndOfHelp
	close(HELP);
exit(1);
# print STDOUT "HSP:$hsp_flg GFF:$gff_flg COM:$comment_flg SPL:$split_flg HLP:$help_flg\n";
}

#
# Get new lines while empty line is not found and append to last line.

sub get_lines {
	my $spc = $_[0];
	my $tmp;
	# local($tmp);
	while (<>) {
		last if /^\s*$/;
		# print STDERR "$_";
		chop;
		s/^\s*/$spc/;
		$tmp .= $_;
	};
	$tmp;
}

#
# Getting scores from scoring vector extracted from HSP record.

sub get_scores {
	my $t = $_[0]; 
	my ($sc, $bt, $ex, $pv, $id);
	my ($frq, $frs) = ('.', '.');
	(($t =~ /Score[^\s]*\s+=\s+\b(\d+)\b\s+\([^,]*,/) || ($t =~ /Score\s+=\s+[^,]*\s+\((\d+)\)[^,]*,/)) &&
		($sc = $1);
	($t =~ /Score[^\s]*\s+=.*[\s\(]([+-]?(\d+\.?\d*|\.\d+))\b \bbits[^,]*,/) && 
		($bt = $1);
	($t =~ /Expect[^\s]*\s+=\s+([+-]?([Ee][+-]?\d+|(\d+\.?\d*|\.\d+)([Ee][+-]?\d+)?))\s*,/) && 
		($ex = $1);
	($ex =~ /^[Ee]/) && ($ex = "1".$ex);
	($t =~ /Sum[^\s]*\s+\bP[^\s]*\s+=\s+([+-]?([Ee][+-]?\d+|(\d+\.?\d*|\.\d+)([Ee][+-]?\d+)?))\s*,/) ?
		($pv = $1) : ($pv = $ex);
	($t =~ /Identities[^\s]*\s+=\s+(\d+)\/(\d+)\s+/) && ($id = $1);
	( $scQ && !$scS) && do {   # BLASTX (translated nucleotides vs protein)
		($t =~ /Frame[^\s]*\s+=\s+(\+|\-)(\d)\b/) && ($frq = $2);
	};
	(!$scQ &&  $scS) && do {   # TBLASTN (protein vs translated nucleotides)
		($t =~ /Frame[^\s]*\s+=\s+(\+|\-)(\d)\b/) && ($frs = $2);
	};
	( $scQ &&  $scS) && do {   # TBLASTX (translated nucleotides vs translated nucleotides)
		($t =~ /Frame[^\s]*\s+=\s+(\+|\-)(\d)\s*\/\s*(\+|\-)(\d)\b/) && ($frq = $2, $frs = $4);
	};
	$sc, $bt, $ex, $pv, $id, $frq, $frs;
}

sub chk_strand {
	my ($first, $last) = @_ ;
    my $kk = $first;
    my $st = "+";
	($first >= $last) && ($first = $last, $last = $kk, $st = "-");
	$first, $last, $st;
}

#
# Formatting output as plain, HSPs or GFF.

sub prt_out {
	my ($n, $nm, $tq, $ts,             # couNter, NaMe, TagQuery, TagSubject
		$sc, $bt, $ex, $pv, $id,       # SCore, BiTscore, EvalueX, IDentityscore
		$qfr, $sfr,                    # QueryFRame, SubjectFRame
		$qst, $sst,                    # QuerySTrand, SubjectSTrand
		$gsc, $prg,                    # GroupSCore, PRoGram
		$hsq, $hss, $heq, $hes,        # HspStartQuery, HspStartSubject, HspEndQuery, HspEndSubject
		$lnq, $lns, $lnmin, $lnmax,    # LeNgthQuery, LeNgthSubject, LeNgthMINqueyxsubject, LeNgthMAXqueyxsubject
		$lnmx, $gpq, $gps, $gpt        # LeNgthMaXhspseq, GaPQuery, GaPSubject, GaPTotal
		);
	if $comment_flg  { 
		print STDOUT $prog_params."\n"."#\n"; 
	} else {
		print STDOUT "#\n# $program $version\n";
	};
	while (@seqlist) {
		$nm = shift(@seqlist);
		(!$hsp_flg) && (print STDOUT "#\n# $prgseq{$nm} :: DB $dbase{$nm} :: MATCH $query{$nm}x$nm :: $cnt{$nm} HSPs\n# DESCR: $desc{$nm}\n#\n");
		($cnt{$nm}>0) && do {
			for ($n = 1; $n <= $cnt{$nm}; $n++) {
				$tq = $nm."query".$n;
				$ts = $nm."sbjct".$n;
				($sc, $bt, $ex, $pv, $id, $qfr, $sfr) = &get_scores($sco{$nm.$n});
				($hsq, $heq, $qst) = &chk_strand($hsp_start{$tq}, $hsp_end{$tq});
				($hss, $hes, $sst) = &chk_strand($hsp_start{$ts}, $hsp_end{$ts});
				$lnq = $heq - $hsq + 1 ;
				$scQ && ($lnq = $lnq / 3) ;
				$lns = $hes - $hss + 1 ; 
				$scS && ($lns = $lns / 3) ;
				$lnmin = ($lnq>$lns) ? $lns : $lnq;
				$lnmax = ($lnq<$lns) ? $lns : $lnq;
				$lnmx = length($hsp_seq{$tq});
				{ my $hh = $hsp_seq{$tq}; $gpq = ($hh =~ s/-/ /g) || 0; };
				{ my $hh = $hsp_seq{$ts}; $gps = ($hh =~ s/-/ /g) || 0; };
				$gpt = $gpq + $gps;
				{
					($ids_flg || $default_flg || $hsp_flg) && # score is Identities divided by minlength
						(($gsc) = eval(($id/$lnmin)*100) =~ /^(\d+(\.\d{0,3})?)/, last);
					$bit_flg && ($gsc = $bt, last);
					$gsc = $sc;
				};
				($prg) = $prgseq{$nm} =~ /^([^\s]+)\s/;
				$gff_flg && do { 
					print STDOUT <<"EndOfGFF";
$query{$nm}\t$prg\thsp\t$hsq\t$heq\t$gsc\t$qst\t$qfr\t$nm\t\# E_value $ex : P_sum $pv
EndOfGFF
}; # gff_flg
				$fullgff_flg && do { 
					print STDOUT <<"EndOfFullGFF";
$query{$nm}\t$prg\thsp\t$hsq\t$heq\t$gsc\t$qst\t$qfr\tTarget \"$nm\"\t$hss\t$hes\tE_value $ex\tStrand $sst\tFrame $sfr
EndOfFullGFF
}; # fullgff_flg
				$aplot_flg && do { 
					print STDOUT <<"EndOfAPLOT";
$query{$nm}:$nm\t$prg\thsp\t$hsq:$hss\t$heq:$hes\t$gsc\t$qst:$sst\t$qfr:$sfr\t$bt:$n\t\# E_value $ex : P_sum $pv
EndOfAPLOT
}; # gff_flg
				$hsp_flg && do {
					print STDOUT <<"EndOfHSPs";
$prg $dbase{$nm} : $id $lnmin $gsc $sc $bt $ex $pv : $query{$nm} $hsq $heq $qst $qfr : $nm $hss $hes $sst $sfr : $desc{$nm}
EndOfHSPs
};
				$default_flg && do {
					print STDOUT <<"EndOfPlain";
SCORE($n): $sc\nBITSC($n): $bt\nEXPEC($n): $ex Psum($pv)
IDENT($n): $id/$lnmin : $gsc \%
T_GAP($n): $gpt\nFRAME($n): $qfr/$sfr\nMXLEN($n): $lnmx
QUERY($n): length $lnq : gaps $gpq : $hsq $heq : $qst : $qfr : $hsp_seq{$tq}
SBJCT($n): length $lns : gaps $gps : $hss $hes : $sst : $sfr : $hsp_seq{$ts}
EndOfPlain
};
			} # for $cnt{$nm}
		} # do if $cnt{$nm}>0
	} # foreach
#	undef %hsp_start; undef %hsp_end; undef %hsp_seq; undef %sco; undef %cnt; # undef %seqlist;
}

###################################################
## Main Loop
###################################################

while (<>) {
	$pt++;
	$err_flg && do {
	  FOO: { my $kk = $pt % 50;
			 ($kk == 0) && ((print STDERR ".[$pt]\n"), last FOO);
			 print STDERR ".";
		 }; 
	};
	# s/\r\n$/\n/; # if your input records finish with "\r\n" (like EMBL).
	next if /^\s*$/;  # /^\s*$/ is equal to AWK /^[ \t]*$/
	my $tmpinput = $_;
	chop;
    # print STDOUT "$. : $_ \n"; # "$." is record number && "$_" is whole record
  CHECK: {
	  /^\s*T?BLAST[PNX]?/        && do { # Starts with "T?BLAST[PNX]?" ?
		  # print STDERR "$_\n";
		  $prt && &prt_out;
		  ($program, $version) = split;
		  # typeQ/typeS: 0 for proteins - 1 for nucleic acids.
		  ($program =~ /^BLASTP$/ ) && ( $scQ = 0, $scS = 0); # Amino Acids vs Amino Acids
		  ($program =~ /^BLASTN$/ ) && ( $scQ = 0, $scS = 0); # Nucleotides vs Nucleotides
		  ($program =~ /^BLASTX$/ ) && ( $scQ = 1, $scS = 0); # Nucleotides vs Amino Acids 
		  ($program =~ /^TBLASTN$/) && ( $scQ = 0, $scS = 1); # Amino Acids vs Nucleotides
		  ($program =~ /^TBLASTX$/) && ( $scQ = 1, $scS = 1); # Nucleotides translated vs Nucleotides translated
		  $prog_params = "#\n# $program $version\n#";
		  $query_name = $db_name = '';
		  $main = 1;
		  $seqflg = $hsp = $fragment = $param = 0;
		  last CHECK; 
	  }; # /^\s*T?BLAST[PNX]?/
	  /^>/                    && do { # Starts with ">" ?: sequences.
		  # print STDERR "$_\n";
		  ($seqname, $descr) = split(/\s+/, $_, 2);
		  $seqname =~ s/^\s*>//;
		  $seqname =~ s/:|\|/_/g;
		  $prgseq{$seqname} = "$program ($version)";
		  $query{$seqname} = $query_name;
		  ($db_name =~ /([^\/]+)$/) && ($dbase{$seqname} = $1);
		  push(@seqlist,$seqname);
		  $desc{$seqname} = join(' ', $descr, &get_lines(' '));
#		  $cnt{$seqname} = 0;
		  $seqflg = 1;
		  $main = $hsp = $fragment = 0;
		  last CHECK;
	  }; # /^\s*>/
	  (/^\s*Score/ && $seqflg) && do { # Starts with "Score" ?: HSPs.
		  # print STDERR "$_\n";
		  ( $score = $_ ) =~ s/^\s*//;
		  $cnt{$seqname}++;
		  $sco{$seqname.$cnt{$seqname}} = join('', $score, &get_lines(', '));
		  # print STDOUT $sco{$seqname.$cnt{$seqname}}."\n";
		  $hsp = 1;
		  $fragment = 0;
		  last CHECK;
	  }; # (/^\s*Score/ && ($seqflg))
	  (/^Query:/ && $hsp)    && do { # Starts with "Query" ?: Fragments.
		  # print STDERR "$_\n";
		  $fragment = 1;          
		  last CHECK;
	  }; # (/^\s*Query/ && ($hsp))
	  (/^\s*(?:Database|Parameters):/ && !$main) && do { # Parameters Section.
		  # print STDERR "$_\n";
		  $prt = $param = 1;
		  $seqflg = $hsp = $fragment = 0;
		  last CHECK;
	  }; # (/^\s*(?:Database|Parameters)/ && ($hsp))
  } # CHECK Block
	# print STDOUT "$. : MAIN=$main SEQFLG=$seqflg HSP=$hsp FRAGMENT=$fragment => SEQNAME=$seqname\n";
  LOAD: {
	  $fragment && do { # We are within a fragment.
		  $txt = '';
		  $tt = $cnt{$seqname};
		  ($txt,$ori,$seq,$end) = split;
		  if ($txt =~ /^Query:/) {
			  # print STDERR "$_\n";
			  ($hsp_start{$seqname."query".$tt}) || ($hsp_start{$seqname."query".$tt} = $ori);
			  $hsp_end{$seqname."query".$tt} = $end;
			  $hsp_seq{$seqname."query".$tt} .= $seq;
		  } # if ($txt =~ /Query/)
		  elsif ($txt =~ /^Sbjct:/) {
			  # print STDERR "$_\n";
			  ($hsp_start{$seqname."sbjct".$tt}) || ($hsp_start{$seqname."sbjct".$tt} = $ori);
			  $hsp_end{$seqname."sbjct".$tt} = $end;
			  $hsp_seq{$seqname."sbjct".$tt} .= $seq;
			  $fragment = 0;
		  } # elsif ($txt =~ /Sbjct/)
		  else { last LOAD; };
	  }; # ($fragment)
	  $main && do { # We are within the blast file header.
		  /^\s*Query= +(.*)\s*$/ && do {
			  # print STDERR "$_\n";
			  ($query_name = $1) =~ s/:|\|/_/g ;
		  };
		  /^\s*Database: +(.*)\s*$/ && do { 
			  # print STDERR "$_\n";
			  $db_name = $1;
			  while (<>) {
				  last if /^(?:.*\bletter.*|\s*)$/;
				  # print STDERR "$_\n";
				  chop;
				  s/^\s*//;
				  s/\s*$//;
				  $db_name .= $_;
 			  }; # while getline
		  }; # /^\s*Database: +(.*)\s*$/
		  last LOAD;
	  }; # ($main)
	  $param && do { # We are within the blast file trailer.
		  if (/^\s*[^\[\<\-]/) {
			  # print STDERR "$_\n";
			  chop;
			  s/^/\n\# /;
			  $prog_params = join('', $prog_params, $_);
		  } else { $param = 0; }
		  last LOAD;
	  }; # ($param)
  } # LOAD Block
	close(ARGV) if (eof);
} # while

print STDERR ".[$pt]\n" if $err_flg;

$prt && &prt_out;

exit(0);
