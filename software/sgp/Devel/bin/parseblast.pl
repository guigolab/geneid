#!/usr/bin/perl 
#
# $Id: parseblast.pl,v 1.4 2000-07-27 02:02:46 jabril Exp $
#

use strict;
use Getopt::Long;

Getopt::Long::Configure("bundling");

my $PROGRAM = "parseblast.pl";
my $VERSION = '$Id: parseblast.pl,v 1.4 2000-07-27 02:02:46 jabril Exp $ ';

my ($hsp_flg, $gff_flg, $fullgff_flg, $aplot_flg, $subject_flg,
	$comment_flg, $nocmmnt_flg, $split_flg, $help_flg, $err_flg,
    $expanded_flg, $pairwise_flg, $aln_flg, $bit_flg, $ids_flg);
my ($prt, $main, $seqflg, $hsp, $fragment, $param) = (0, 0, 0, 0, 0, 0);
my ($prog_params, $program, $version, $seqname);
my ($scQ, $scS);   # reSCale_Query/Subject : 1 to re-scalate lengths.
my ($query_name, $db_name, $score, $descr,
	$ori, $end, $seq, $txt, $tt, $pt);
my (@seqlist, %prgseq, %dbase, %query, %cnt,
	%desc, %sco, %hsp_start, %hsp_end, %hsp_seq);

###################################################
## Getting command-line options.
###################################################

GetOptions( "G|gff"            => \$gff_flg      ,
			"F|fullgff"        => \$fullgff_flg  ,
			"A|aplot"          => \$aplot_flg    ,
            "S|subject"        => \$subject_flg  ,
            "X|extended"       => \$expanded_flg ,
			"P|pairwise"       => \$pairwise_flg ,
			"b|bit-score"      => \$bit_flg      ,
			"i|identity-score" => \$ids_flg      ,
			"c|comments"       => \$comment_flg  ,
			"n|no-comments"    => \$nocmmnt_flg  ,
			"v|verbose"        => \$err_flg      ,
			"h|help|\?"        => \$help_flg      );
#			"a|align-score "   => \$aln_flg     ,
#			"s|split-output"   => \$split_flg   ,
#   -a, --align-score    : set <score> field to Alignment Score.
#	-s, --split-output : output each sequence match in a separate file in the current directory.

{   # first choose disables any other command-line option.
    $pairwise_flg
		&& (($expanded_flg,$gff_flg,$fullgff_flg,$aplot_flg,$hsp_flg)=(0,0,0,0,0), last);
	$gff_flg
		&& (($expanded_flg,$pairwise_flg,$fullgff_flg,$aplot_flg,$hsp_flg)=(0,0,0,0,0), last);
	$fullgff_flg
		&& (($expanded_flg,$pairwise_flg,$gff_flg,$aplot_flg,$hsp_flg)=(0,0,0,0,0), last);
	$aplot_flg
		&& (($expanded_flg,$pairwise_flg,$gff_flg,$fullgff_flg,$hsp_flg)=(0,0,0,0,0), last);
	$expanded_flg
		&& (($pairwise_flg,$gff_flg,$fullgff_flg,$aplot_flg,$hsp_flg)=(0,0,0,0,0), last);
	$hsp_flg = 1; $subject_flg = 0;
};

{   # first choose disables any other command-line option.  
    $bit_flg && (($aln_flg,$ids_flg) = (0,0), last);
    $ids_flg && (($aln_flg,$bit_flg) = (0,0), last);
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

USAGE:	parseblast.pl [options] <results.from.blast>

COMMAND-LINE OPTIONS:

  "$PROGRAM" prints output in "HSP" format by default (see below).
  It takes input from <STDIN> or single/multiple files, and writes
  its output to <STDOUT>, so user can redirect to a file but
  he also could use the program as a filter within a pipe. 
  Available options are:

    -G, --gff            : prints output in GFF format.
    -F, --fullgff        : prints output in GFF "alignment" format.
    -A, --aplot          : prints output in APLOT "GFF" format.
    -S, --subject        : projecting by SUBJECT in GFF output (default by QUERY).
    -X, --expanded       : expanded output (producing multiline output records).
    -b, --bit-score      : set <score> field to Bit Score (default Alignment Score).
    -i, --identity-score : set <score> field to Identity Score (default Alignment).
    -c, --comments       : include parameters from blast program as comments.
    -n, --no-comments    : do not print "#" lines (raw output without comments).
    -v, --verbose        : warnings sent to <STDERR>.
    -h, --help           : shows this help pages.

OUTPUT FORMATS:

  "S_" stands for "Subject_Sequence" and "Q_" for "Query_Sequence".
  <Program> name is taken from input blast file.
  <Strands> are calculated from <start> and <end> positions on original blast file.
  <Frame> is obtained from the blast file if is present else is set to ".".
  <SCORE> is set to Alignment Score by default, you can change it with "-b" and "-i".
  If "-S" or "--subject" options are given, then QUERY fields are referred to SUBJECT
  and SUBJECT fields are relative to QUERY (this only available for GFF output records).
  Dots ("...") mean that record description continues in the following line,
  but such record is printed as a single line by "$PROGRAM".

[HSP]  <- (This is the DEFAULT OUTPUT FORMAT)
 <Program> <DataBase> : ...
   ... <IdentityMatches> <Min_Length> <IdentityScore> <AlignmentScore> <BitScore> <E_Value> <P_Sum> : ...
   ... <Q_Name> <Q_Start> <Q_End> <Q_Strand> <Q_Frame> : ...
   ... <S_Name> <S_Start> <S_End> <S_Strand> <S_Frame> : <S_FullDescription>

[GFF]
 <Q_Name> <Program> hsp <Q_Start> <Q_End> <SCORE> <Q_Strand> <Q_Frame> <S_Name>

[FULL GFF]  <- (GFF showing alignment data)
 <Q_Name> <Program> hsp <Q_Start> <Q_End> <SCORE> <Q_Strand> <Q_Frame> ...
   ... Target "<S_Name>" <S_Start> <S_End> E_value <E_Value> Strand <S_Strand> Frame <S_Frame>

[APLOT]  <- (GFF format enhanced for APLOT program)
 <Q_Name>:<S_Name> <Program> hsp <Q_Start>:<S_Start> <Q_End>:<S_End> <SCORE> ...
   ... <Q_Strand>:<S_Strand> <Q_Frame>:<S_Frame> <BitScore>:<HSP_Number> \# E_value <E_Value>

[EXPANDED]
 MATCH(<HSP_Number>): <Q_Name> x <S_Name>
 SCORE(<HSP_Number>): <AlignmentScore>
 BITSC(<HSP_Number>): <BitScore>
 EXPEC(<HSP_Number>): <E_Value> Psum(<P_Sum>)
 IDENT(<HSP_Number>): <IdentityMatches>/<Min_Length> : <IdentityScore> \%
 T_GAP(<HSP_Number>): <TotalGaps(BothSeqs)>
 FRAME(<HSP_Number>): <Q_Frame>/<S_Frame>
 STRND(<HSP_Number>): <Q_Strand>/<S_Strand>
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
	my ($qfr, $sfr) = ('.', '.');
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
		($t =~ /Frame[^\s]*\s+=\s+(\+|\-)(\d)\b/) && ($qfr = $2);
	};
	(!$scQ &&  $scS) && do {   # TBLASTN (protein vs translated nucleotides)
		($t =~ /Frame[^\s]*\s+=\s+(\+|\-)(\d)\b/) && ($sfr = $2);
	};
	( $scQ &&  $scS) && do {   # TBLASTX (translated nucleotides vs translated nucleotides)
		($t =~ /Frame[^\s]*\s+=\s+(\+|\-)(\d)\s*\/\s*(\+|\-)(\d)\b/) && ($qfr = $2, $sfr = $4);
	};
	$sc, $bt, $ex, $pv, $id, $qfr, $sfr;
}

sub chk_strand {
	my ($first, $last) = @_ ;
    my $st = "+";
	($first >= $last) && (($first, $last) = ($last, $first), $st = "-");
	$first, $last, $st;
}

#
# Formatting output as plain, HSPs or GFF.

sub prt_out {
	my ($n, $nm, $tq, $ts,             # couNter, NaMe, TagQuery, TagSubject
		$sc, $bt, $ex, $pv, $id,       # SCore, BiTscore, EvalueX, IDentityscore
		$frq, $frs,                    # QueryFRame, SubjectFRame
		$stq, $sts,                    # QuerySTrand, SubjectSTrand
		$gsc, $prg,                    # GroupSCore, PRoGram
		$hsq, $hss, $heq, $hes,        # HspStartQuery, HspStartSubject, HspEndQuery, HspEndSubject
		$lnq, $lns, $lnmin, $lnmax,    # LeNgthQuery, LeNgthSubject, LeNgthMINqueyxsubject, LeNgthMAXqueyxsubject
		$lnmx, $gpq, $gps, $gpt        # LeNgthMaXhspseq, GaPQuery, GaPSubject, GaPTotal
		);
	if ($comment_flg)  { 
		print STDOUT $prog_params."\n"."#\n" if !$nocmmnt_flg; 
	} else {
		print STDOUT "#\n# $program $version\n" if !$nocmmnt_flg;
	};
	while (@seqlist) {
		$nm = shift(@seqlist);
		(!$hsp_flg && !$nocmmnt_flg) && (print STDOUT "#\n# $prgseq{$nm} :: DB $dbase{$nm} :: $cnt{$nm} HSPs for $query{$nm}x$nm \n# DESCR: $desc{$nm}\n#\n");
		($cnt{$nm}>0) && do {
			for ($n = 1; $n <= $cnt{$nm}; $n++) {
				$tq = $nm."query".$n;
				$ts = $nm."sbjct".$n;
				($sc, $bt, $ex, $pv, $id, $frq, $frs) = &get_scores($sco{$nm.$n});
				($hsq, $heq, $stq) = &chk_strand($hsp_start{$tq}, $hsp_end{$tq});
				($hss, $hes, $sts) = &chk_strand($hsp_start{$ts}, $hsp_end{$ts});
				$lnq = $heq - $hsq + 1 ;
				$scQ && ($lnq = $lnq / 3) ;
				$lns = $hes - $hss + 1 ; 
				$scS && ($lns = $lns / 3) ;
				$lnmin = ($lnq>$lns) ? $lns : $lnq;
				$lnmax = ($lnq<$lns) ? $lns : $lnq;
				$lnmx = length($hsp_seq{$tq});
				{
					my $hh = $hsp_seq{$tq};
					$gpq = ($hh =~ s/-/ /g) || 0; 
					$hh = $hsp_seq{$ts};
					$gps = ($hh =~ s/-/ /g) || 0;
				};
				$gpt = $gpq + $gps;
				{
					($ids_flg || $expanded_flg || $hsp_flg) && # score is Identities divided by minlength
						(($gsc) = eval(($id/$lnmin)*100) =~ /^(\d+(\.\d{0,3})?)/, last);
					$bit_flg && ($gsc = $bt, last);
					$gsc = $sc;
				};
				($prg) = $prgseq{$nm} =~ /^([^\s]+)\s/;
				# GFF format
			  PRINT:
				{
					$subject_flg && do {
						$gff_flg && do { 
							print STDOUT <<"EndOfGFF";
$nm\t$prg\thsp\t$hss\t$hes\t$gsc\t$sts\t$frs\t$query{$nm}\t\# E_value $ex : P_sum $pv
EndOfGFF
last PRINT;
						}; # gff_flg
						$fullgff_flg && do { 
							print STDOUT <<"EndOfFullGFF";
$nm\t$prg\thsp\t$hss\t$hes\t$gsc\t$sts\t$frs\tTarget \"$query{$nm}\"\t$hsq\t$heq\tE_value $ex\tStrand $stq\tFrame $frq
EndOfFullGFF
last PRINT;
						}; # fullgff_flg
						$aplot_flg && do { 
							print STDOUT <<"EndOfAPLOT";
$nm:$query{$nm}\t$prg\thsp\t$hss:$hsq\t$hes:$heq\t$gsc\t$sts:$stq\t$frs:$frq\t$bt:$n\t\# E_value $ex : P_sum $pv
EndOfAPLOT
last PRINT;
						}; # gff_flg
						last SBJCT;
					}; # $subject_flg
					$gff_flg && do { 
						print STDOUT <<"EndOfGFF";
$query{$nm}\t$prg\thsp\t$hsq\t$heq\t$gsc\t$stq\t$frq\t$nm\t\# E_value $ex : P_sum $pv
EndOfGFF
last PRINT;
					}; # gff_flg
					$fullgff_flg && do { 
						print STDOUT <<"EndOfFullGFF";
$query{$nm}\t$prg\thsp\t$hsq\t$heq\t$gsc\t$stq\t$frq\tTarget \"$nm\"\t$hss\t$hes\tE_value $ex\tStrand $sts\tFrame $frs
EndOfFullGFF
last PRINT;
					}; # fullgff_flg
					$aplot_flg && do { 
						print STDOUT <<"EndOfAPLOT";
$query{$nm}:$nm\t$prg\thsp\t$hsq:$hss\t$heq:$hes\t$gsc\t$stq:$sts\t$frq:$frs\t$bt:$n\t\# E_value $ex : P_sum $pv
EndOfAPLOT
last PRINT;
					}; # gff_flg					
				# GFF format
					$expanded_flg && do {
						print STDOUT <<"EndOfPlain";
MATCH($n): $query{$nm} x $nm
SCORE($n): $sc\nBITSC($n): $bt\nEXPEC($n): $ex Psum($pv)
IDENT($n): $id/$lnmin : $gsc \%
T_GAP($n): $gpt\nFRAME($n): $frq/$frs\nSTRND($n): $stq/$sts\nMXLEN($n): $lnmx
QUERY($n): length $lnq : gaps $gpq : $hsq $heq : $stq : $frq : $hsp_seq{$tq}
SBJCT($n): length $lns : gaps $gps : $hss $hes : $sts : $frs : $hsp_seq{$ts}
EndOfPlain
last PRINT;
					};
					$pairwise_flg && do {
						print STDOUT <<"EndOfALIGN";
$tq $hsp_seq{$tq}
$ts $hsp_seq{$ts}
EndOfALIGN
last PRINT;
					};
					$hsp_flg && do {
						print STDOUT <<"EndOfHSPs";
$prg $dbase{$nm} : $id $lnmin $gsc $sc $bt $ex $pv : $query{$nm} $hsq $heq $stq $frq : $nm $hss $hes $sts $frs : $desc{$nm}
EndOfHSPs
};
				} # PRINT
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

$err_flg && print STDERR ".[$pt]\n";

$prt && &prt_out;

exit(0);
