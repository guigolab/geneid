#!/usr/bin/perl
#
# GetSRsAln.pl - Obtaining Similarity Regions and its sequence from HSPs.
#
# $Id: GetSRsAln.pl,v 1.10 2000-08-11 21:21:29 jabril Exp $
#

my $PROGRAM = "GetSRsAln";
my $VERSION = "Version 1.0";
my $Start = time;

use strict;
use IPC::Open2;
use Getopt::Long;


##############################################################################
##                       Getting Comand-line Options                        ##
##############################################################################

Getopt::Long::Configure("bundling","pass_through");

my ( $verbose_flg, $bit_flg, $ids_flg,  # GETOPTs FLAGS
	 $aln_flg, $aplot_flg, $hsp_flg, $hsp_too_flg,
	 $only_Q_flg, $only_S_flg, $to_file, $to_file_flg,
	 $help_flg, $debug_flg, $cutoff_flg, $stdin_flg
	 ) = (0,0,0,0,0,0,0,0,0,0,0,0,0,"-",1);

GetOptions( "Q|query-only"      => \$only_Q_flg  ,
			"S|subject-only"    => \$only_S_flg  ,
			"A|aplot"           => \$aplot_flg   ,
			"P|pairwise"        => \$aln_flg     , # show only alignment not GFF
			"H|print-hsps:s"    => \$hsp_flg     ,
			"b|bit-score"       => \$bit_flg     ,
			"i|identity-score"  => \$ids_flg     ,
			"C|cutoff=f"        => \$cutoff_flg  ,
			"W|write-to-file:s" => \$to_file     ,
			"v|verbose"         => \$verbose_flg ,
			"D|debugging"       => \$debug_flg   ,
 			"h|help|\?"         => \$help_flg    ,
			);

&prt_help if $help_flg; 

my $ARGV = my $scoreopt = my $vrbopt = ""; 
$verbose_flg = 1 if $debug_flg;
$vrbopt = " -v" if $verbose_flg;
$scoreopt = " -i" if $ids_flg; # 'identity' score from BLAST
$scoreopt = " -b" if $bit_flg; # 'bit' score from BLAST
$only_S_flg = 0 if $only_Q_flg;
$aln_flg = 0 if $aplot_flg;
$hsp_flg = 1 if $hsp_flg eq ""; # Print also HSPs from input
$hsp_too_flg = 1 if $hsp_flg;
$to_file = 1 if $to_file eq ""; # Write output to file
$to_file_flg = 1 if $to_file;


##############################################################################
##                          Global Variables                                ##
##############################################################################

my $parseblast = "parseblast"; # $blastgff = "blast2gff -g"
my $pb_PROG   = "$parseblast$vrbopt$scoreopt -nFQ";
   # 'n'o comments # 'F'ullgff better then 'G'ff # Append se'Q'uence to GFF
   # only SUBJECT "-nGbS" # only ALN "-nPW"
my ( $blast_file, $blast_alnQ, $blast_alnS, $hsp_file, @data, $FLH);
my ( %hsp, %sr, @project_query, @project_subject,   # HSPs RECORDS
     @stack, %opened, $coord, $score, $index, $codeQ, $codeS,
	 );
my ( $sr_count, $ori, $end, $sco, $idx, $rec, $dif, # SRs RECORDS
	 $which_proj, $prj_flg, $Q_ori, $Q_end,
	 $S_ori, $S_end, $Q_seq, $S_seq, $srs, %sr_ctr,
	 );
my ( $current_coord, $last_coord, $current_score,   # Building SRs
	 $last_score, $current_index, $last_index, $rcd, $w_sr, $lines,
	 $op, $n, $c, $stv, $t_score, $t_index, %tfxs, $no_stdout,
	 );
my %STR = ( "+" => 0,
			"-" => 3 );
my %FRxST = ( 1 => "+1",
			  2 => "+2",
			  3 => "+3",
			  4 => "-1",
			  5 => "-2",
			  6 => "-3",
			  );
my $max_frame_num = keys %FRxST ;


##############################################################################
##                          Global Subroutines                              ##
##############################################################################

# Print help to STDERR.
#
sub prt_help {
        open(HELP, "| more") ;
        print HELP <<EndOfHelp;
PROGRAM:
        $PROGRAM
        $VERSION

USAGE:  $PROGRAM [options] <results.from.blast>

REQUIRES:

	This program needs to find "parseblast" in the path.
  (Don\'t you have your ~/bin dir in the path?)............ ;^) 
  
COMMAND-LINE OPTIONS:

    "$PROGRAM" prints all the Similarity Regions (SRs) that finds 
  from the HSPs given by the blast input (now works just with TBLASTX).
  It takes input from <STDIN> or single/multiple files, and writes
  its output to <STDOUT>, so user can redirect to a file but
  he also could use the program as a filter within a pipe. It also 
  can write to files with the "-W" option. The output then is splited
  onto two files with the ".alnQ" and ".alnS" extensions. By default
  SRs for both Query and Subject are generated, you can change this
  with "-Q" and "-S" options, see below. The basic output is in GFF, 
  you can switch to pairwise alignment output not in GFF.

    -Q, --query-only     : just print QUERY SRs (default both).
    -S, --subject-only   : just print SUBJECT SRs (default both).
    -A, --aplot          : prints output in APLOT "GFF" format.
    -P, --pairwise       : print pairwise alignment for each SR.
    -H, --print-hsps     : also includes HSPs in output (in GFF format).
                           As "-W" option, see below, a non-mandatory
                           parameter can be especified to send this
                           output to a file named "<file_name>.hsp".
    -b, --bit-score      : set <score> field to Bits (default Alignment Score).
    -i, --identity-score : set <score> field to Identities (default Alignment).
    -C <f>, --cutoff <f> : <f> is a real number for lower HSP score cutoff.
    -W, --write-to-file  : write output to separate files
                                 + for QUERY:   "<file_name>.alnQ"
                                 + for SUBJECT: "<file_name>.alnS"
                           you can provide "<file_name>" as parameter for
                           this option, as example, if you provide "-W results"
                           you send output to "results.alnQ" and "results.alnS".
                           If "<file_name>" is not provided, input filename is
                           set by default as "<file_name>".
    -v, --verbose        : warnings sent to <STDERR>.
    -D, --debugging      : extended report for debugging sent to <STDERR>.
    -h, --help           : show this help pages.

BUGS:    Report any problem to: abril\@imim.es

AUTHOR:  $PROGRAM is under GNU-GPL (C) 2000 - Josep F. Abril

EndOfHelp
        close(HELP);
exit(1);
} # END_SUB: prt_help
 
# Reporting IN/OUT progress.
#
sub prt_progress {
    $verbose_flg && do {
		print STDERR ".";
		(($_[0] % 50) == 0) && print STDERR "[".&fill_left($_[0],6,"0")."]\n";
	};
} # END_SUB: prt_progress
#
sub prt_foeprg {
    $verbose_flg && ((($_[0] % 50) != 0) && print STDERR "[".&fill_left($_[0],6,"0")."]\n" );
} # END_SUB: prt_foeprg

# Get a fixed length string from a given string and filling char/s.
#
sub fill_right { $_[0].($_[2] x ($_[1] - length($_[0]))) }
#
sub fill_left  { ($_[2] x ($_[1] - length($_[0]))).$_[0] }

# returns the max value from input array
#
sub max { my ($z) = shift @_; my $l; foreach $l (@_) { $z = $l if $l > $z ; }; $z; }

# Open a file or STDOUT
#
sub open_file_handle {
	my ($flg, $fl) = @_;
	$flg && do {
		open(HF,"> $fl");
		$FLH=*HF;
		return 1;
	};
	$FLH=*STDOUT;
	return 0;
} # END_SUB: open_file_handle
 
# Timing.
#
sub get_exec_time {
	$verbose_flg && do {
		my $End = @_;
        my ($c,$s,$m,$h,$r);
        $r = $End - $Start;
        $s = $r % 60;
        $r = ($r - $s) / 60;
        $m = $r % 60;
        $r = ($r - $m) / 60;
        $h = $r % 24;
		($s,$m,$h) = (&fill_left($s,2,"0"),&fill_left($m,2,"0"),&fill_left($h,2,"0"));
print STDERR <<EOF;
##
##########################################################
## \"$PROGRAM\"  Execution Time:  $h:$m:$s
##########################################################
EOF
    };
} # END_SUB: get_exec_time

# Defining output filenames.
#
sub set_input_file {
	$blast_file = "STDIN";
	if ($ARGV ne '-') {
		$blast_file = $ARGV;
		$stdin_flg = 0;
		-e $blast_file || die "Can't open the file '$blast_file' : $!\n";
	};
}
sub set_output_files {
	$blast_alnQ = $blast_alnS = "STDOUT";
	$to_file_flg && do {
		if ($to_file == 1) {
			$blast_alnQ = "$blast_file.alnQ"; 
			$blast_alnS = "$blast_file.alnS";
		} else {
			$blast_alnQ = "$to_file.alnQ"; 
			$blast_alnS = "$to_file.alnS";	
		};
	};
}
sub set_hsp_file {
	$hsp_file = "STDOUT";
	($to_file_flg && $hsp_too_flg) && do {
		if ($hsp_flg == 1) {
			$hsp_file = "$blast_file.hsp";
		} else {
			$hsp_file = "$hsp_flg.hsp"; 
		};
	};
}	
sub define_files {
	my ($hstr,$qstr,$sstr,$scocf);

	&set_input_file;

	&set_output_files;

	&set_hsp_file;

	$hstr = "\n\#\#     Writing INPUT HSPs has been disabled\n\#\#";
	$hstr = "\n\#\#     HSPs from input \"$blast_file\" : STDOUT\n\#\#" if $hsp_too_flg;
	$hstr = "\n\#\#     HSPs from input \"$blast_file\": $hsp_file\n\#\#"
		if ($to_file_flg && $hsp_too_flg);

	$qstr = "Alignment from QUERY   SRs: $blast_alnQ";
	$qstr = "Writing for QUERY has been disabled (\"-S\" option)." if $only_S_flg; 

	$sstr = "Alignment from SUBJECT SRs: $blast_alnS";
	$sstr = "Writing for SUBJECT has been disabled (\"-Q\" option)." if $only_Q_flg; 

	$scocf = "\n\#\# No SCORE CUTOFF was set.\n\#\#";
	$scocf = "\n\#\# SCORE CUTOFF was set to: $cutoff_flg\n\#\#"
		if ($cutoff_flg ne "-" && $_[5]<$cutoff_flg);

print STDERR <<EOF if $verbose_flg;
##########################################################
##
## $PROGRAM - $VERSION
##
## Reading Input from: $blast_file
##
## Writing Results to following files.
##$hstr
##     $qstr
##
##     $sstr
##$scocf
##########################################################
##
EOF
} # END_SUB: define_files


##############################################################################
##                           MAIN Subroutines                               ##
##############################################################################

# Printing hsp records
#
sub print_hsp {
	my $record;
print STDERR <<EndOfPrt;
##
##########################################################
## $index HSP records were Read: LISTING
##########################################################
EndOfPrt

	foreach my $k (sort { $a <=> $b } keys %hsp) {
		$record = \%{ $hsp{$k} };
print STDERR <<EndOfPrt;
#
# HSP_INDEX: $k
# SCORE    : $record->{SCORE}
# E_VALUE  : $record->{E_VALUE}

  QUERY    : $record->{QUERY}    
  START_Q  : $record->{START_Q}  
  END_Q    : $record->{END_Q}    
  STRAND_Q : $record->{STRAND_Q} 
  FRAME_Q  : $record->{FRAME_Q}  
  SEQ_Q    : $record->{SEQ_Q}

  |  SUBJECT  : $record->{SUBJECT}
  |  START_S  : $record->{START_S}
  |  END_S    : $record->{END_S}
  |  STRAND_S : $record->{STRAND_S}
  |  FRAME_S  : $record->{FRAME_S}
  |  SEQ_S    : $record->{SEQ_S}

#
EndOfPrt
	};
} # END_SUB: print_hsp
#
sub new_hsp {
	$hsp{$index} = {
		QUERY    => $_[0],  SUBJECT  => $_[9],
		START_Q  => $_[3],  START_S  => $_[10],
		END_Q    => $_[4],  END_S    => $_[11],
		STRAND_Q => $_[6],  STRAND_S => $_[15],
		FRAME_Q  => $_[7],  FRAME_S  => $_[17],
		SEQ_Q    => $_[19], SEQ_S    => $_[21],
		SCORE    => $_[5],  E_VALUE  => $_[13],
		INDEX    => $index,
	};
	# Here are defined the auxiliary lists needed to obtain SRs
	!$only_S_flg && do {
		push @{ $project_query[$codeQ] },   [ $hsp{$index}{START_Q}, 
											  $hsp{$index}{START_S}, 
											  $hsp{$index}{SCORE},
											  $index ] ;
		push @{ $project_query[$codeQ] },   [ $hsp{$index}{END_Q},
											  $hsp{$index}{END_S},
											  $hsp{$index}{SCORE},
											  $index ] ;
	};
	!$only_Q_flg && do {
		push @{ $project_subject[$codeS] }, [ $hsp{$index}{START_S},
											  $hsp{$index}{START_Q},
											  $hsp{$index}{SCORE},
											  $index ] ;
		push @{ $project_subject[$codeS] }, [ $hsp{$index}{END_S},
											  $hsp{$index}{END_Q},
											  $hsp{$index}{SCORE},
											  $index ] ;
	};
} # END_SUB: new_hsp

# Read HSPs from file or STDIN, and defining data structures.
#
sub read_HSPs {
print STDERR <<EndOfPrt if $verbose_flg;
##########################################################
## FILTERING INPUT WITH \"$pb_PROG\"
##########################################################
##
EndOfPrt

  READING: {
	  $stdin_flg && do {
		  open2(*QUERY,*INPUT,"$pb_PROG");
		  print INPUT <ARGV>;
		  close(INPUT);
		  last READING;
	  };
	  open(QUERY, "$pb_PROG $blast_file |");
  }

	$no_stdout = 0;
	$hsp_too_flg && do {
		$no_stdout = &open_file_handle($to_file_flg,"$hsp_file");
		print { $FLH } "\#\#\n\#\# HSPs from $blast_file\n\#\#\n";
	};

	$index = 0;
	while (<QUERY>){ 
		next if /^\#|^\s*$/;
		chomp;
		print { $FLH } "@_\n" if $hsp_too_flg;
		split;
		next if ($cutoff_flg ne "-" && $_[5]<$cutoff_flg);
		$_[9] =~ s/\"//og;
		$codeQ = $STR{$_[6]}  + $_[7] ;
		$codeS = $STR{$_[15]} + $_[17];
		&new_hsp;
		$index++;
	}; 

	close($FLH) if $no_stdout;
	close(QUERY);

	&print_hsp if $debug_flg;
} # END_SUB: read_HSPs

# Printing Lists of Projected Points
#
sub print_proj_records {
	my ($w_sr, $array) = @_;
	my ($nx, $nxt, $rcd);
	my $m_sr = $w_sr eq "SUBJECT" ? "Query" : "Subject";
print STDERR <<EndOfPrt;
##
##########################################################
## Projection Lists were Made for $w_sr
##########################################################
EndOfPrt

	foreach my $c (1..$max_frame_num) {
		defined(@{ $array->[$c] }) && do {
			$nx = $#{ $array->[$c] };
			print STDERR "*** $w_sr *** STRAND/FRAME: $FRxST{$c} : ".($nx+1)." points ***\n";
			foreach my $n (0..$nx) {
				$rcd = \@{ $array->[$c][$n] };
print STDERR <<EndOfPrt;
Position: $rcd->[0]  At$m_sr: $rcd->[1]  Score: $rcd->[2]  Index: $rcd->[3]  HSP: $hsp{$rcd->[3]}{INDEX} $hsp{$rcd->[3]}{QUERY}_x_$hsp{$rcd->[3]}{SUBJECT}
EndOfPrt
            };
            $nxt += $nx + 1;
        };
	};
    return $nxt;
} # END_SUB: print_proj_records
#
sub print_project { 
	my ($tfxsq,$tfxss) = (0,0); 
	!$only_S_flg && ( $tfxsq = &print_proj_records(" QUERY ", \@project_query) );
	!$only_Q_flg && ( $tfxss = &print_proj_records("SUBJECT", \@project_subject) );

print STDERR <<EndOfPrt;
##########################################################
## HSPs: $index
## Total Projected Points: $tfxsq (QUERY) x $tfxss (SUBJECT)
##########################################################
EndOfPrt
} # END_SUB: print_project

# Sorting list of projected points for each hsp coords by position
#
sub sort_by_fields {
	my ($w_sr, $array) = @_;
    my $c;
	print STDERR "\#\#\n\#\# $w_sr\n" if $verbose_flg;
 SORT: foreach $c (1..$max_frame_num) {
	  defined(@{ $array->[$c] }) && do {
		  print STDERR "\#\#     SORTING *** STRAND/FRAME: $FRxST{$c} ***\n" if $verbose_flg;
		  $array->[$c] = [ map { $_->[4] }
						   sort { $a->[0] <=> $b->[0] # sorting by position
					                      ||
					   	          $b->[2] <=> $a->[2] # reverse sorting by score
					                      ||
						          $a->[1] <=> $b->[1] # sorting by SUBJECT position
		   			                      ||
							      $a->[3] <=> $b->[3] } # sorting by global index
						   map { [ $_->[0] , $_->[1] , $_->[2] , $_->[3] , $_ ] }
						   @{ $array->[$c] } ];
		  next SORT;
	  }; # defined(@{ $array[$c] })
	  print STDERR "\#\#     SORTING *** STRAND/FRAME: $FRxST{$c} - No Coords Defined ***\n" if $verbose_flg;
  }; # foreach my $c (1..$max_frame_num)
} # END_SUB: sort_by_fields
#
sub sort_projections {
print STDERR <<EndOfPrt if $verbose_flg;
##
##########################################################
## SORTING LISTs of PROJECTED HSPs COORDs by POSITION
##########################################################
EndOfPrt

	!$only_S_flg && &sort_by_fields(" QUERY ", \@project_query);
    !$only_Q_flg && &sort_by_fields("SUBJECT", \@project_subject);

    &print_project if $debug_flg;

} # END_SUB: sort_projections

# Printing SR records
#
sub prt_Q_fullgff {
print { $FLH } <<"EndOfFullGFF";
$srs->{QUERY}\t$PROGRAM\tsr\t$srs->{START_Q}\t$srs->{END_Q}\t$srs->{SCORE}\t$srs->{STRAND_Q}\t$srs->{FRAME_Q}\tTarget \"$srs->{SUBJECT}\"\t$srs->{START_S}\t$srs->{END_S}\tE_value $srs->{E_VALUE}\tStrand $srs->{STRAND_S}\tFrame $srs->{FRAME_S}\t\#Projection $srs->{PROJECTION} \#Seq-Query: $srs->{SEQ_Q} \#Seq-Subject: $srs->{SEQ_S}
EndOfFullGFF
} # END_SUB: prt_Q_fullgff
#
sub prt_S_fullgff {
print { $FLH } <<"EndOfFullGFF";
$srs->{SUBJECT}\t$PROGRAM\tsr\t$srs->{START_S}\t$srs->{END_S}\t$srs->{SCORE}\t$srs->{STRAND_S}\t$srs->{FRAME_S}\tTarget \"$srs->{QUERY}\"\t$srs->{START_Q}\t$srs->{END_Q}\tE_value $srs->{E_VALUE}\tStrand $srs->{STRAND_Q}\tFrame $srs->{FRAME_Q}\t\#Projection $srs->{PROJECTION} \#Seq-Subject: $srs->{SEQ_S} \#Seq-Query: $srs->{SEQ_Q}
EndOfFullGFF
} # END_SUB: prt_S_fullgff
#
sub prt_Q_aplot {
print { $FLH } <<"EndOfAPLOT";
$srs->{QUERY}:$srs->{SUBJECT}\t$PROGRAM\tsr\t$srs->{START_Q}:$srs->{START_S}\t$srs->{END_Q}:$srs->{END_S}\t$srs->{SCORE}\t$srs->{STRAND_Q}:$srs->{STRAND_S}\t$srs->{FRAME_Q}:$srs->{FRAME_S}\t$srs->{INDEX_HSP}:$srs->{INDEX_SR}\t\#E_value $srs->{E_VALUE} \#Projection $srs->{PROJECTION} \#Seq-Query: $srs->{SEQ_Q} \#Seq-Subject: $srs->{SEQ_S}
EndOfAPLOT
} # END_SUB: prt_Q_aplot
# 
sub prt_S_aplot {
print { $FLH } <<"EndOfAPLOT";
$srs->{SUBJECT}:$srs->{QUERY}\t$PROGRAM\tsr\t$srs->{START_S}:$srs->{START_Q}\t$srs->{END_S}:$srs->{END_Q}\t$srs->{SCORE}\t$srs->{STRAND_S}:$srs->{STRAND_Q}\t$srs->{FRAME_S}:$srs->{FRAME_Q}\t$srs->{INDEX_HSP}:$srs->{INDEX_SR}\t\#E_value $srs->{E_VALUE} \#Projection $srs->{PROJECTION} \#Seq-Subject: $srs->{SEQ_S} \#Seq-Query: $srs->{SEQ_Q}
EndOfAPLOT
} # END_SUB: prt_S_aplot
# 
sub prt_pairwise {
	my ($ml,$a,$b,$x,$y,$hsq,$heq,$hss,$hes);
	($a,$b) = ($srs->{QUERY},$srs->{SUBJECT});
	$ml = max(length($a),length($b));
	($a,$b) = (&fill_right($a,$ml," "),&fill_right($b,$ml," "));
	($hsq,$heq,$hss,$hes) = ($srs->{START_Q},$srs->{END_Q},$srs->{START_S},$srs->{END_S});
	$ml = &max(length($hsq),length($heq),length($hss),length($hes));
	($x,$y) = ( &fill_left($hsq,$ml," ")." ".&fill_left($heq,$ml," ")." $srs->{STRAND_Q} $srs->{FRAME_Q}",
				&fill_left($hss,$ml," ")." ".&fill_left($hes,$ml," ")." $srs->{STRAND_S} $srs->{FRAME_S}");

	$a .= " Q $x $srs->{SEQ_Q}\n";
	$b .= " S $y $srs->{SEQ_S}\n";
	$which_proj eq "QUERY" && do {
		print { $FLH } "\#\n$a$b";
		return;
	};
	print { $FLH } "\#\n$b$a";
} # END_SUB: prt_pairwise
#
sub print_SRs {
	my $record;
print STDERR <<EndOfPrt;
##
##########################################################
## $sr_count SRs records were Build: LISTING
##########################################################
EndOfPrt

	foreach my $k (sort { $a <=> $b } keys %sr) {
		$record = \%{ $sr{$k} };
print STDERR <<EndOfSRs
#
#   SR_INDEX: $k 
#
# PROJECTION: $record->{PROJECTION}
#  HSP_INDEX: $record->{INDEX_HSP}
#      SCORE: $record->{SCORE}
#    E_VALUE: $record->{E_VALUE}

  QUERY    : $record->{QUERY}    
  START_Q  : $record->{START_Q}  
  END_Q    : $record->{END_Q}    
  STRAND_Q : $record->{STRAND_Q} 
  FRAME_Q  : $record->{FRAME_Q}  
  SEQ_Q    : $record->{SEQ_Q}

  |  SUBJECT  : $record->{SUBJECT}
  |  START_S  : $record->{START_S}
  |  END_S    : $record->{END_S}
  |  STRAND_S : $record->{STRAND_S}
  |  FRAME_S  : $record->{FRAME_S}
  |  SEQ_S    : $record->{SEQ_S}

#
EndOfSRs
	};
} # END_SUB: print_SRs

# Defining NEW SR records
#
sub new_SR {
	$sr{$sr_count} = {
		PROJECTION => $which_proj,
		QUERY     => $rec->{QUERY},    SUBJECT  => $rec->{SUBJECT},
		START_Q   => $Q_ori,           START_S  => $S_ori,
		END_Q     => $Q_end,           END_S    => $S_end,
		STRAND_Q  => $rec->{STRAND_Q}, STRAND_S => $rec->{STRAND_S},
		FRAME_Q   => $rec->{FRAME_Q},  FRAME_S  => $rec->{FRAME_S},
		SEQ_Q     => $Q_seq,           SEQ_S    => $S_seq,
		SCORE     => $sco,             E_VALUE  => $rec->{E_VALUE},
		INDEX_HSP => $idx,
		INDEX_SR  => $sr_count,
	};
} # END_SUB: new_SR
#
sub get_sr { # GETS SRstart SRend HSPscore HSPindex
	my $seq_ori;
	($ori,$end,$sco,$idx) = @_;
	$rec = \%{ $hsp{$idx} };
	print STDERR "##GREP##*****\n##GREP##********** @_ :: $stv \n##GREP##*****\n" if $debug_flg;

  MKVARS: {
	  $which_proj eq "QUERY" && do {
		  ($Q_ori,$Q_end) = ($ori,$end);
		  $seq_ori = $Q_ori   - $rec->{START_Q};
		  $S_ori   = $seq_ori + $rec->{START_S};
		  $dif     = $end - $ori + 1;
		  $S_end   = $dif - 1 + $rec->{START_S};
		  last MKVARS;
	  }; # else $which_proj eq "SUBJECT" 
	  ($S_ori,$S_end) = ($ori,$end);
	  $seq_ori = $S_ori   - $rec->{START_S};
	  $Q_ori   = $seq_ori + $rec->{START_Q};
	  $dif     = $end - $ori + 1;
	  $Q_end   = $dif - 1 + $rec->{START_Q};
  };
	$Q_seq = substr($rec->{SEQ_Q},$seq_ori,$dif);
	$S_seq = substr($rec->{SEQ_S},$seq_ori,$dif);
	&new_SR($rec);
	
    $srs = \%{ $sr{$sr_count} };
  PRTOUT: {
	  $aln_flg    && (&prt_pairwise, last PRTOUT);
	  $which_proj eq "QUERY" && do {
		  $aplot_flg  && (&prt_Q_aplot, last PRTOUT);
		  &prt_Q_fullgff; last PRTOUT;
	  }; # else $which_proj eq "SUBJECT"
	  $aplot_flg  && (&prt_S_aplot, last PRTOUT);
	  &prt_S_fullgff;
  };
	$sr_ctr{$which_proj}{$stv}++;
	$sr_count++;
} # END_SUB: get_sr

# Building SRs for each Coord on HSPs.
#
sub chkvars {
	print STDERR "##GREP## @_ ## COORD:$last_coord/$current_coord ## SCORE:$last_score/$current_score ## INDEX:$last_index/$current_index\n" if $debug_flg;
} # END_SUB: chkvars
#
sub get_from_stack {
	@stack>0 && do {
		do { # Cleaning stack for all elements already closed.
			($t_score, undef, $t_index) = @{ shift @stack };
		} until (exists($opened{$t_index}) || @stack==0);
		exists($opened{$t_index}) && do { # must be exists($opened{$t_index}) from above do/until.
			($last_coord, $last_score, $last_index) = ($current_coord+1, $t_score, $t_index);
			return;
		};
		############################################################ NEW # $last_coord=$current_coord;
	}; # @stack>0 
} # END_SUB: get_from_stack
#
sub check_GT_closed_score {
	push @stack, [ $current_score, $current_coord, $current_index ];
  CLSC: {
	  # $current_coord>$last_coord && !opened{$current_index} && $current_score>$last_score
	  $current_score>$last_score && do { 
		  &chkvars("#01($op,".@stack.")");
		  &get_sr($last_coord, $current_coord-1, $last_score, $last_index);
		  ($last_coord,$last_score,$last_index) = ($current_coord,$current_score,$current_index);
		  last CLSC;
	  };
	  # $current_coord>$last_coord && !opened{$current_index} && $current_score<$last_score
	  $current_score<$last_score && do {
		  &chkvars("#02($op,".@stack.")");
		  last CLSC;
	  };
	  # $current_coord>$last_coord && !opened{$current_index} && $current_score==$last_score
	  # &get_sr($last_coord, $current_coord, $last_score, $last_index);
	  &chkvars("#25($op,".@stack.")");
	  exists($opened{$last_index}) && do {
		  &chkvars("#03($op,".@stack.")");
		  $current_index == $last_index && do {
			  &chkvars("#27($op,".@stack.")");
			  &get_sr($last_coord, $current_coord, $last_score, $last_index);
			  delete($opened{$current_index}); $op--;
			  ($last_coord,$last_index) = ($current_coord+1,$current_index);
		  };
		  last CLSC;
		  # shift @stack; ########################################################################### NEW
	  };
	  &chkvars("#04($op,".@stack.")");
	  $opened{$current_index} = $n; $op++;
	  ($last_coord,$last_index) = ($current_coord+1,$current_index);
	  last CLSC;
  }; # CLSC
	# Sorting temporary stack for opened-index scores.
	&chkvars("#05($op,".@stack.")");
	my @t_stack = ( map  { $_->[0] }
					sort { $b->[1] <=> $a->[1] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3] }
					map  { [ $_, $_->[0] , $_->[1] , $_->[2] ] } @{ @stack } ) ;
	@stack = @t_stack;
} # END_SUB: check_GT_closed_score
#
sub check_GT_opened_score {
  OPSC: {
	  # $current_coord>$last_coord && opened{$current_index} && $current_score>$last_score
	  $current_score>$last_score && do { 
		  &chkvars("#06($op,".@stack.")");
		  &get_sr($last_coord, $current_coord-1, $last_score, $last_index);
		  ($last_coord,$last_score,$last_index) = ($current_coord,$current_score,$current_index);
		  last OPSC;
	  };
	  # $current_coord>$last_coord && opened{$current_index} && $current_score<$last_score
	  $current_score<$last_score && do {
		  &chkvars("#07($op,".@stack.")");
		  !exists($opened{$last_index}) && do { 
			  &chkvars("#08($op,".@stack.")");
			  ($last_score,$last_index) = ($current_score,$current_index);
			  &get_sr($last_coord, $current_coord, $last_score, $last_index);
			  $last_coord = $current_coord;
		  };
		  &chkvars("#09($op,".@stack.")");
		  delete($opened{$current_index}); $op--;
		  last OPSC;
	  };
	  # $current_coord>$last_coord && opened{$current_index} && $current_score==$last_score
	  &chkvars("#10($op,".@stack.")");
	  $current_index == $last_index && do {
		  &chkvars("#26($op,".@stack.")");	  
		  &get_sr($last_coord, $current_coord, $last_score, $last_index);
		  # shift @stack; ########################################################################### NEW
		  ($last_coord,$last_score,$last_index) = ($current_coord,$current_score,$current_index);
		  &get_from_stack;
	  };
	  delete($opened{$current_index}); $op--;
  }; # OPSC
} # END_SUB: check_GT_opened_score
#
sub check_EQ_coord {
	# $current_coord==$last_coord && $current_score==$last_score
	$current_score==$last_score &&  do { ### replace eq by ==
		&chkvars("#11($op,".@stack.")");
	  CURRIDX: {
		  # $current_index == $last_index
		  $current_index == $last_index && do {
			  &chkvars("#12($op,".@stack.")");
			  &get_sr($last_coord, $current_coord, $last_score, $last_index);
			  last CURRIDX;
		  };
		  # ELSE $current_index != $last_index
		  &chkvars("#13($op,".@stack.")");
		  push @stack, [ $current_score, $current_coord, $current_index ];
	  };
	};
	# $current_coord==$last_coord && $current_score>$last_score
	$current_score>$last_score && do { 
		&chkvars("#14($op,".@stack.")");
		($last_score,$last_index) = ($current_score,$current_index);
	};
# General stuff when $current_coord==$last_coord
	# !exists($opened{$current_index}) && $current_coord==$last_coord
	&chkvars("#15($op,".@stack.")");
	!exists($opened{$current_index}) && do {
		&chkvars("#16($op,".@stack.")");
		$opened{$current_index} = $n; $op++;
		return;
	};
	# ELSE exists($opened{$current_index}) && $current_coord==$last_coord
	&chkvars("#17($op,".@stack.")");
	delete($opened{$current_index}); $op--;
} # END_SUB: check_EQ_coord
#
sub check_coords {
	# $current_coord<$last_coord
	$current_coord<$last_coord && do {
		&chkvars("#18($op,".@stack.")");
		&get_from_stack;
		exists($opened{$current_index}) && do { ############################################### NEW
			&chkvars("#19($op,".@stack.")");
			delete($opened{$current_index}); $op--; ########################################### NEW
		}; #################################################################################### NEW
		return; ############################################################################### NEW
	};
	# $current_coord>$last_coord
	$current_coord>$last_coord && do {
		&chkvars("#20($op,".@stack.")");
		# !exists($opened{$current_index}) && $current_coord>$last_coord
		!exists($opened{$current_index}) && do {
			&chkvars("#21($op,".@stack.")");
			&check_GT_closed_score;
			$opened{$current_index} = $n; $op++;
			return;
		};
		# ELSE exists($opened{$current_index}) && $current_coord>$last_coord
		&chkvars("#22($op,".@stack.")");
		&check_GT_opened_score;
		return;
	}; # $current_coord>$last_coord
	# ELSE $current_coord==$last_coord
	&chkvars("#23($op,".@stack.")");
	&check_EQ_coord;
} # END_SUB: check_coords
#
sub test_projection {
# 	$tfxs{$w_sr} = $#{ $lines->[$c] };
	($w_sr, $lines) = @_;
	my $lnc;
  MAIN: foreach $c (1..$max_frame_num) {
	  my $yc = 0;
	  ($which_proj = $w_sr) =~ s/\s//og; #QUERY/SUBJECT
	  $stv = $FRxST{$c};
	  $sr_ctr{$which_proj}{$stv} = 0;
	  defined(@{ $lines->[$c] }) && do {
		  $lnc = $#{ $lines->[$c] };
		  # print STDERR "*** $w_sr *** SRs on STRAND/FRAME: $stv *** (checking ".&fill_left(($lnc+1),3," ")." points) ***\n" if $verbose_flg;
		  print STDERR "##GREP## ".("*"x73)."\n##GREP## *** $w_sr *** SRs on STRAND/FRAME: $stv *** (checking ".&fill_left(($lnc+1),3," ")." points) ***\n##GREP## ".("*"x73)."\n" if $debug_flg;
		  $op = 0;
		LINES: for $n (0..$lnc) {
			# &prt_progress(++$yc);
			$rcd = \@{ $lines->[$c][$n] };
			$current_coord = $rcd->[0];
			$current_score = $rcd->[2];
			$current_index = $rcd->[3];

			# $n>0 && $op>0 :: n_ary element to be checked && HSPs opened
			($op>0) && do { 
				&chkvars(("-"x73)."\n##GREP## #24($op,".@stack.")");
				&check_coords;
				next LINES;
			};
			# ELSE $n==0 || $op==0 :: First element to be checked || No HSPs opened
			&chkvars("#25($op,".@stack.")");
			$opened{$current_index} = $n; $op++;
			@stack = ();
			push @stack, [ $current_score, $current_coord, $current_index ];
			($last_coord,$last_score,$last_index) = ($current_coord,$current_score,$current_index);

		} # for :LINES: $n (0..$lnc) 
		  # &prt_foeprg($yc);
		  print STDERR "##GREP## \n" if $debug_flg;
		  next MAIN;
	  }; # defined(@{ $lines[$c] })
	  # print STDERR "*** $w_sr *** SRs on STRAND/FRAME: $stv *** No Coords Defined ***\n" if $verbose_flg;	
	  print STDERR "##GREP## *** $w_sr *** SRs on STRAND/FRAME: $stv *** No Coords Defined ***\n" if $debug_flg;	
  }; # foreach :MAIN: $c (1..$max_frame_num)
} # END_SUB: test_projection 
#
sub summary_of_SRs {
print STDERR <<EndOfPrt;
##
################################################ SUMMARY #
EndOfPrt

	my @g = ();
    foreach my $p (values %{ $sr_ctr{QUERY} }, values %{ $sr_ctr{SUBJECT} }) { 
		push @g, length($p);
	};
    my $f = &max(@g);
	my $hd = 0;
    foreach my $r ("QUERY", "SUBJECT") {
		print STDERR "\#\#\n\#\#  $r\n";
		my $hc = 0;
		foreach my $j (sort keys %{ $sr_ctr{$r} }) {
			$hc += $sr_ctr{$r}{$j};
			print STDERR "\#\#    STRAND/FRAME:  $j  ->  ".&fill_left($sr_ctr{$r}{$j},$f," ")."  SRs found.\n";
		};
		print STDERR "\#\#    ------------ SUBTOTAL: ".&fill_left($hc,$f," ")."  SRs for $r.\n";
		$hd += $hc;
	};
		print STDERR "\#\#\n\#\#".&fill_left(" TOTAL SRs: $hd #",56,"\#")."\n";
} # END_SUB: sumary_of_SRs

# Computing SRs from projected HSP coords.
#
sub obtain_SRs {
print STDERR <<EndOfPrt if $verbose_flg;
##
##########################################################
## OBTAINING SRs from PROJECTED HSPs COORDs
##########################################################
EndOfPrt
	
    $sr_count = 0;
    !$only_S_flg && do {
		$no_stdout = &open_file_handle($to_file_flg,$blast_alnQ);
		print { $FLH } "\#\#\n\#\# SRs from $blast_file : QUERY Projection\n\#\#\n";
		&test_projection( " QUERY ", \@project_query );
		close($FLH) if $no_stdout;
	};
	
    !$only_Q_flg && do {
		$no_stdout = &open_file_handle($to_file_flg,$blast_alnS);
		print { $FLH } "\#\#\n\#\# SRs from $blast_file : SUBJECT Projection\n\#\#\n";
		!$only_Q_flg && &test_projection( "SUBJECT", \@project_subject );
		close($FLH) if $no_stdout;
	};
	
    &print_SRs if $debug_flg;	
	
	&summary_of_SRs if $verbose_flg;

} # END_SUB: obtain_SRs


##############################################################################
##                               MAIN LOOP                                  ##
##############################################################################

unshift(@ARGV,'-') unless @ARGV;

while ($ARGV = shift @ARGV){
	
	&define_files;
	
	&read_HSPs;
	
	&sort_projections;

	&obtain_SRs;

};

&get_exec_time(time);

exit(0);

# GetSRsAln -Db -H -W -- U66875_U34801.tbx 2>&1 | grep "##GREP##" | perl -ne 's/^\#\#GREP\#\#//o; print +(" " x (4 - length($.))).($.)." :".$_;' | head -350 | enscript -rf Courier5 --columns=3 --margins=25:25:-50:25 
