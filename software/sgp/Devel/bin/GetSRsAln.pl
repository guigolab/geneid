#!/usr/bin/perl
#
# GetSRsAln.pl - Obtaining Similarity Regions and its sequence from HSPs.
#
# $Id: GetSRsAln.pl,v 1.2 2000-08-10 02:15:35 jabril Exp $
#

my $PROGRAM = "GetSRsAln.pl";
my $VERSION = '$Version:$ ';
my $Start = time;

use strict; # 'refs';
use IPC::Open2;
use Getopt::Long;


##############################################################################
##                       Getting Comand-line Options                        ##
##############################################################################

Getopt::Long::Configure("bundling","pass_through");

my ( $verbose_flg, $aln_flg, $bit_flg, $ids_flg,
	 $only_Q_flg, $only_S_flg, $to_file,
	 $help_flg, $prt_aln_flg, $debug_flg, $stdin_flg
	 ) = (0,0,0,0,0,0,0,0,0,0,1);

GetOptions( "Q|query-only"     => \$only_Q_flg  ,
			"S|subject-only"   => \$only_S_flg  ,
			"A|alignment"      => \$prt_aln_flg , # show only alignment not GFF
			"b|bit-score"      => \$bit_flg     ,
			"i|identity-score" => \$ids_flg     ,
			"W|write-to-file"  => \$to_file     ,
			"v|verbose"        => \$verbose_flg ,
			"D|debugging"      => \$debug_flg   ,
 			"h|help|\?"        => \$help_flg    ,
			);

my $ARGV = my $scoreopt = my $vrbopt = "";
$verbose_flg = 1 if $debug_flg;
$vrbopt = " -v" if $verbose_flg;
# $aln_flg = 0 if ($bit_flg || $ids_flg);
$scoreopt = " -i" if $ids_flg;
$scoreopt = " -b" if $bit_flg;
$VERSION =~ s/\$//og;

$help_flg && &prt_help; 


##############################################################################
##                          Global Variables                                ##
##############################################################################

my $parseblast = "parseblast"; # $blastgff = "blast2gff -g"
my $pb_QUERY   = "$parseblast$vrbopt$scoreopt -nFQ";
                 # 'n'o comments.
                 # 'F'ullgff better then 'G'ff.
                 # Append se'Q'uence to GFF.
                 # 'bit' score from BLAST .
# my $pb_SUBJECT = "$parseblast$vrbopt$scoreopt -nGbS";
# my $pb_ALN     = "$parseblast$vrbopt$scoreopt -nPW";
my ( $blast_file, $blast_alnQ, $blast_alnS, @data );

my ( %hsp, %sr, @project_query, @project_subject, $rcd, $last_rcd,
     @stack, %opened, $coord, $score, $index,
	 $current_coord, $last_coord, $current_score, $last_score,
	 $current_index, $last_index, $index, $sr_count, # $stack_index,
	 $op, $n, $c, $stv, $t_score, $t_index,
	 %tfxs, $s, $e, $dif, $sQ, $sS, 
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
  (Don't you have your ~/bin dir in the path?)............ ;^) 
  
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
    -A, --alignment      : print pairwise alignment for each SR.
    -b, --bit-score      : set <score> field to Bits (default Alignment Score).
    -i, --identity-score : set <score> field to Identities (default Alignment).
    -W, --write-to-file  : write output to separate files
                                 + for QUERY:   "<input_name>.alnQ"
                                 + for SUBJECT: "<input_name>.alnS"
                           you can provide "<input_name>" as parameter for
                           this option, as example, if you provide "-W results"
                           you send output to "results.alnQ" and "results.alnS".
    -v, --verbose        : warnings sent to <STDERR>.
    -D, --debuggingv     : extended report for debugging sent to <STDERR>.
    -h, --help           : show this help pages.

BUGS:    Report any problem to: abril\@imim.es

AUTHOR:  $PROGRAM is under GNU-GPL (C) 2000 - Josep F. Abril

EndOfHelp
        close(HELP);
exit(1);
}
 
# Reporting IN/OUT progress.
#
sub prt_progress {
    $verbose_flg && do {
                print STDERR ".";
                (($_[0] % 50) == 0) && print STDERR "[".&fill_left($_[0],6,"0")."]\n";
        };
}
#
sub prt_foeprg {
    $verbose_flg && ((($_[0] % 50) != 0) && print STDERR "[".&fill_left($_[0],6,"0")."]\n" );
}

# Get a fixed length string from a given string and filling char/s.
#
sub fill_right { $_[0].($_[2] x ($_[1] - length($_[0]))) }
#
sub fill_left  { ($_[2] x ($_[1] - length($_[0]))).$_[0] }

# Timing.
#
sub get_exec_time {
	$verbose_flg && do {
        my $End = time;
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
sub define_files {
	$blast_file = "STDIN";
	if ($ARGV ne '-') {
		$blast_file = $ARGV;
		$stdin_flg = 0;
		-e $blast_file || die "Can't open the file '$blast_file' : $!\n";
	}
	$blast_alnQ = "$blast_file.alnQ"; 
	$blast_alnS = "$blast_file.alnS";

print STDERR <<EOF if $verbose_flg;
##########################################################
##
## $PROGRAM
##
## Reading Input from: $blast_file
##
## Writing Results to following files.
##   Alignment from QUERY    SRs: $blast_alnQ
##   Alignment from SUBJECTS SRs: $blast_alnS
##
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

# Read HSPs from file or STDIN, and defining data structures.
#
sub read_HSPs {
	# open(ALN,     "$pb_ALN     $blast_file | ");
	# open(SUBJECT, "$pb_SUBJECT $blast_file | $blastgff |");
	# close(SUBJECT);
	# close(ALN);

  READING: {
	  $stdin_flg && do {
		  open2(*QUERY,*INPUT,"$pb_QUERY");
		  print INPUT <ARGV>;
		  close(INPUT);
		  last READING;
	  };
	  open(QUERY, "$pb_QUERY $blast_file |");
  }

    my ( $codeQ, $codeS, $id, $rec); # , @tsq, @teq, @tss, @tes);
	$index = 0;
	while (<QUERY>){ 
		next if /^\#|^\s*$/;
		chomp;
		split;
		$_[9] =~ s/\"//og;
		$codeQ = $STR{$_[6]}  + $_[7] ;
		$codeS = $STR{$_[15]} + $_[17];
        $id = $index; # &fill_left($index,8,"0");
		$hsp{$id} = {
			QUERY    => $_[0],  SUBJECT  => $_[9],
			START_Q  => $_[3],  START_S  => $_[10],
			END_Q    => $_[4],  END_S    => $_[11],
			STRAND_Q => $_[6],  STRAND_S => $_[15],
			FRAME_Q  => $_[7],  FRAME_S  => $_[17],
			SEQ_Q    => $_[19], SEQ_S    => $_[21],
			SCORE    => $_[5],
			INDEX    => $id,
		};
		
		$rec = \%{ $hsp{$id} };
		push @{ $project_query[$codeQ] },   [ $hsp{$id}{START_Q}, 
											  $hsp{$id}{START_S}, 
											  $hsp{$id}{SCORE},
											  $id, $rec ] ;
		push @{ $project_query[$codeQ] },   [ $hsp{$id}{END_Q},
											  $hsp{$id}{END_S},
											  $hsp{$id}{SCORE},
											  $id, $rec ] ;
		push @{ $project_subject[$codeS] }, [ $hsp{$id}{START_S},
											  $hsp{$id}{START_Q},
											  $hsp{$id}{SCORE},
											  $id, $rec ] ;
		push @{ $project_subject[$codeS] }, [ $hsp{$id}{END_S},
											  $hsp{$id}{END_Q},
											  $hsp{$id}{SCORE},
											  $id, $rec ] ;

		$index++;
	}; 

	close(QUERY);

	&print_hsp if $verbose_flg;
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
Position: $rcd->[0]  At$m_sr: $rcd->[1]  Score: $rcd->[2]  Index: $rcd->[3]  HSP: ${$rcd->[4]}{INDEX} ${$rcd->[4]}{QUERY}_x_${$rcd->[4]}{SUBJECT}
EndOfPrt
            };
            $nxt += $nx + 1;
        };
	};
    return $nxt;
} # END_SUB: print_proj_records
#
sub print_project { 

	my $tfxsq = &print_proj_records(" QUERY ", \@project_query);
	my $tfxss = &print_proj_records("SUBJECT", \@project_subject);

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
 SORT: foreach $c (1..$max_frame_num) {
	  defined(@{ $array->[$c] }) && do {
		  print STDERR "*** $w_sr *** SORTING STRAND/FRAME: $FRxST{$c} ***\n" if $verbose_flg;
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
	  print STDERR "*** $w_sr *** SORTING STRAND/FRAME: $FRxST{$c} - No Coords Defined ***\n" if $verbose_flg;
  }; # foreach my $c (1..$max_frame_num)
} # END_SUB: sort_by_fields
#
sub sort_projections {
print STDERR <<EndOfPrt;
##
##########################################################
## SORTING LISTs of PROJECTED HSPs COORDs by POSITION
##########################################################
##
EndOfPrt

	&sort_by_fields(" QUERY ", \@project_query);
    &sort_by_fields("SUBJECT", \@project_subject);

    &print_project if $verbose_flg;
} # END_SUB: sort_projections

# Printing SR records
#
sub print_SRs {
print STDERR <<EndOfPrt;
##
##########################################################
## $sr_count SRs records were Build: LISTING
##########################################################
EndOfPrt

	foreach my $k (sort { $a <=> $b } keys %sr) {
print STDERR <<EndOfSRs
#
# SR_INDEX: $k
# SCORE   : $sr{$k}{SCORE}

  QUERY    : $sr{$k}{QUERY}    
  START_Q  : $sr{$k}{START_Q}  
  END_Q    : $sr{$k}{END_Q}    
  STRAND_Q : $sr{$k}{STRAND_Q} 
  FRAME_Q  : $sr{$k}{FRAME_Q}  
  SEQ_Q    : $sr{$k}{SEQ_Q}

  |  SUBJECT  : $sr{$k}{SUBJECT}
  |  START_S  : $sr{$k}{START_S}
  |  END_S    : $sr{$k}{END_S}
  |  STRAND_S : $sr{$k}{STRAND_S}
  |  FRAME_S  : $sr{$k}{FRAME_S}
  |  SEQ_S    : $sr{$k}{SEQ_S}

#
EndOfSRs
	};
} # END_SUB: print_SRs

# Defining NEW SR records
#
sub new_SR {
	my @lstrcd = @_;
	$sr_count++;
	$sr{$last_index} = {
		QUERY    => $lstrcd[4]{QUERY},    SUBJECT  => $lstrcd[4]{SUBJECT},
		START_Q  => $s,                   START_S  => $lstrcd[4]{START_S}-$dif+1,
		END_Q    => $e,                   END_S    => $lstrcd[4]{START_S}+$dif-1,
		STRAND_Q => $lstrcd[4]{STRAND_Q}, STRAND_S => $lstrcd[4]{STRAND_S},
		FRAME_Q  => $lstrcd[4]{FRAME_Q},  FRAME_S  => $lstrcd[4]{FRAME_S},
		SEQ_Q    => $sQ,                  SEQ_S    => $sS,
		SCORE    => $last_score,
		INDEX    => $last_index,
	};
} # END_SUB: new_SR
sub prt_sr {
	print STDOUT "@_ :: $stv\n";
	print STDERR "\n********************* @_ :: $stv ## COORD:$last_coord/$current_coord ## SCORE:$last_score/$current_score ## INDEX:$last_index/$current_index\n";
}
#							  ($s, $e) = ($last_coord, ($current_coord - 1));
#							  $dif = ($e-$s)+1;
#							  $sQ = substr($last_rcd[4]{SEQ_Q}, ($last_rcd[4]{START_Q}), $dif );
#							  $sS = substr($last_rcd[4]{SEQ_S}, ($last_rcd[4]{START_Q}), $dif );
#							  &new_SR(@rcd);

# Building SRs for each Coord on HSPs.
#
sub chkvars {
	print STDERR "@_";
}
#
sub get_from_stack {
	@stack>0 && do {
		&chkvars("#19($op,".@stack.")\n");
		do { # Cleaning stack for all elements already closed.
			($t_score, undef, $t_index) = @{ shift @stack };
			&chkvars("#--->($op,".@stack.",$t_score,$t_index,".(exists($opened{$t_index})?"TRUE":"FALSE").")\n");
		} until (exists($opened{$t_index}) || @stack==0);
		&chkvars("#20($op,".@stack.",$t_score,$t_index,".(exists($opened{$t_index})?"TRUE":"FALSE").")");
		exists($opened{$t_index}) && do { # must be exists($opened{$t_index}) from above do/until.
			&chkvars("#21($op,".@stack.")");
			($last_coord, $last_score, $last_index) = ($current_coord+1, $t_score, $t_index);
			return;
		};
		$last_coord=$current_coord;
	}; # @stack>0 
} # END_SUB: get_from_stack
#
sub check_GT_closed_score {
	push @stack, [ $current_score, $current_coord, $current_index ];
  CLSC: {
	  # $current_coord>$last_coord && !opened{$current_index} && $current_score>$last_score
	  $current_score>$last_score && do { 
		  &chkvars("#11($op,".@stack.")");
		  &prt_sr($last_coord, $current_coord-1, $last_score, $last_index);
		  ($last_coord,$last_score,$last_index) = ($current_coord,$current_score,$current_index);
		  last CLSC;
	  };
	  # $current_coord>$last_coord && !opened{$current_index} && $current_score<$last_score
	  $current_score<$last_score && do {
		  &chkvars("#12($op,".@stack.")");
		  last CLSC;
	  };
	  # $current_coord>$last_coord && !opened{$current_index} && $current_score==$last_score
	  &chkvars("#13($op,".@stack.")");
	  # &prt_sr($last_coord, $current_coord, $last_score, $last_index);
	  exists($opened{$last_index}) && do {
		  &chkvars("#14($op,".@stack.")");
		  &prt_sr($last_coord, $current_coord, $last_score, $last_index);
	  };
	  ($last_coord,$last_index) = ($current_coord+1,$current_index);
  }; # CLSC
	# Sorting temporary stack for opened-index scores.
	my @t_stack = ( map  { $_->[2] }
					sort { $b->[0] <=> $a->[0] || $a->[1] <=> $b->[1] }
					map  { [ $_->[0] , $_->[1] , $_ ] } @{ @stack } ) ;
	@stack = @t_stack;
} # END_SUB: check_GT_closed_score
#
sub check_GT_opened_score {
  OPSC: {
	  # $current_coord>$last_coord && opened{$current_index} && $current_score==$last_score
	  $current_score>$last_score && do { 
		  &chkvars("#15($op,".@stack.")");
		  &prt_sr($last_coord, $current_coord-1, $last_score, $last_index);
		  ($last_coord,$last_score,$last_index) = ($current_coord,$current_score,$current_index);
		  last OPSC;
	  };
	  # $current_coord>$last_coord && opened{$current_index} && $current_score<$last_score
	  $current_score<$last_score && do {
		  &chkvars("#16($op,".@stack.")");
		  !exists($opened{$last_index}) && do { 
			  &chkvars("#17($op,".@stack.")");
			  ($last_score,$last_index) = ($current_score,$current_index);
			  &prt_sr($last_coord+1, $current_coord, $last_score, $last_index);
			  $last_coord = $current_coord;
		  };
		  delete($opened{$current_index}); $op--;
		  last OPSC;
	  };
	  # $current_coord==$last_coord && opened{$current_index} && $current_score==$last_score
	  &chkvars("#18($op,".@stack.")"); # $current_score>$last_score not possible if last_score is max.
	  &prt_sr($last_coord, $current_coord, $last_score, $last_index);
	  delete($opened{$current_index}); $op--;
	  ($last_coord,$last_score,$last_index) = ($current_coord,$current_score,$current_index);
	  &get_from_stack;
  }; # OPSC
} # END_SUB: check_GT_opened_score
#
sub check_coords {
	$current_coord<$last_coord && &get_from_stack;
	# $current_coord>$last_coord
	$current_coord>$last_coord && do {
		# !exists($opened{$current_index}) && $current_coord>$last_coord
		&chkvars("#03($op,".@stack.")");
		!exists($opened{$current_index}) && do {
			&chkvars("#04($op,".@stack.")");
			&check_GT_closed_score;
			$opened{$current_index} = $n; $op++;
			return;
		};
		# ELSE exists($opened{$current_index}) && $current_coord>$last_coord
		&chkvars("#05($op,".@stack.")");
		&check_GT_opened_score;
		return;
	}; # $current_coord>$last_coord
	# ELSE $current_coord==$last_coord
	&chkvars("#06($op,".@stack.")");
	($current_score==$last_score && $current_index eq $last_index) &&  do { ### replace eq by ==
		&chkvars("#07($op,".@stack.")");
		&prt_sr($last_coord, $current_coord, $last_score, $last_index);
	};
	$current_score>$last_score && do { 
		&chkvars("#08($op,".@stack.")");
		($last_score,$last_index) = ($current_score,$current_index);
	};
	# !exists($opened{$current_index}) && $current_coord==$last_coord
	!exists($opened{$current_index}) && do {
		&chkvars("#09($op,".@stack.")");
	    $opened{$current_index} = $n; $op++;
		return;
	};
	# ELSE exists($opened{$current_index}) && $current_coord==$last_coord
	&chkvars("#10($op,".@stack.")");
	delete($opened{$current_index}); $op--;
} # END_SUB: check_coords
#
sub test_projection {
# 	$tfxs{$w_sr} = $#{ $lines->[$c] };
	my ($w_sr, $lines) = @_;
	my $lnc;
  MAIN: foreach $c (1..$max_frame_num) {
	  my $yc = 0;
	  defined(@{ $lines->[$c] }) && do {
		  $stv = $FRxST{$c};
		  $lnc = $#{ $lines->[$c] };
		  print STDERR "\n*** $w_sr *** SRs on STRAND/FRAME: $stv *** (checking ".&fill_left(($lnc+1),3," ")." points) ***\n" if $verbose_flg;
		  $op = 0;
		LINES: for $n (0..$lnc) {
			# &prt_progress(++$yc);
			$rcd = \@{ $lines->[$c][$n] };
			$current_coord = $rcd->[0];
			$current_score = $rcd->[2];
			$current_index = $rcd->[3];

			&chkvars("\n#00($op,".@stack.") Pos: $last_coord/$current_coord Sco:$last_score/$current_score Ind:$last_index/$current_index\n");
	  
			# $n>0 && $op>0 :: n_ary element to be checked && HSPs opened
			($op>0) && do { 
				&chkvars("#02($op,".@stack.")"); # %opened has 1 or more 
				&check_coords;
				&chkvars("\n");
				# foreach $kk (0..$#stack) { print STDERR &fill_right("#####STACK",$kk+12,">").$stack[$kk][0]." ".$stack[$kk][1]." ".$stack[$kk][2]."\n"; };
				next LINES;
			};
			# ELSE $n==0 || $op==0 :: First element to be checked || No HSPs opened
			&chkvars("#01($op,".@stack.")"); # $op>0
			$opened{$current_index} = $n; $op++;
			@stack = ();
			push @stack, [ $current_score, $current_coord, $current_index ];
			($last_coord,$last_score,$last_index) = ($current_coord,$current_score,$current_index);

		} # for :LINES: $n (0..$lnc) 
		  &chkvars("\n");
		  # &prt_foeprg($yc);
		  next MAIN;
	  }; # defined(@{ $lines[$c] })
	  print STDERR "*** $w_sr *** SRs on STRAND/FRAME: $stv *** No Coords Defined ***\n" if $verbose_flg;	
  }; # foreach :MAIN: $c (1..$max_frame_num)
} # END_SUB: test_projection 

# Computing SRs from projected HSP coords.
#
sub obtain_SRs {
print STDERR <<EndOfPrt;
##
##########################################################
## OBTAINING SRs from PROJECTED HSPs COORDs
##########################################################
##
EndOfPrt
	
    $sr_count = 0;
    &test_projection( " QUERY ", \@project_query );
    &test_projection( "SUBJECT", \@project_subject );
	
#	&print_SRs if $verbose_flg;	
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

&get_exec_time;

exit(0);

# gawk '{printf "%10s %10s %4s %s\n", $4,$5,$6,$7$8}' tt3 | enscript -1f Courier9 
# time ../GetSRsAln.pl -v TBLASTX.ex > tt 2> k 

