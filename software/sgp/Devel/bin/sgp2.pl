#!/usr/bin/perl -w
#

=head1 B<sgp2> : $Id: sgp2.pl,v 1.2 2000-10-06 17:59:25 jabril Exp $

=cut Back to the compiler...

my $PROGRAM = "sgp2";
my @tmp_ver = split / +/, ' $Id: sgp2.pl,v 1.2 2000-10-06 17:59:25 jabril Exp $ ';
my $VERSION = "v$tmp_ver[3] [$tmp_ver[4] $tmp_ver[5] $tmp_ver[7]]";
my $Start = time;

=head1 Development Options

=over 4

=item C<perl -w> 

Prints all sorts of useful and interesting warning messages at compile time. 

=item C<use strict;>

Restrict unsafe constructs, like attempting to use missed symbolic references ('I<refs>'), undeclared variables ('I<vars>') or not predeclared subroutine ('I<subs>').

=back

=cut Back to the compiler...

use strict;

=head1 Self Documenting



=cut Back to the compiler...

use Pod::Text;

=head1 Reading Command Line Options

=over 4

=item C<use Getopt::Long;>

The C<Getopt::Long> module implements an extended getopt function called C<GetOptions()>. This function adheres to the POSIX syntax for command line options, with GNU extensions. In general, this means that options having long names instead of single letters are introduced with a double dash "--". 

=item C<Getopt::Long::Configure qw(> I<bundling> I<pass_through> C<);>

C<GetOptions> can be configured by calling subroutine C<Getopt::Long::Configure>. This subroutine takes a list of quoted strings, each specifying a configuration option to be set. Options can be reset by prefixing with C<no_>. 

=over 4

=item I<bundling>

Support for bundling of command line options, as was the case with the more traditional single-letter approach (introduced with a single dash "-"), is provided but not enabled by default. 

=item I<pass_through>

Unknown options are passed through in @ARGV instead of being flagged as errors. This makes it possible to write wrapper scripts that process only part of the user supplied options, and passes the remaining options to some other program.

=back

=back

=cut Back to the compiler...

use Getopt::Long;
Getopt::Long::Configure qw/ bundling pass_through /;

=head2 C<Which_Options()>

This function parses input options, checking whether files exist, 


=cut Back to the compiler...

sub Which_Options() {

	my ($help_flg);
	
	GetOptions( 
				"1"        => \$seq1         , # seqfile_1
				"2"        => \$seq2         , # seqfile_2
				"g"        => \$geneid_opt   , # geneid options      
				"P"        => \$geneid_param , # geneid parameter file 
				"o"        => \$blast_opt    , # tblastx options 
				"c"        => \$score_cutoff , # tblastx score cutoff 
				"s"        => \$ , # shrink hsp's by
				"t"        => \$ , # read tblastx from file
				"f"        => \$ , # read HSP files in directory
				"k"        => \$ , # intermediate filename
				"p"        => \$ps_output    , # postscript output 
				"v"        => \$verbose_flg  , # verbose    
				"h|help|?" => \$help_flg     , 
				);
	
	&prt_Help if $help_flg;

}; # sub Which_Options

=head2 C<prt_Help()>



=cut Back to the compiler...

sub prt_Help() {
	open(HELP, "| more");
	print HELP <<"EndOfHelp";
PROGRAM:  $PROGRAM $VERSION

NAME:
    $PROGRAM - Improving Gene Prediction with Sinteny.

SYNOPSIS:
    $PROGRAM [-hv] [-o \'options\'] [-g \'options\'] \
             [-P filename] [-p filename] [-k filename] \
             [-c value] [-s value] -1 seqfile_1 -2 seqfile_2

DESCRIPTION:

OPTIONS:
 
  -1 seqfile_1   : input file for first species.
  -2 seqfile_2   : input file for second species.
  -g             : geneid options
  -o             : tblastx options
  -c value       : tblastx score cuttof
  -s value       : shrink hsp\'s by value
  -t filename    : read tblastx file
     -f prefix   : read hsp gff files with in directory
                   prefix and extension .hsp-rs
  -k prefix      : keep intermediate files with prefix
  -p filename    : ps output in filename file 
  -P filename    : geneid parameter file
  -v             : verbose mode
  -h             : produces this message

FILES:

DIAGNOSTICS:

REQUIRES:

BUGS:
    Report any problem to: <jabril\@imim.es>

AUTHORS:
    Roderic Guigo  <rguigo\@imim.es>
    Josep F. Abril <jabril\@imim.es>

    $PROGRAM is under GNU-GPL (C) 2000

EndOfHelp
	close(HELP);
	exit(1);
} # sub prt_Help

=head1 Main Loop



=cut Back to the compiler...

	&Which_Options();
