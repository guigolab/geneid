#!/usr/bin/perl -w
#

=head1 B<sgp2> : $Id: sgp2.pl,v 1.1 2000-10-05 16:05:33 jabril Exp $

=cut Back to the compiler...

my $PROGRAM = "sgp2";
my $VERSION = ' $Revision: 1.1 $:$Date: 2000-10-05 16:05:33 $:$State: Exp $ ';
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
Getopt::Long::Configure qw( bundling pass_through );

=head2 C<Which_Options()>



=cut Back to the compiler...

sub Which_Options() {

	my ($ret,$help_flg);
	
	$ret = GetOptions( 
					   "h|help|?"         => \$help_flg    , 
					   );

	

};

