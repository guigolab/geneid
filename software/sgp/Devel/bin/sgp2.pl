#!/usr/bin/perl -w
#line 73 "./sgp2.nw"
use strict;

my $PROGRAM = "sgp2";
my @tmp_ver = split / +/, ' $Id: sgp2.pl,v 1.3 2000-10-09 17:59:59 jabril Exp $ ';
my $VERSION = "v$tmp_ver[3] [$tmp_ver[4] $tmp_ver[5] $tmp_ver[7]]";
my $Start = time;

#line 124 "./sgp2.nw"
use Getopt::Long;
Getopt::Long::Configure qw/ bundling pass_through /;
#line 166 "./sgp2.nw"
use Pod::Text;
#line 127 "./sgp2.nw"
# GetOptions Variables
my ( $seq1,$seq2,$geneid_opt,$geneid_param,$blast_opt,$score_cutoff,
     $shrink,$tbx,$hsp,$ofn,$ps_output,$verbose_flg,$help_flg );


#line 169 "./sgp2.nw"
# Prints help 
sub prt_Help() {
    my $tmp_pod_file = "/tmp/sgp.pod";
    open(KI, "> $tmp_pod_file");
    while (<DATA>) {
        s/\$PROGRAM/$PROGRAM/g ;
        s/\$VERSION/$tmp_ver[3]/g ;
        print KI $_ ;
    };
    close(KI);
    pod2text($tmp_pod_file);
    unlink($tmp_pod_file) or die "Can't delete $tmp_pod_file: $!\n";
    exit(1);
} # sub prt_Help
#line 196 "./sgp2.nw"
# Reporting IN/OUT progress.
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
sub fill_right { $_[0].($_[2] x ($_[1] - length($_[0]))) }
sub fill_left  { ($_[2] x ($_[1] - length($_[0]))).$_[0] }

# returns the max value from input array
sub max { my ($z) = shift @_; my $l; foreach $l (@_) { $z = $l if $l > $z ; }; $z; } 

# Timing.
sub get_exec_time {
    $verbose_flg && do {
        my $End = $_[0];
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
#line 131 "./sgp2.nw"
sub Which_Options() {
    GetOptions( 
                "1"        => \$seq1         , # seqfile_1
                "2"        => \$seq2         , # seqfile_2
                "g"        => \$geneid_opt   , # geneid options      
                "P"        => \$geneid_param , # geneid parameter file 
                "o"        => \$blast_opt    , # tblastx options 
                "c"        => \$score_cutoff , # tblastx score cutoff 
                "s"        => \$shrink       , # shrink hsp's by
                "t"        => \$tbx          , # read tblastx from file
                "f"        => \$hsp          , # read HSP files in directory
                "k"        => \$ofn          , # intermediate filename
                "p"        => \$ps_output    , # postscript output 
                "v"        => \$verbose_flg  , # verbose    
                "h|help|?" => \$help_flg     , # print help
                );
    &prt_Help if $help_flg;
}; # sub Which_Options




#line 108 "./sgp2.nw"
&Which_Options();


&get_exec_time(time);

#line 92 "./sgp2.nw"
exit(1);

### EOF ###

#line 271 "./sgp2.nw"
__DATA__
=head1 NAME

    
$PROGRAM ($VERSION) - Improving Gene Prediction with Sinteny.

=head1 SYNOPSIS

    
    $PROGRAM [-hv] [-o \'options\'] [-g \'options\'] \
             [-P filename] [-p filename] [-k filename] \
             [-c value] [-s value] -1 seqfile_1 -2 seqfile_2

=head1 DESCRIPTION

=head1 OPTIONS

    

=over 4

=item B<-1> I<seqfile_1>

input file for first species.

=item B<-2> I<seqfile_2>

input file for second species.

=item B<-g>

geneid options

=item B<-o>

tblastx options

=item B<-c> I<value>

tblastx score cuttof

=item B<-s> I<value>

shrink hsp\'s by value

=item B<-t> I<filename>

read tblastx file

=item B<-f> I<prefix>

read hsp gff files with in directory prefix and extension .hsp-rs

=item B<-k> I<prefix>

keep intermediate files with prefix

=item B<-p> I<filename>

ps output in filename file 

=item B<-P> I<filename>

geneid parameter file

=item B<-v>

verbose mode

=item B<-h>

produces this message

=back

=head1 FILES

=head1 DIAGNOSTICS

=head1 REQUIRES

=head1 BUGS

    
Report any problem to: B<jabril@imim.es>

=head1 AUTHOR

    
Roderic Guigo  B<rguigo@imim.es>

Josep F. Abril B<jabril@imim.es>

$PROGRAM is under GNU-GPL (C) 2000
