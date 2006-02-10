#!/usr/local/bin/perl -w

use strict;
use CGI;

my $cgi = new CGI;
my $matrix = $cgi->param ('matrix');
my $threshold = $cgi->param ('threshold');

# print "Content-type: text/html\n\n";

# print "hello!<br><br>";
# print "being implemented...<br>";

print "Content-type: image/png\n\n";
# my $picture = qx/cat clustering.png/;

my $path_to_script = "/home/ug/gmaster/projects/moby/prod/scripts/workflows_implementations";
my $path_to_data   = "/usr/local/Install/apache2/cgi-bin/moby";

my $picture = qx/$path_to_script\/GenesClustering_FASTA.pl -x 3 -c $path_to_script\/config.txt -f $path_to_data\/Homo_sapiens.fa/;

print $picture;
