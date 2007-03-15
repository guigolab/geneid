#!/usr/local/bin/perl -w

use strict;
use EnvServices;

BEGIN {
  load_environment();
}

use CGI;

use HTML::Template;

my $APACHE_ROOT = $ENV{'APACHE_ROOT'};
my $FRAME  = "$APACHE_ROOT/htdocs/software/geneid/Plantilla.html";
my $FRAME2 = "$APACHE_ROOT/htdocs/software/geneid/Plantilla2.html";

my $_debug = 0;

my $cgi    = new CGI;
my @params = $cgi->param;

my $_template_page = "$APACHE_ROOT/cgi-bin/moby/template/workflow_result.tmpl";
my $_path_to_html_page = "$APACHE_ROOT/htdocs/webservices/workflows/results";
my $jobid;

if (defined($cgi->param ('jobid'))) {
    $jobid = $cgi->param ('jobid');
}
else {
  print_error ("Error, no job identifier specified!\n");
}

if ($_debug && defined $jobid) {
  print STDERR "jobid, $jobid\n";
}

my $outfile = $_path_to_html_page . "/GENE_CLUSTERING." . $jobid . ".html";

if ($_debug) {
  print STDERR "outifle: $outfile\n";
}

# Check if the HTML out file has been generated
# if yes, return this page
# if not, return the retrieval page

my $template;

if (-f $outfile && ! -z $outfile) {
  $template = HTML::Template->new (filename => $outfile);
}
else {
  # open the html template
  $template = HTML::Template->new(filename => $_template_page);
  # fill in some parameters
  $template->param(JOBID => $jobid);
}

##############################
#
# Generate an asynchronous retrieval HTML page
#
##############################

if ($_debug) {
    print STDERR "Return the HTML page...\n";
}
     
print "Content-type: text/html\n\n";
print $template->output;
          
##############################
#
# The End
#
##############################

sub print_error {
# el parametro es un mensaje de error
    my @mess = @_;
    
    # imprimiendo inicio de la plantilla    
    # open(OUTPRINT,$FRAME);
    # while (my $line = <OUTPRINT>) {
	# print "$line";
    # }
    # close(OUTPRINT); 
    
    # Instead
    print "<html><body>";
    
    print "<br><center><font size=6 color=red>Genes Clustering server: ERROR management</font></center><br>";
    print "</FONT>";
    print "</TD>";
    print "</TR>";
    print "</TABLE>";
    print "<P>";
    print @mess;
    print "<hr>";
    print "<P>";
    
    # imprimiendo final de la plantilla  
    open(OUTPRINT,$FRAME2);
    while (my $line = <OUTPRINT>) {
	print "$line";
    }  
    close(OUTPRINT);
} 

