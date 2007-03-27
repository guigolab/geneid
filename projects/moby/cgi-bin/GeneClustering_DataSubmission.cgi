#!/usr/local/bin/perl -w

# Gene Clustering data submission page
# it executes in asynchronous mode the workflow
# it uses ajax to add the data input HTML submission form

use strict;
use EnvServices;

BEGIN {
  load_environment();
}

use CGI;
use CGI::Ajax;

use HTML::Template;

my $_debug = 0;

my $APACHE_ROOT = $ENV{'APACHE_ROOT'};

if (!defined $APACHE_ROOT) {
  # MACOSX
  # Also i have added a symb limk htdocs -> Documents
  $APACHE_ROOT = "/Library/WebServer";  
}
else {
  if ($_debug) {
    print STDERR "APACHE_ROOT: $APACHE_ROOT\n";
  }
}

my $FRAME  = "$APACHE_ROOT/htdocs/software/geneid/Plantilla.html";
my $FRAME2 = "$APACHE_ROOT/htdocs/software/geneid/Plantilla2.html";

my $_template_top_page  = $APACHE_ROOT . "/cgi-bin/moby/template/submission_top.tmpl";
my $_template_down_page = $APACHE_ROOT . "/cgi-bin/moby/template/submission_down.tmpl";

my $cgi = new CGI;

if ($_debug) {
    my @params = $cgi->param;
    print STDERR "data submission param names, @params\n";
}

# CGI script URL to give to the template on the fly
my $url = $cgi->url(-relative => 1);

my $ajax = CGI::Ajax->new( check_datatype => \&check_datatype);
# $ajax->JSDEBUG(1);
print $ajax->build_html( $cgi, \&main );

# Without Ajax
# main ();

exit 0;

##############################
#
# The End
#
##############################

sub check_datatype {
  my ($input_type) = @_;
  my $data_input_html;

  if ($_debug && defined $input_type) {
    print STDERR "checking datatype, $input_type\n";
  }
  
  if (!defined $input_type) {
    if ($_debug) {
      print STDERR "the input type has not been defined\n";
    }
    my $_template_page = $APACHE_ROOT . "/cgi-bin/moby/template/datatype_submission.tmpl";
    
    if (! -f $_template_page) {
      print STDERR "Error, can't find template page, $_template_page\n";
    }
    
    my $_template = HTML::Template->new(filename => $_template_page);
    $data_input_html = $_template->output;
  }
  elsif ($input_type eq "FASTA") {
    my $_template_page = $APACHE_ROOT . "/cgi-bin/moby/template/fasta_submission.tmpl";
    my $_template = HTML::Template->new(filename => $_template_page);
    $data_input_html = $_template->output;
  }
  elsif ($input_type eq "LIST") {
    my $_template_page = $APACHE_ROOT . "/cgi-bin/moby/template/list_submission.tmpl";
    my $_template = HTML::Template->new(filename => $_template_page);
    $data_input_html = $_template->output;
  }
  else {
    print STDERR "don't know anything about this input type, $input_type\n";
    my $_template_page = $APACHE_ROOT . "/cgi-bin/moby/template/gene_clustering_entrypoint.tmpl";
    my $_template = HTML::Template->new(filename => $_template_page); 
    $data_input_html = $_template->output;
  }
  
  return "$data_input_html";
}

sub main {

my $input_type = $cgi->param ('datatype');

if ((! -f $_template_top_page) || (! -f $_template_down_page)) {
  print STDERR "can't find the template page(s), $_template_top_page, $_template_down_page\n";
  print "Content-type: text/html\n\n";
  print_error("<b>ERROR> Internal Error");
  exit 1;
}

# open the html template
my $first_template = HTML::Template->new(filename => $_template_top_page);
my $third_template = HTML::Template->new(filename => $_template_down_page);
     
my $data_input_html = check_datatype ($input_type);     

if ($_debug) {
  print STDERR "Check done\n";
} 
    
# Without Ajax
# print "Content-type: text/html\n\n";
# my $template = HTML::Template->new (filename => "$APACHE_ROOT/cgi-bin/moby/template/gene_clustering_entrypoint.tmpl");
# print $template->output;

# Ajax
my $html;
$html = $first_template->output;
$html .= $data_input_html;
$html .= $third_template->output;

return "$html";

} # End main
          
          
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

