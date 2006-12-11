#!/usr/local/bin/perl -w

use strict;
use EnvServices;

BEGIN {
  load_environment();
}

use CGI;

use File::Temp qw/tempfile/;
use File::Temp qw/tempdir/;

# Benchmark Module
use Benchmark;

use HTML::Template;

my $t1 = Benchmark->new ();

my $_debug = 0;

my $_path_to_script = "/home/ug/gmaster/projects/moby/prod/scripts";
my $_script_name = "runGeneClustering_Workflow.pl";

my $APACHE_ROOT = $ENV{'APACHE_ROOT'};
my $FRAME  = "$APACHE_ROOT/htdocs/software/geneid/Plantilla.html";
my $FRAME2 = "$APACHE_ROOT/htdocs/software/geneid/Plantilla2.html";

my $cgi    = new CGI;
my @params = $cgi->param;

##############################
#
# OUTPUT
#
##############################

my ($out_fh, $outfile);
eval {
    ($out_fh, $outfile) = tempfile("/usr/local/Install/apache2/htdocs/webservices/workflows/results/GENE_CLUSTERING.XXXXXX", UNLINK => 0);
};
if ($@) {
  print STDERR "Error, can't create a temporary file\n";
    warn $@;
}
if (defined $out_fh) {
  close $out_fh;
}
unlink $outfile;
# Add the HTML extension
$outfile .= ".html";

if ($_debug) {
  print STDERR "output html file, $outfile\n";
  print STDERR "parsing the job id...\n";
}

$outfile =~ /GENE_CLUSTERING\.([^\.]+)\.html/;
my $jobid = $1;

print STDERR "Job ID: $jobid\n";

##############################
#
# INPUT
#
##############################

# Get the sequences input and store it into a temprary file

my ($seq_fh, $seqfile);
eval {
    ($seq_fh, $seqfile) = tempfile("/tmp/GENE_CLUSTERING_INPUT.XXXXXX", UNLINK => 0);
    # for testing benchmarking using NFS!
    # ($seq_fh, $seqfile) = tempfile("/usr/local/Install/apache2/htdocs/test/GENE_CLUSTERING_INPUT.XXXXXX", UNLINK => 0);
};

if ($_debug) {
    print STDERR "temporary input file, $seqfile\n";
}

if ($_debug) {
    print STDERR "clustering param names, @params\n";
}

if (defined($cgi->param ('sequences')) && (($cgi->param ('sequences') =~ />/) || (length ($cgi->param ('sequences')) > 0))) {

    if ($_debug) {
	print STDERR "parsing pasted sequences...\n";
    }
    
    my $sequences = $cgi->param ('sequences');
    $sequences =~ s/[\r]//g;
    
    print $seq_fh $sequences;
}
else {

    if ($_debug) {
	print STDERR "check for upfile param...\n";
    }
    
    if (defined ($cgi->param ('upfile'))) {	
	# Copy the uploaded data into a temporary file
	
	if ($_debug) {
	    print STDERR "uploading...\n";
	}
	
	# my $upload_filehandle = $cgi->upload('upfile');
	# while ( my $line = <$upload_filehandle> ) {
	    # print $seq_fh $line;
	# }
	
	# Seems faster that way!
	
	my $upfilename = $cgi->param ('upfile');
	binmode($seq_fh);
	while (my $bytesread = read($upfilename, my $buffer, 1024)) {
	    print $seq_fh $buffer || die "cannot write to file, $seqfile ($!)\n";
	}
    }
    else {
	close $seq_fh;
	unlink $seqfile;

	print "Content-type: text/html\n\n";
	print_error("<b>ERROR> No data were submitted.</b><br><br>Please, fill the textarea in or select a file");
	exit 1;
    }
}
close $seq_fh;

if (-z $seqfile) {
    
    if ($_debug) {
	print STDERR "Empty input file, $seqfile...\n";
    }
    unlink $seqfile;
    
    print "Content-type: text/html\n\n";
    print_error("<b>ERROR> No data were submitted.</b><br><br>Please, fill the textarea in or select a file");

    exit 1;
}

my $t2 = Benchmark->new ();
if ($_debug) {
    # print STDERR "\nGene clustering Input Uploading : ", timestr (timediff ($t2, $t1)), "\n";
}

# Get The type of input, so we know if we have to deal with FASTA sequences or a list of genes

my $input_type = $cgi->param ('input_type');

if ($_debug) {
    print STDERR "input type, $input_type\n";
}

if ((!defined $input_type) || (($input_type ne "FASTA") && ($input_type ne "LIST"))) {
    print STDERR "Error, the input type has not been set up properly!\n";
    if (defined $input_type) {
	print STDERR "input type, $input_type\n";
    }
    print "Content-type: text/html\n\n";
    print_error("<b>ERROR> Internal Error");
    exit 1;
}

# The script name that will be called, depending of the input type
my $script_name;

# Extraction parameters
my $species;
my $upstream_length;
my $downstream_length;

if ($input_type eq "FASTA") {
    
    $script_name = "GenesClustering_FASTA.pl";
    
    my $nb_sequences = qx/grep -c ">" $seqfile/;
    
    if ($_debug) {
	print STDERR "number of input sequences, $nb_sequences\n";
    }
    
    if ($nb_sequences > 60) {
	print "Content-type: text/html\n\n";
	print_error("<b>ERROR> Too many sequences have been submitted, there is a limit of 60!");
	exit 1;
    }
    
    my $nb_bases = qx/grep -v ">" $seqfile | awk '{l+=length($1)}END{print l}'/;
    
    if ($_debug) {
	print STDERR "number of input bases, $nb_bases\n";
    }
    
    if ($nb_bases > 80000) {
	print "Content-type: text/html\n\n";
	print_error("<b>ERROR> Too long sequences have been submitted, the overall limit is 80,000 bp!");
	exit 1;
    }
}
else {
    
    $script_name = "GenesClustering_LIST.pl";
    
    my $nb_genes = qx/grep -c "" $seqfile/;
    
    if ($nb_genes > 20) {
	print "Content-type: text/html\n\n";
	print_error("<b>ERROR> Too many genes have been submitted, there is a limit of 20!");
	exit 1;
    }
    
    # Parse the promoter sequence extraction parameters

    $species           = $cgi->param ('species');
    $upstream_length   = $cgi->param ('upstream');
    $downstream_length = $cgi->param ('downstream');
    
    my $nb_bases = $nb_genes * ($upstream_length + $downstream_length);
    
    if ($nb_bases > 80000) {
	print "Content-type: text/html\n\n";
	print_error("<b>ERROR> Too long sequences have been required, the overall limit is 80 000bp!");
	exit 1;
    }
}
   
##############################
#
# CGI PARAMETERS
#
##############################

my $matrix           = $cgi->param ('matrix');
my $threshold        = $cgi->param ('threshold');

my $alpha            = $cgi->param ('alpha');
my $lambda           = $cgi->param ('lambda');
my $mu               = $cgi->param ('mu');

my $nj_method        = $cgi->param ('method') || "nearest";

my $iteration_number = $cgi->param ('iterations');
my $cluster_number   = $cgi->param ('clusters');

my $gamma            = $cgi->param ('gamma');
my $non_colinear     = $cgi->param ('noncol');

# Parameters Validation

if (!defined $cluster_number || ! ($cluster_number =~ /\d+/)) {
    print STDERR "number of clusters not defined or is not numerical, $cluster_number!\n";
    print "Content-type: text/html\n\n";
    print_error("<b>ERROR> Error parsing the parameters!");
    exit 1;
}

if (! ($iteration_number =~ /\d+/)) {
    print STDERR "number of k-means iterations is not numerical, $iteration_number!\n";
}

if ($_debug) {
    print STDERR "NJ method, $nj_method\n";
}

if ($_debug) {
    print STDERR "alpha, lambda, mu, gamma, noncol, $alpha, $lambda, $mu, $gamma, $non_colinear\n";
}

##############################
#
# Workflow Execution
#
##############################

if ($_debug) {
    print STDERR "executing the gene clustering workflow...\n";
}

# Make the arguments line

my $args = "-e $input_type -o $outfile -d $matrix -t $threshold -a $alpha -l $lambda -u $mu -m $nj_method -n $cluster_number -i $iteration_number -g $gamma -r $non_colinear -f $seqfile";

if ($input_type eq "LIST") {
    $args .= " -s $species -v $upstream_length -w $downstream_length";
}

if ($_debug) {
    print STDERR "executing the following command: $_path_to_script\/$_script_name $args\n";
}

my $failure = 0;
my @commands = ("$_path_to_script\/$_script_name $args &");
system (@commands);

# if ($_debug && defined $failure) {
    # print STDERR "workflow submission result, $failure\n";
# }

if ($_debug) {
    print STDERR "submission done\n\n";
}

if ($failure) {
    print STDERR "workflow submission has failed!\n";

    print "Content-type: text/html\n\n";
    print_error("<b>ERROR> The submission of the genes clustering workflow has failed!");
    exit 1;
}

##############################
#
# Generate an asynchronous retrieval HTML page
#
##############################

if ($_debug) {
    print STDERR "Make a HTML page...\n";
}

my $_template_page = "/usr/local/Install/apache2/cgi-bin/moby/template/workflow_result.tmpl";

# open the html template
my $template = HTML::Template->new(filename => $_template_page);
# fill in some parameters
$template->param(JOBID => $jobid);
     
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

