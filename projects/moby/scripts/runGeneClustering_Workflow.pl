#!/usr/local/bin/perl -w

# This script is called by the Gene clustering Workflow CGI script
# The cgi script delegates to it the execution of the workflow and the generation of the output HTML page.

# To do..............................
# Add input options handling code

use strict;

use File::Temp qw/tempfile/;
use File::Temp qw/tempdir/;

# Benchmark Module
use Benchmark;

my $t1 = Benchmark->new ();

my $_debug = 0;

my $_path_to_script = "/home/ug/gmaster/projects/moby/prod/scripts/workflows_implementations";

my $APACHE_ROOT = $ENV{'APACHE_ROOT'};
my $FRAME  = "$APACHE_ROOT/htdocs/software/geneid/Plantilla.html";
my $FRAME2 = "$APACHE_ROOT/htdocs/software/geneid/Plantilla2.html";

##############################
#
# OUTPUT HTML FILE
#
##############################

# input argument
my $output_html_file;
# at the last moment, copy the result from the temp file to the visible file, $output_html_file
# Because asynchronous test is based on whether the file exists or not !
my ($out_fh, $temp_output_html_file) = tempfile("/tmp/GENE_CLUSTERING_HTML_OUTPUT.XXXXXX", UNLINK => 0);

##############################
#
# INPUT
#
##############################

# Get the sequences input and store it into a temporary file

my ($seq_fh, $seqfile);

if ($_debug) {
    print STDERR "temporary input file, $seqfile\n";
}

if (-z $seqfile) {
    
    if ($_debug) {
	print STDERR "Empty input file, $seqfile...\n";
    }
    unlink $seqfile;
    
    print "Content-type: text/html\n\n";
    print_error("<b>ERROR> No data were submitted.</b><br><br>Please, fill the textarea in or select a file");

    close $out_fh;
    exit 1;
}

my $t2 = Benchmark->new ();
if ($_debug) {
    print STDERR "\nGene clustering Input Uploading : ", timestr (timediff ($t2, $t1)), "\n";
}

# Get The type of input, so we know if we have to deal with FASTA sequences or a list of genes

my $input_type;

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
    close $out_fh;
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
    
}
else {
    
    $script_name = "GenesClustering_LIST.pl";
    
}
   
##############################
#
# CGI PARAMETERS
#
##############################

my $matrix;
my $threshold;

my $alpha;
my $lambda;
my $mu;

my $nj_method;

my $iteration_number;
my $cluster_number;

my $gamma;
my $non_colinear;

# Parameters Validation

if (! ($cluster_number =~ /\d+/)) {
    print STDERR "number of clusters is not numerical, $cluster_number!\n";
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

my $gene_clustering_output_dir = tempdir( "/usr/local/Install/apache2/htdocs/webservices/workflows/results/GENE_CLUSTERING_OUTPUT.XXXXXX" );
my $gene_clustering_output_dirname = $gene_clustering_output_dir;
$gene_clustering_output_dirname =~ s/\/usr\/local\/Install\/apache2\/htdocs\/webservices\/workflows\/results\///;

# Make the arguments line

my $args = "-x 2 -c $_path_to_script\/workflow.config -d $matrix -t $threshold -a $alpha -l $lambda -u $mu -m $nj_method -n $cluster_number -i $iteration_number -g $gamma -r $non_colinear -f $seqfile -o $gene_clustering_output_dir";

if ($input_type eq "LIST") {
    $args .= " -s $species -v $upstream_length -w $downstream_length";
}

if ($_debug) {
    print STDERR "executing the following command,\n";
    print STDERR "$_path_to_script\/$script_name $args\n";
}

my $failure;
if ($_debug) {
    $failure = qx/$_path_to_script\/$script_name $args/;
}
else {
    $failure = qx/$_path_to_script\/$script_name $args >& \/dev\/null/;
}

if ($_debug) {
    print STDERR "workflow execution result, $failure\n";
}

if ($_debug) {
    print STDERR "execution done\n\n";
}

if ($failure) {
    print STDERR "workflow execution has failed!\n";

    print "Content-type: text/html\n\n";
    print_error("<b>ERROR> The execution of the genes clustering workflow has failed!");
    close $out_fh;
    exit 1;
}

##############################
#
# Workflow results processing
#
##############################

if ($_debug) {
    print STDERR "Make an archive...\n";
}

# Make an archive

my $archive_path = "/usr/local/Install/apache2/htdocs/webservices/workflows/results";
my $output_dir_name = $gene_clustering_output_dir;
$output_dir_name =~ s/\/usr\/local\/Install\/apache2\/htdocs\/webservices\/workflows\/results\///;
my $archive_filename = $output_dir_name . ".zip";
my $archive_URL = "http://genome.imim.es/webservices/workflows/results/" . $archive_filename;

my $result = qx/cd $archive_path; zip -r $archive_path\/$archive_filename $output_dir_name/;

# Display in HTML the results
# First the archive link

print $out_fh "Content-type: text/html\n\n";
print $out_fh "<html><head><title>Gene clustering results</title></head>\n<body>";

print $out_fh "There is an <a href=\"$archive_URL\">archive</a> available to download all the results\n";

if ($_debug) {
    print STDERR "Archive done!\n";
    print STDERR "zip result status, $result\n";
}

if ($_debug) {
    print STDERR "process the cluster results\n";
}

# Now, get the cluster results

print $out_fh "<ul>\n";

my $cluster_index = 1;
while ($cluster_index <= $cluster_number) {
    
    # The results for current cluster
    
    my $cluster_directory_name = $cluster_index . ".cluster_results";
    
    if ($_debug) {
	print STDERR "processing data from directory, $cluster_directory_name...\n";
    }
    
    if (-d "$gene_clustering_output_dir/$cluster_directory_name") {
	
	if ($_debug) {
	    print STDERR "processing cluster directory, $cluster_directory_name...\n";
	}
	
	# The gene members
	
	my $genes = qx/cat $gene_clustering_output_dir\/$cluster_directory_name\/$cluster_index.ids.lst/;
	my @genes = grep {$_ ne ""} split ("\n", $genes);
	
	if ($_debug) {
	    print STDERR "genes in cluster $cluster_index, @genes\n";
	}
	
	my $nb_genes = @genes;
	
	print $out_fh "<li><h3>cluster $cluster_index:</h3>";
	# print $out_fh "</ul>\n";
	print $out_fh join ("<br>", @genes);
	print $out_fh "<br><br>\n";
	
	if ($nb_genes > 1) {
	    
	    # JPEG
	    
	    my $cluster_image_filepath = "/usr/local/Install/apache2/htdocs/webservices/workflows/results/$gene_clustering_output_dirname/$cluster_directory_name/" . $cluster_index . ".TFBSs_maps.jpg";
	    my $cluster_image_webpath = "/webservices/workflows/results/$gene_clustering_output_dirname/$cluster_directory_name/" . $cluster_index . ".TFBSs_maps.jpg";
	    
	    if ($_debug) {
		print STDERR "gff maps picture file path for cluster $cluster_index, $cluster_image_filepath\n";
	    }
	    
	    if (! -f $cluster_image_filepath) {
		print STDERR "can't find the gff maps picture for cluster $cluster_index - file path is $cluster_image_filepath\n";
	    }
	    else {
		
		# MatScan
		
		my $matscan_results = "";
		
		# ...
		
		# MMeta
		
		my $mmeta_results = "";
		my $mmeta_filename = "$gene_clustering_output_dir/$cluster_directory_name/" . $cluster_index . ".MultiMeta.txt";
		
		if ($_debug) {
		    print STDERR "mmeta filename for cluster $cluster_index, $mmeta_filename\n";
		}
		
		if (-f $mmeta_filename) {
		    $mmeta_results = qx/cat $mmeta_filename/;
		}
		else {
		    print STDERR "can't find mmeta results for cluster $cluster_index - filename is $mmeta_filename\n";
		}
		
		# HTML
		
		print $out_fh "<font color=blue><b>Graphical representation of the TF-map alignment:</b><br><br></font>\n";
		print $out_fh "<img src=\"$cluster_image_webpath\"><br>\n";
		
		print $out_fh "<TABLE border=0 cellpadding=0 cellspacing=0 width=100%>\n";
		print $out_fh "<TR>\n";
		print $out_fh "<TD class='section'>\n";
		print $out_fh "<FONT size=5 class='K'>\n";
		
		print $out_fh "<code>Multiple meta-alignment</code> predictions for this cluster are:</FONT></TD></TR></TABLE><P><pre><tt>\n";
		
		# Multiple meta-alignment results
		print $out_fh "$mmeta_results\n";
		
		print $out_fh "<hr>\n";
		
		# MatScan results
		print $out_fh "$matscan_results\n";
		
		print $out_fh "</pre></tt>\n";
	    }
	}
	
	if ($_debug) {
	    print STDERR "cluster processing done\n";
	}
	
    } # End processing current cluster

    $cluster_index++;
} # End processing all cluster results

print $out_fh "</ul>\n";
print $out_fh "</body></html>\n";

close $out_fh;

qx/cp $temp_output_html_file $output_html_file/;

if ($_debug) {
    print STDERR "processing of the cluster results done\n";
}

##############################
#
# Cleaning out
#
##############################

if (!$_debug) {
    # get rid of the input file, $seqfile, and output file
    unlink $seqfile;
    unlink $temp_output_html_file;
}

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
	# print $out_fh "$line";
    # }
    # close(OUTPRINT); 
    
    # Instead
    print $out_fh "<html><body>";
    
    print $out_fh "<br><center><font size=6 color=red>Genes Clustering server: ERROR management</font></center><br>";
    print $out_fh "</FONT>";
    print $out_fh "</TD>";
    print $out_fh "</TR>";
    print $out_fh "</TABLE>";
    print $out_fh "<P>";
    print $out_fh @mess;
    print $out_fh "<hr>";
    print $out_fh "<P>";
    
    # imprimiendo final de la plantilla  
    open(OUTPRINT,$FRAME2);
    while (my $line = <OUTPRINT>) {
	print $out_fh "$line";
    }  
    close(OUTPRINT);
    
} 

