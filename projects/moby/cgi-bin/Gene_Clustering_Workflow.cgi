#!/usr/local/bin/perl -w

use strict;
use CGI;

use File::Temp qw/tempfile/;
use File::Temp qw/tempdir/;

my $_debug = 0;

my $_path_to_script = "/home/ug/gmaster/projects/moby/prod/scripts/workflows_implementations";

my $APACHE_ROOT = $ENV{'APACHE_ROOT'};
my $FRAME  = "$APACHE_ROOT/htdocs/software/geneid/Plantilla.html";
my $FRAME2 = "$APACHE_ROOT/htdocs/software/geneid/Plantilla2.html";

my $cgi    = new CGI;
my @params = $cgi->param;

##############################
#
# INPUT
#
##############################

# Get the sequences input and store it into a temprary file

my ($seq_fh, $seqfile);
eval {
    ($seq_fh, $seqfile) = tempfile("/tmp/GENE_CLUSTERING_INPUT.XXXXXX", UNLINK => 0);
};

if ($_debug) {
    print STDERR "clustering param names, @params\n";
}

if (defined($cgi->param ('sequences')) && ($cgi->param ('sequences') =~ />/)) {

    if ($_debug) {
	print STDERR "sequences...\n";
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
	my $upload_filehandle = $cgi->upload('upfile');
	
	# Copy the uploaded data into a temporary file
	while ( my $line = <$upload_filehandle> ) {
	    print $seq_fh $line;
	}
    }
    else {
	close $seq_fh;
	unlink $seqfile;

	print "Content-type: text/html\n\n";
	print_error("<b>ERROR> A DNA sequence has not been submitted.</b><br><br>Please, fill the textarea in or select a file for submitting a DNA sequence");
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
    print_error("<b>ERROR> A DNA sequence has not been submitted.</b><br><br>Please, fill the textarea in or select a file for submitting a DNA sequence");

    exit 1;
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

my $nj_method        = $cgi->param ('method');

my $iteration_number = $cgi->param ('iterations');
my $cluster_number   = $cgi->param ('clusters');

# Parameters Validation

if (! ($cluster_number =~ /\d+/)) {
    print STDERR "number of clusters is not numerical, $cluster_number!\n";
}

if (! ($iteration_number =~ /\d+/)) {
    print STDERR "number of k-means iterations is not numerical, $iteration_number!\n";
}

if ($_debug) {
    print STDERR "alpha, lambda, mu,  $alpha, $lambda, $mu\n";
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

if ($_debug) {
    print STDERR "executing the following command,\n";
    print STDERR "$_path_to_script\/GenesClustering_FASTA.pl -x 2 -c $_path_to_script\/workflow.config -d $matrix -t $threshold -a $alpha -l $lambda -u $mu -m $nj_method -n $cluster_number -i $iteration_number -f $seqfile -o $gene_clustering_output_dir\n";
}

my $result;
if ($_debug) {
    $result = qx/$_path_to_script\/GenesClustering_FASTA.pl -x 2 -c $_path_to_script\/workflow.config -d $matrix -t $threshold -a $alpha -l $lambda -u $mu -m $nj_method -n $cluster_number -i $iteration_number -f $seqfile -o $gene_clustering_output_dir/;
}
else {
    $result = qx/$_path_to_script\/GenesClustering_FASTA.pl -x 2 -c $_path_to_script\/workflow.config -d $matrix -t $threshold -a $alpha -l $lambda -u $mu -m $nj_method -n $cluster_number -i $iteration_number -f $seqfile -o $gene_clustering_output_dir >& \/dev\/null/;
}

if ($_debug) {
    print STDERR "execution done\n\n";
}

if (length ($result) > 1) {
    print STDERR "workflow execution has failed!\n";
    print STDERR "execution result, $result\n";

    print "Content-type: text/html\n\n";
    print_error("<b>ERROR> The execution of the genes clustering workflow has failed!");
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
$output_dir_name =~ s/\/tmp\///;
my $archive_filename = $output_dir_name . ".zip";
my $archive_URL = "http://genome.imim.es/webservices/workflows/results/" . $archive_filename;

$result = qx/cd \/tmp; zip -r $archive_path\/$archive_filename $output_dir_name/;

# Display in HTML the results
# First the archive link

print "Content-type: text/html\n\n";
print "<html><head><title>Gene clustering results</title></head>\n<body>";

print "There is an <a href=\"$archive_URL\">archive</a> available to download all the results\n";

if ($_debug) {
    print STDERR "Archive done!\n";
    print STDERR "zip result status, $result\n";
}

if ($_debug) {
    print STDERR "process the cluster results\n";
}

# Now, get the cluster results

print "<ul>\n";

my $cluster_index = 1;
while ($cluster_index <= $cluster_number) {
    
    # The results for current cluster
    
    my $cluster_directory_name = $cluster_index . ".cluster_results";
    
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
	
	print "<li><h3>cluster $cluster_index:</h3>";
	# print "</ul>\n";
	print join ("<br>", @genes);
	print "<br><br>\n";
	
	if ($nb_genes > 1) {
	    
	    # JPEG
	    
	    my $cluster_image_filepath = "/webservices/workflows/results/$gene_clustering_output_dirname/$cluster_directory_name/" . $cluster_index . ".TFBSs_maps.jpg";
	    
	    if ($_debug) {
		print STDERR "gff maps picture file path for cluster $cluster_index, $cluster_image_filepath\n";
	    }
	    
	    if (! -f $cluster_image_filepath) {
		print STDERR "can't find the gff maps picture for cluster $cluster_index - file path is $cluster_image_filepath\n";
	    }
	    
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
	    
	    print "<font color=blue><b>Graphical representation of the TF-map alignment:</b><br><br></font>\n";
	    print "<img src=\"$cluster_image_filepath\"><br>\n";
	    
	    print "<TABLE border=0 cellpadding=0 cellspacing=0 width=100%>\n";
	    print "<TR>\n";
	    print "<TD class='section'>\n";
	    print "<FONT size=5 class='K'>\n";
	    
	    print "<code>Multiple meta-alignment</code> predictions for this cluster are:</FONT></TD></TR></TABLE><P><pre><tt>\n";
	    
	    # Multiple meta-alignment results
	    print "$mmeta_results\n";
	    
	    print "<hr>\n";
	    
	    # MatScan results
	    print "$matscan_results\n";
	    
	    print "</pre></tt>\n";
	}
	
	if ($_debug) {
	    print STDERR "cluster processing done\n";
	}
	
    } # End processing current cluster

    $cluster_index++;
} # End processing all cluster results

print "</ul>\n";
print "</body></html>\n";

if ($_debug) {
    print STDERR "processing of the cluster results done\n";
}

##############################
#
# Cleaning out
#
##############################

if (!$_debug) {
    # get rid of output directory, $gene_clustering_output_dir, and input file, $seqfile
    unlink $seqfile;
    # ...
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
    print "<p>";
    print "<b>List of incompatibilities and suggestions:</b>";
    print "<ul>";
    print "<li>A DNA sequence in FASTA format must be always provided either by cut&paste or a file";
    print "<li>To obtain a graphical representation of the predictions, please set the format field to <tt>gff</tt>";
    print "<li>Submitted sequences must be lower than 100 Kbps when a graphical representation must be developed";
    print "<li> Assembling mode requires entering some experimental evidences in GFF format";
    print "</ul>";
    print "<P>";
    
    # imprimiendo final de la plantilla  
    open(OUTPRINT,$FRAME2);
    while (my $line = <OUTPRINT>) {
	print "$line";
    }  
    close(OUTPRINT);
} 

