#!/usr/local/bin/perl -w

use strict;
use CGI;

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
    # for testing benchmarking using NFS!
    # ($seq_fh, $seqfile) = tempfile("/usr/local/Install/apache2/htdocs/test/GENE_CLUSTERING_INPUT.XXXXXX", UNLINK => 0);
};

if ($_debug) {
    print STDERR "temporary input file, $seqfile\n";
}

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
	print_error("<b>ERROR> No DNA sequences were submitted.</b><br><br>Please, fill the textarea in or select a file for submitting promoter sequences");
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
    print_error("<b>ERROR> No DNA sequences were submitted.</b><br><br>Please, fill the textarea in or select a file for submitting promoter sequences");

    exit 1;
}

my $t2 = Benchmark->new ();
if ($_debug) {
    print STDERR "\nGene clustering Input Uploading : ", timestr (timediff ($t2, $t1)), "\n";
}

# Get The type of input, so we know if we have to deal with FASTA sequences or a list of genes

my $input_type = $cgi->param ('input_type');

if ($_debug) {
    print STDERR "input type, $input_type\n";
}

if (!defined $input_type || $input_type ne "FASTA" || $input_type ne "LIST") {
    print STDERR "Error, the input type has not been set up properly!\n";
    if (defined $input) {
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
    
    if ($nb_sequences > 30) {
	print "Content-type: text/html\n\n";
	print_error("<b>ERROR> Too many sequences have been submitted, there is a limit of 30!");
	exit 1;
    }
    
    my $nb_bases = qx/grep -v ">" $seqfile | awk '{l+=length($1)}END{print l}'/;
    
    if ($_debug) {
	print STDERR "number of input bases, $nb_bases\n";
    }
    
    if ($nb_bases > 80000) {
	print "Content-type: text/html\n\n";
	print_error("<b>ERROR> Too long sequences have been submitted, the overall limit is 80 000bp!");
	exit 1;
    }
}
else {
    
    $script_name = "GenesClustering_FASTA.pl";
    
    my $nb_genes = qx/wc -l $seqfile/;
    
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

# ...

if ($_debug) {
    print STDERR "executing the following command,\n";
    print STDERR "$_path_to_script\/$script_name -x 2 -c $_path_to_script\/workflow.config -d $matrix -t $threshold -a $alpha -l $lambda -u $mu -m $nj_method -n $cluster_number -i $iteration_number -g $gamma -r $non_colinear -f $seqfile -o $gene_clustering_output_dir\n";
}

my $failure;
if ($_debug) {
    $failure = qx/$_path_to_script\/$script_name -x 2 -c $_path_to_script\/workflow.config -d $matrix -t $threshold -a $alpha -l $lambda -u $mu -m $nj_method -n $cluster_number -i $iteration_number -g $gamma -r $non_colinear -f $seqfile -o $gene_clustering_output_dir/;
}
else {
    $failure = qx/$_path_to_script\/$script_name -x 2 -c $_path_to_script\/workflow.config -d $matrix -t $threshold -a $alpha -l $lambda -u $mu -m $nj_method -n $cluster_number -i $iteration_number -g $gamma -r $non_colinear -f $seqfile -o $gene_clustering_output_dir >& \/dev\/null/;
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
	
	print "<li><h3>cluster $cluster_index:</h3>";
	# print "</ul>\n";
	print join ("<br>", @genes);
	print "<br><br>\n";
	
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
		
		print "<font color=blue><b>Graphical representation of the TF-map alignment:</b><br><br></font>\n";
		print "<img src=\"$cluster_image_webpath\"><br>\n";
		
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
    # get rid of the input file, $seqfile
    unlink $seqfile;
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
    print "<P>";
    
    # imprimiendo final de la plantilla  
    open(OUTPRINT,$FRAME2);
    while (my $line = <OUTPRINT>) {
	print "$line";
    }  
    close(OUTPRINT);
} 

