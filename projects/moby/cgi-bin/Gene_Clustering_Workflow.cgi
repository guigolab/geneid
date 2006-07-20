#!/usr/local/bin/perl -w

use strict;
use CGI;

use File::Temp qw/tempfile/;
use File::Temp qw/tempdir/;

my $_debug = 0;

my $APACHE_ROOT = $ENV{'APACHE_ROOT'};
my $FRAME  = "$APACHE_ROOT/htdocs/software/geneid/Plantilla.html";
my $FRAME2 = "$APACHE_ROOT/htdocs/software/geneid/Plantilla2.html";

my $cgi = new CGI;

# Get the sequences input and store it into a temprary file

my ($seq_fh, $seqfile);
eval {
    ($seq_fh, $seqfile) = tempfile("/tmp/GENE_CLUSTERING_INPUT.XXXXXX", UNLINK => 0);
};

my @params = $cgi->param;
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
	    if ($_debug) {
		print STDERR "copying data...\n";
		# print STDERR "line, $line\n";
	    }
	    
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

# Get the parameters

my $matrix           = $cgi->param ('matrix');
my $threshold        = $cgi->param ('threshold');

my $alpha            = $cgi->param ('alpha');
my $lambda           = $cgi->param ('lambda');
my $mu               = $cgi->param ('mu');

my $nj_method        = $cgi->param ('method');

my $iteration_number = $cgi->param ('iterations');
my $cluster_number   = $cgi->param ('clusters');

if (! ($cluster_number =~ /\d+/)) {
    print STDERR "number of clusters is not numerical, $cluster_number!\n";
}

if (! ($iteration_number =~ /\d+/)) {
    print STDERR "number of k-means iterations is not numerical, $iteration_number!\n";
}

if ($_debug) {
    print STDERR "alpha, lambda, mu,  $alpha, $lambda, $mu\n";
}

my $path_to_script = "/home/ug/gmaster/projects/moby/prod/scripts/workflows_implementations";

if ($_debug) {
    print STDERR "executing the gene clustering workflow...\n";
}

my $gene_clustering_output_dir = tempdir( "/tmp/GENE_CLUSTERING_OUTPUT.XXXXXX" );

if ($_debug) {
    print STDERR "executing the following command,\n";
    print STDERR "$path_to_script\/GenesClustering_FASTA.pl -x 2 -c $path_to_script\/workflow.config -d $matrix -t $threshold -a $alpha -l $lambda -u $mu -m $nj_method -n $cluster_number -i $iteration_number -f $seqfile -o $gene_clustering_output_dir\n";
}

my $result = qx/$path_to_script\/GenesClustering_FASTA.pl -x 2 -c $path_to_script\/workflow.config -d $matrix -t $threshold -a $alpha -l $lambda -u $mu -m $nj_method -n $cluster_number -i $iteration_number -f $seqfile -o $gene_clustering_output_dir/;

if ($_debug) {
    print STDERR "execution done\n";
}

if ((-f "$gene_clustering_output_dir/clustering_tree.png") && !(-z "$gene_clustering_output_dir/clustering_tree.png")) {
    my $picture = qx/cat $gene_clustering_output_dir\/clustering_tree.png/;
    
    if ($_debug) {
	print STDERR "got a picture!\n";
    }
    
    print "Content-type: image/png\n\n";
    print $picture;
}
else {
    
    # Get the clusters
    
    my @clusters = ();
    
    opendir (CLUSTERDIR, "$gene_clustering_output_dir/K-means_clusters");
    my @cluster_files = grep { $_ ne '.' and $_ ne '..' } readdir CLUSTERDIR;
    
    closedir CLUSTERDIR;
    
    if (@cluster_files > 0) {
	
	# Make a n archive
	my $archive_path = "/usr/local/Install/apache2/htdocs/webservices/workflows/results";
	my $output_dir_name = $gene_clustering_output_dir;
	$output_dir_name =~ s/\/tmp\///;
	my $archive_filename = $output_dir_name . ".zip";
	
	print STDERR "output_dir_name, $output_dir_name\n";
	
	qx/cd \/tmp; zip -r $archive_path\/$archive_filename $output_dir_name/;
	
	print "Content-type: text/html\n\n";
	
	print "<html><head><title>Gene clustering results</title></head>\n<body>";
	
	my $archive_URL = "http://genome.imim.es/webservices/workflows/results/" . $archive_filename;
	
	print "There is an <a href=\"$archive_URL\">archive</a> available to download all the results\n";
	
	print "<ul>\n";
	my $index = 1;
	foreach my $cluster_file (@cluster_files) {
	    print STDERR "parsing file, $cluster_file\n";
	    my $genes = qx/cat $gene_clustering_output_dir\/K-means_clusters\/$cluster_file/;
	    $genes =~ s/\n/<br>/g;
	    
	    print "<li><h3>cluster $index:</h3><br>";
	    print "$genes<br>";
	    
	    $index++;
	}
	
	print "</ul>\n";
	print "</body></html>\n";
    }
    else {
	print STDERR "no clusters found, genes clustering failed!!\n";
	
	print "Content-type: text/html\n\n";
	print_error("<b>ERROR> The execution of the genes clustering workflow has failed!");
    }
}

if (!$_debug) {
    # get rid of output directory, $gene_clustering_output_dir, and input file, $seqfile
    unlink $seqfile;
    # ...
}

#################################################

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

