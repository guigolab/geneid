#!/usr/local/bin/perl -w

use strict;
use CGI;

use File::Temp qw/tempfile/;

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

# print "Content-type: text/html\n\n";
# print "hello!<br><br>";
# print "being implemented...<br>";
# exit 0;

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
}

# Get the parameters

my $matrix         = $cgi->param ('matrix');
my $threshold      = $cgi->param ('threshold');
my $method         = $cgi->param ('method');
my $cluster_number = $cgi->param ('clusters');

if (! ($cluster_number =~ /\d+/)) {
    print STDERR "number of clusters is not numerical, $cluster_number!\n";
}

my $path_to_script = "/home/ug/gmaster/projects/moby/prod/scripts/workflows_implementations";

print STDERR "executing the gene clustering workflow...\n";

# my $picture = qx/$path_to_script\/GenesClustering_FASTA.pl -x 2 -c $path_to_script\/workflow.config -d $matrix -t $threshold -m $method -n $cluster_number -f $seqfile >& \/dev\/null/;
my $result = qx/$path_to_script\/GenesClustering_FASTA.pl -x 2 -c $path_to_script\/workflow.config -d $matrix -t $threshold -m $method -n $cluster_number -f $seqfile -o \/tmp\/output_clustering/;

print STDERR "execution done\n";

if ((-f "/tmp/output_clustering/clustering_tree.png") && !(-z "/tmp/output_clustering/clustering_tree.png")) {
    my $picture = qx/cat \/tmp\/output_clustering\/clustering_tree.png/;
    
    if ($_debug) {
	print STDERR "got a picture!\n";
    }
    
    if (!$_debug) {
	unlink $seqfile;
    }
    print "Content-type: image/png\n\n";
    print $picture;
}
else {
    
    # Get the clusters
    
    my @clusters = ();

    #...
    
    if (@clusters < 1) {
	
	print STDERR "no clusters found, genes clustering failed!!\n";
	
	if (!$_debug) {
	    unlink $seqfile;
	}
	print "Content-type: text/html\n\n";
	print_error("<b>ERROR> The execution of the genes clustering workflow has failed!");
    }
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

    exit(1);
} 

