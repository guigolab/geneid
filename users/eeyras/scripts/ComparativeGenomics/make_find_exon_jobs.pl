#!/usr/local/bin/perl -w

use strict;


my $jobs_file = "ALL_JOBS";

open (JOBS_FILE, ">$jobs_file");

my $job_counter = 1;
my $counter = 0;
my @job_lines = ();
while (<>){
    chomp;
    my @e = split;
    
    next unless ($e[0]=~/ENS/ && $e[1]=~/ENS/);
    my $human_id = $e[0];
    my $mouse_id = $e[1];
    
    my $job = "/home/ug/eeyras/bin/scripts/ComparativeGenomics/find_new_coding_exons_from_gene_pair_using_e2t.pl $human_id $mouse_id\n";
    push ( @job_lines, $job );
    $counter ++;

    if ( $counter%30 == 0 ){
	
	#print "counter: $counter\n";
	
	my $job_name = "coding_job_using_e2t".$job_counter;
	open (JOB,">$job_name") or die ("cannot open $job_name");    
        foreach my $line (@job_lines){
	    print JOB $line;
	}
	close(JOB);
	@job_lines = ();
	$job_counter++;
	my $stderr = $job_name.".err";
	my $stdout = $job_name.".out";
	print JOBS_FILE "qsub -q main -e $stderr -o $stdout $job_name\n";
    }
    
}

if (@job_lines){
    my $job_name = "coding_job_using_e2t".$job_counter;
    open (JOB,">$job_name") or die ("cannot open $job_name");    
    print JOB @job_lines;
    close(JOB);
    @job_lines = ();
    $job_counter++;
    my $stderr = $job_name.".err";
    my $stdout = $job_name.".out";
    print JOBS_FILE "qsub -q main -e $stderr -o $stdout $job_name\n";
}


close(JOBS_FILE);
system("chmod u+x $jobs_file");
