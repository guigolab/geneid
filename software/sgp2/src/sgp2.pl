######################################################################
#                               sgp2                                 #
######################################################################
#
#     Sinteny based Gene Prediction tool.
#
#     Copyright (C) 2003 -          Genis PARRA FARRE
#                          Josep Francesc ABRIL FERRANDO
#                                 Roderic GUIGO SERRA
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
######################################################################


use strict;
use Getopt::Long;

##########################################################################
##                                                                      ##
##                     VARIABLES  AND DEFINITIONS                       ##
##                                                                      ##
##########################################################################

# MAIN VARIABLES 
my $PROGRAM = "SGP2";
my $VERSION = "1.0";
 
my $SGP2      = defined($ENV{'SGP2'}) ? $ENV{'SGP2'} : '.';
my $BIN       = defined($ENV{'SGP2'}) ? "$SGP2/bin" : '.';
my $SGP2param = defined($ENV{'SGP2'}) ? "$SGP2/param" : "../param";
my $TMP       = defined($ENV{'SGP2TMP'}) ? $ENV{'SGP2TMP'} : '/tmp';
my $TMPROOT   = "sgp2_$$";
my $SGPTMP    = "$TMP/$TMPROOT";

# PROGRAMS DEFINITION 
my $PARSEBLAST = "parseblast";
my $GENEID     = "geneid-sgp";
my $HSP2SR     = "blast2gff";

# VARIABLES USED TO MODIFY SR SCORES
my $S_CUTOFF = 26;          # hsp wiht lower score discarded
my $SCF = 0;                # substract to tblastx scores S_CUTOFF - SCF
my $DSC = $S_CUTOFF - $SCF; # this value is substracted from the hsp score
my $SHSP  = 0;              # shrink hsp by $SHSP
my $WTBX  = 0.19;           # weigth of tblastx score

my $GeneidParam   = "$SGP2param/human3iso.param.sgp";
my $GeneidOptions = "";
my $verbose_str = "";

##########################################################################
##                                                                      ##
##                          READING OPTIONS                             ##
##                                                                      ##
##########################################################################

# print help and exit if any option is given
&prt_Help 
	unless scalar(@ARGV); 

# getOptions Variables
my ( $Seq1, $Seq2, $geneid_opt, $geneid_param, 
     $blast_opt, $score_cutoff, $shrink, $tbx,
     $hsp, $ofn, $ps_output, $savefiles_flg, $verbose_flg, 
     $quiet_flg, $temp_flg, $help_flg
     ) = ( undef, undef, undef, undef, undef, undef,
           undef, undef, undef, undef, undef, 0, 0, 0, 0, 0 );

# sequence identifiers
my ($Id1, $Id2);

# file flags
my ($blast_flg, $hsps_flg, $geneid_flg, $plots_flg) = (1, 1, 1, 1);

# processes
my $status;

# reading options
&Which_Options();


##########################################################################
##                                                                      ##
##                          INITIAL CHECKINGS                           ##
##                                                                      ##
##########################################################################

# Checking if all programs exists and are executable
(-e "$BIN/$PARSEBLAST" && -x _) || do {
   &go_to_die("$PARSEBLAST is not found or executable in $BIN")
};
(-e "$BIN/$HSP2SR" && -x _) || do {
   &go_to_die("$HSP2SR is not found or executable in $BIN")
};
(-e "$BIN/$GENEID" && -x _) || do {
   &go_to_die("$GENEID is not found or executable in $BIN")
};


##########################################################################
##                                                                      ##
##                          COMPUTING DATA                              ##
##                                                                      ##
##########################################################################

## Running parseblast
&prt_Header("RUNNING $PARSEBLAST");

&run_parseblast("$TMP/$TMPROOT.query.hsp.gff","")
	if (defined($Seq1));

&run_parseblast("$TMP/$TMPROOT.subjct.hsp.gff","--subject")
	if (defined($Seq2));

## Running  blast2gff
&prt_Header("RUNNING $HSP2SR");

&run_hsp2sr("$TMP/$TMPROOT.query.hsp.gff","$TMP/$TMPROOT.query.sr.gff")
	if (defined($Seq1));

&run_hsp2sr("$TMP/$TMPROOT.subjct.hsp.gff","$TMP/$TMPROOT.subjct.sr.gff")
	if (defined($Seq2));

## Modifying the score of the SR hits 
&prt_Header("FILTERING SR");

&filter_sr("$TMP/$TMPROOT.query.sr.gff","$TMP/$TMPROOT.query.mod.sr.gff")
	if (defined($Seq1));

&filter_sr("$TMP/$TMPROOT.subjct.sr.gff","$TMP/$TMPROOT.subjct.mod.sr.gff")
	if (defined($Seq2));

## Running geneid
&prt_Header("RUNNING GENE PREDICTION");

&run_geneid("$TMP/$TMPROOT.query.mod.sr.gff","$Seq1")
	if (defined($Seq1));

&run_geneid("$TMP/$TMPROOT.subjct.mod.sr.gff","$Seq2")
	if (defined($Seq2));


##########################################################################
##                                                                      ##
##                              ENDING                                  ##
##                                                                      ##
##########################################################################

## Removing temporary files
&clean_tmp() unless $temp_flg;
exit(0);

##########################################################################
##                                                                      ##
##                              THE END                                 ##
##                                                                      ##
##########################################################################
##########################################################################


##########################################################################
##                                                                      ##
##                            SUBROUTINES                               ##
##                                                                      ##
##########################################################################

# Checking fasta format sequence files
sub check_fasta_format() {
    my $file = $_[0];
    my ($n, $c) = (undef, 0);
    open(TMP,"< $file");
    while (<TMP>) {
        next unless /^>/;
        />(\S+)\b/ && do {
            $n = $1;
            $c++;
            next unless $c>1;
        };
        &go_to_die("FATAL ERROR !!! Multiple locus names found.\n  File \'$file\' must contain only one sequence definition.");
    }
    &go_to_die("FATAL ERROR !!! There is no '>' line, locus name not found.\n  Please, verify your fasta file \'$file\'") unless defined($n);
    return $n;
}

# Deleting temporary files on TMP
sub clean_tmp() {
    # Obtaining the list of temporary files
    opendir(DIR,"$TMP");
    my @files = map { "$TMP/$_" } grep { /^$TMPROOT/ } readdir(DIR);
    closedir(DIR);
    # Unlinking the temporary files if they exist
    ($#files >= 0) and 
        ( unlink(@files)
          and &prt_Header("********* Temporary files were deleted *********")
          or &prt_Header("********* Can't unlink @files : $! *********")
        ) or 
          &prt_Header("********* There are no temporary files in $TMP *********");
}

# Checking input sequence files
sub exists_file() {
    my @files = @_;
    my ($n, $r) = (' ', 0);
    foreach $n (@files) {
        $r++ if (-e "$n");
    };
    return $r;
}

# Get a fixed length string from a given string and filling char/s.
sub fill_mid() { 
   my $l = length($_[0]); 
   my $k = int(($_[1] - $l)/2); 
   return ($_[2] x $k).$_[0].($_[2] x ($_[1] - ($l+$k)));
}

# writing die messages to STDERR and clean_tmp before exit.
sub go_to_die() { 
   (print STDERR "\n   @_ \n\n") && &clean_tmp(); 
   exit(1) 
}

# section headers to STDERR
sub prt_Header() { 
   print STDERR ("*" x 80)."\n** ".&fill_mid("@_",74," ")." **\n".("*" x 80)."\n" unless $quiet_flg;
}
# print the command to be executed
sub prt_Command() { 
	print STDERR ("\n\n** ".&fill_mid("@_",74," ")." **\n\n");
}

# the following routines run parseblast,hsp2sr, filter the sr and run geneid 
sub run_parseblast { 
	my ($outputfile, $option) = @_;

	my $parse = "$BIN/$PARSEBLAST --bit-score --gff $option $verbose_str $tbx | gawk '{OFS=\"\\t\"; 
      print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,0}' > $outputfile;";
    &prt_Command("$parse") if $verbose_flg;     

	$status = system($parse);
	if ($status != 0) {
		&go_to_die("FATAL ERROR !!! $? $!\n
          Unsuccessful $PARSEBLAST command :\n  $parse")
	};   
}

sub run_hsp2sr() {
	my ($inputfile,$outputfile) = @_;

	my $projection = "$BIN/$HSP2SR $verbose_str -g $inputfile | sort +3n > $outputfile;";
    &prt_Command($projection) if $verbose_flg;     

	$status = system $projection;
	if ($status != 0) {
		&go_to_die("FATAL ERROR !!! $? $!\n
          Unsuccessful $HSP2SR  command :\n  $projection")
	};
}   

sub filter_sr(){
	my ($inputfile,$outputfile) = @_;

	if (!open(Srfile,"< $inputfile")) {
		&go_to_die ("SGP2.pl: $inputfile: no such temporary file");
	}

	open(Srmodfile,"> $outputfile");

	while (<Srfile>) {
		next if /^\#/;
		next if /^[ \t]*$/;
		chomp;	
		my @gff_record=split(/\s+/,$_);
		$gff_record[5]=(($gff_record[5] - $DSC) * $WTBX);
		print Srmodfile join("\t",@gff_record)."\n";
	};

	close(Srmodfile);
	close(Srfile);
}

sub run_geneid(){
	my ($srfile,$fastafile) = @_;

	my $prediction = "$BIN/$GENEID $verbose_str  $GeneidOptions -P $GeneidParam  -S $srfile $fastafile;";
    &prt_Command("$prediction") if $verbose_flg;     

	$status = system $prediction;
	if ($status != 0) {
		&go_to_die("FATAL ERROR !!! $? $!\n
          Unsuccessful $GENEID command :\n  $prediction");
	};   
}

# get a fixed length string from a given string and filling char/s.
sub Which_Options() {
    GetOptions( 
             "1|query=s"   => \$Seq1         , # seqfile_1
             "2|sbjct=s"   => \$Seq2         , # seqfile_2
             "g=s"         => \$geneid_opt   , # geneid options      
             "P=s"         => \$geneid_param , # geneid parameter file 
             "c=f"         => \$score_cutoff , # tblastx score cutoff
             "s=f"         => \$shrink       , # shrink hsp's by
             "t|tblastx=s" => \$tbx          , # read tblastx from file
             "f=s"         => \$hsp          , # read HSP files in directory
             "v|verbose"   => \$verbose_flg  , # verbose    
             "q|quiet"     => \$quiet_flg    , # quiet mode
             "tmp"         => \$temp_flg     , # keep temporary files
             "h|help|?"    => \$help_flg     , # print help
                ) || &prt_Help();
    &prt_Help if $help_flg;
    &prt_Header("Processing Command-Line Options");

   # do seq files exists
   &go_to_die("FATAL ERROR!!! --query or --sbjct and --tblastx options are mandatory.") 
	   unless ((defined($Seq1) || defined($Seq2)) && defined($tbx));
   
   my $file_number = &exists_file(defined($Seq1) ? $Seq1 : $Seq2, $tbx);
   &go_to_die("FATAL ERROR!!! Input files do not exist.")
	   unless $file_number == 2; 
     
   # check if files are provided in fasta format and get locus names
   $Id1 = &check_fasta_format($Seq1) if defined($Seq1);
   $Id2 = &check_fasta_format($Seq2) if defined($Seq2);

   #&go_to_die("FATAL ERROR!!! 
   #    Identifier1($Id1) have the same name as Indentifier2($Id2).\n  
   #    Sequences \'$Seq1\' and \'$Seq2\' 
   #    must have different locus names.\n")
   #   if $Id1 eq $Id2;
   # setting other variables

   &go_to_die("FATAL ERROR!!! --verbose and --quiet are incompatible options.")
	   if ($verbose_flg && $quiet_flg);

   $GeneidOptions = "-"."$geneid_opt"   if defined($geneid_opt);
   $GeneidParam   = "$geneid_param" if defined($geneid_param);
   $verbose_str = " -v " if $verbose_flg;
} 

# prints help 
sub prt_Help() {
    open(HELP, "| cat") ;
    print HELP <<"+++EndOfHelp+++";
                                                                     $PROGRAM


PROGRAM:
                        $PROGRAM - $VERSION

               Syntenic Gene Prediction tool 

USAGE:
        
    $PROGRAM [-vhq][-g option] <[-1|-2] fasta_sequence> <-t tblastx_aligment>

DESCRIPTION:

    $PROGRAM is a program to predict genes by comparing anonymous
    genomic sequences from two different species. It combines tblastx,
    a sequence similarity search program, with geneid, an "ab initio"
    gene prediction program.

REQUIRES:

    $PROGRAM needs standard Perl distribution installed in your system.

ENVIRONMENT VARIABLES:

    You can specify the path where $PROGRAM can find the default files
    with the shell variable \"$PROGRAM\". 

    You also can specify the path for the tempotary files with the
    shell variable \"SGP2TMP\". Default value is /tmp.

	Setting those vars in Bourne-shell and C-shell:

     o Using a Bourne-Shell (e.g. bash):
           export SGP2="path"
           export SGP2TMP="path"

     o Using a C-Shell:
           setenv SGP2 "path"
           setenv SGP2TMP "path"

 
COMMAND-LINE OPTIONS:

	 Available options and a short description are listed here:

     -1, --query     fasta file of the query sequence in the tblastx.
     -2, --sbjct     fasta file of the subject sequence in the tblastx.
     -t, --tblastx   tblastx results file.
     -g, --geneid <option>   
                    geneid output options:
                     G - Use GFF format to print predictions. 
                     X - Use extended-format to print gene predictions 
                     D - Output genomic sequence of exons in predicted genes

     -v, --verbose   Verbose mode, a report is sent to standard error 
                      (default is set to showing only WARNINGS).
     -q, --quiet     Quiet mode, do not show any message/warning
                     to standard error (only ERRORS are reported).
     -h, --help      Shows this help.

     --tmp           Keep temporary files.

BUGS:    
  
    Report any problem to 'SGP2\@imim.es'.

AUTHOR:  

    $PROGRAM has been developed by Genis Parra, Josep Francesc Abril
    and Roderic Guigo.



GNU-GPL (C)                         July 2003                         $PROGRAM
+++EndOfHelp+++
    close(HELP);
    exit(1);
} # prt_help

