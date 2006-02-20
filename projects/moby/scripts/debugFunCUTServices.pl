#!/usr/bin/perl -w
 
use strict;

use MOBY::Client::Central; 
use SOAP::Lite + 'trace'; 
use MOBY::Client::Service; 

use Time::HiRes qw(time);
 
# be prepare for command-line options/arguments
use Getopt::Std;

sub help {
return <<"END_HELP";
Description: Debug the PDG Services (FunCUT or ISS). 
Usage:

debugPDGServices.pl [-h] -x {Moby Central} -s {Service Name} [-r {Rounds}] [-f {Filters}] [-d {Database}] [-e {e-value}] [-m {Maximum Searches to ISS}] [-l {Cutlen}] [-c {Autoarc}] [-w {Framework}] [-a {FASTA file}] [-i {Input File}]
	-h help
	-x MOBY Central: Chirimoyo, Mobydev, Inab or BioMoby
		<1> or Chirimoyo
		<2> or Mobydev
		<3> or Inab
		<4> or BioMoby
	-s Service Name:
		<1> or runFunCUT
		<2> or runISS
		<3> or runISSComplete
		<4> or runNCut
		<5> or parserISS_Output_Into_Text_Xml
		<6> or parserISS_Output_Into_NCBI_Blast_Text
		<7> or parserISS_Output_Into_NCut_Input
		<8> or runOFunCUT
		<9> or fromFunCUTtoGFF
	-r Number of rounds (from 2 until 6). By default is 3 rounds.
	-f Filters: [XNU,SEG,COILS]. It is possible the combination of every one or none of them.
	-d Database: [Swiss-Prot,TrEMBL,SWALL]. You have to choose one of them. By default is Swiss-Prot.
	-e Rest of parameters: e-value
	-m Maximum Searches used by ISS method
	-l Cutlen: Minimum length of intermediate sequences (optional, defval 35)
	-c Autoarc: E-value given to an autoarc (optional, defval 1e-20)
	-w Framework: (0->local Blast+raw FASTA file, 1->local Blast+raw SW file, 2->local Blast+SQUID, 3->Paracel+SQUID, 4->MOBY) (optional, defval 3)
	
	-a File which contains multiple aminoacid sequence in FASTA format (by default it uses the "P82197" sequence). This input is usually used by runFunCUT, runISS, and runISSComplete services.
	-p File in XML that contains an bioMoby object (ISS_Output). This input is used by parser services.
	-i File in XML that contains all information necessary to execute OFunCUT service (by default it uses the "P82197" sequence).

Examples using some combinations:
	perl debugPDGServices.pl -x Inab -s runFunCUT -r 5 -f XNU,SEG -d Swiss-Prot -a aaSequence.faa
	perl debugPDGServices.pl -x 2 -s 2 
	perl debugPDGServices.pl -x Mobydev -s runFunCUT -f COILS -d TrEMBL
	perl debugPDGServices.pl -x Inab -s runFunCUT -i multimpleAA.seq
	
END_HELP

}


BEGIN {
	
	# Determines the options with values from program
	use vars qw/$opt_h $opt_x $opt_s $opt_r $opt_f $opt_d $opt_e $opt_m $opt_l $opt_c $opt_w $opt_a $opt_p $opt_i/;
	   
	# these are switches taking an argument (a value)
	my $switches = 'hxsrfdemlcwapi';
	   
	# Get the switches
	getopt($switches);
	
	# If the user does not write nothing, skip to help
	if (defined($opt_h) || !defined($opt_x) || !defined($opt_s)){
		print help;
		exit 0;
	}
	
}

# Global variables
my $URL;
my $URI;
my $PROXY;

#----------------------------------------------
#	ASSIGN THE MOBY URI AND MOBY SERVER
#----------------------------------------------
if (defined($opt_x)) {

	# Delete spaces
	$opt_x =~ s/\s//g;

	# Assign the MOBY Server and MOBY URI
	if (($opt_x == 1) || ($opt_x eq 'Chirimoyo')) {
	
		# export MOBY_URI
		# export MOBY_SERVER
		$URL = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://chirimoyo.ac.uma.es/cgi-bin/MOBY-Central.pl';
		$URI = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://chirimoyo.ac.uma.es/MOBY/Central';
		$PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';
	
	}elsif (($opt_x == 2) || ($opt_x eq 'Mobydev')) {
	
		# export MOBY_URI
		# export MOBY_SERVER
		$URL = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://moby-dev.inab.org/cgi-bin/MOBY-Central.pl';
		$URI = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://moby-dev.inab.org/MOBY/Central';
		$PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';

	}elsif (($opt_x == 3) || ($opt_x eq 'Inab')) {

		# export MOBY_URI
		# export MOBY_SERVER
		$URL = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'http://www.inab.org/cgi-bin/MOBY-Central.pl';
		$URI = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'http://www.inab.org/MOBY/Central';
		$PROXY = $ENV{MOBY_PROXY}?$ENV{MOBY_PROXY}:'No Proxy Server';
	
	}elsif (($opt_x == 4) || ($opt_x eq 'BioMoby')) {

		# export MOBY_URI
		# export MOBY_SERVER
		#$URL = $ENV{MOBY_SERVER}?$ENV{MOBY_SERVER}:'';
		#$URI = $ENV{MOBY_URI}?$ENV{MOBY_URI}:'';
		#$ENV{MOBY_URI}='http://mobycentral.icapture.ubc.ca/MOBY/Central';
		#$ENV{MOBY_SERVER}='http://mobycentral.icapture.ubc.ca/cgi-bin/MOBY05/mobycentral.pl';

	}else {
		print help;
		exit 0;
	}

}else {
	print help;
	exit 0;
}

#---------------------------------------
# GET THE SERVICE NAME
#---------------------------------------

# Declare variable which contains service name
my $serviceName;

# Check the option of service name.
if (($opt_s == 1) || ($opt_s eq 'runFunCUT')) {
	$serviceName = 'runFunCUT';
	
} elsif (($opt_s == 2) || ($opt_s eq 'runISS')) {
	$serviceName = 'runISS';
	
} elsif (($opt_s == 3) || ($opt_s eq 'runISSComplete')) {
	$serviceName = 'runISSComplete';
	
} elsif (($opt_s == 4) || ($opt_s eq 'runNCut')) {
	$serviceName = 'runNCut';
	
} elsif (($opt_s == 5) || ($opt_s eq 'parserISS_Output_Into_Text_Xml')) {
	$serviceName = 'parserISS_Output_Into_Text_Xml';
	
} elsif (($opt_s == 6) || ($opt_s eq 'parserISS_Output_Into_NCBI_Blast_Text')) {
	$serviceName = 'parserISS_Output_Into_NCBI_Blast_Text';
	
} elsif (($opt_s == 7) || ($opt_s eq 'parserISS_Output_Into_NCut_Input')) {
	$serviceName = 'parserISS_Output_Into_NCut_Input';
	
} elsif (($opt_s == 8) || ($opt_s eq 'runOFunCUT')) {
	$serviceName = 'runOFunCUT';
	
} elsif (($opt_s == 9) || ($opt_s eq 'fromFunCUTtoGFF')) {
	$serviceName = 'fromFunCUTtoGFF';
	
} else {

	print "The service name is not correct\n";
	exit 0;
}


#---------------------------------------
# GET THE XML INPUT LIST
#---------------------------------------

# Declare XML list 
my @XMLinputlist;

#-------------------------------------------------------------------------------------------------------
#		EXECUTE FunCUT USING ONE SEQUENCE
#-------------------------------------------------------------------------------------------------------
# Check if the user wants to test using one sequence.
if (defined($opt_s) && !defined($opt_a) && !defined($opt_i) && !defined($opt_p)) {

	# String that stores the AminoacidSequence object
	my $AAobject = "
	<moby:AminoAcidSequence namespace='UniProt' id='P82197'>
	<moby:Integer namespace='' id='' articleName='Length'>312</moby:Integer>
	<moby:String namespace='' id='' articleName='SequenceString'>MEEECRVLSIQSHVVRGYVGNRAATFPLQVLGFEVDAVNSVQFSNHTGYSHWKGQVLNSDELQELYDGLKLNHVNQYDYVLTGYTRDKSFLAMVVDIVQELKQQNPRLVYVCDPVMGDQRNGEGAMYVPDDLLPVYREKVVPVADIITPNQFEAELLTGRKIHSQEEALEVMDMLHSMGPDTVVITSSNLLSPRGSDYLMALGSQRTRAPDGSVVTQRIRMEMHKVDAVFVGTGDLFAAMLLAWTHKHPNNLKVACEKTVSAMHHVLQRTIKCAKAKSGEGVKPSPAQLELRMVQSKKDIESPEIVVQATVL</moby:String>
	</moby:AminoAcidSequence>";
		
	# Declare XML varibles 
	my @XML;
	
	# Add the AAsequence object
	push(@XML, '');
	push(@XML, $AAobject);
	
	# Add filters into XML list
	if (defined($opt_f)){

		# Delete spaces
		$opt_f =~ s/\s*//;
		
		# Obtain every filters.
		my @filters = split(/,/, $opt_f);
		
		# Add filters into XML list
		push(@XML, 'filters');
		push(@XML, "<Value>@filters</Value>");
	}
		
	# Add rounds into XML list
	if (defined($opt_r)) {
		
		# Add filters into XML list
		push(@XML, 'rounds');
		push(@XML, "<Value>$opt_r</Value>");
	}
	
	# Add database into XML list
	if (defined($opt_d)) {
		
		# Add filters into XML list
		push(@XML, 'database');
		push(@XML, "<Value>$opt_d</Value>");
		
	}
	
	if (defined($opt_e) || defined($opt_m) || defined($opt_l) || defined($opt_c) || defined($opt_w)) {
		
		push(@XML, 'restparameters');
		my $restparameters = "<Value>";
		
		# Add e-value into XML list
		if (defined($opt_e)) {$restparameters .= "evalue=$opt_e ";}
		
		# Add maxsearches into XML list
		if (defined($opt_m)) {$restparameters .= "maxsearches=$opt_m ";}
		
		# Add cutlen into XML list
		if (defined($opt_l)) {$restparameters .= "cutlen=$opt_l ";}
		
		# Add autoarc into XML list
		if (defined($opt_c)) {$restparameters .= "autoarc=$opt_c ";}
		
		# Add framework into XML list
		if (defined($opt_w)) {$restparameters .= "framework=$opt_w ";}
		
		$restparameters .= "</Value>";
		push(@XML, $restparameters);
		
	}
	
	# Add the input and secondary parameters
	push(@XMLinputlist, [@XML]);

	# Case 2: ERROR-The object is not son of the correct one.
	my $AAobject2 = "
	<moby:GenericSequence namespace='UniProt' id='P82197'>
	<moby:String namespace='' id='' articleName='SequenceString'>MEEECRVLSIQSHVVRGYVGNRAATFPLQVLGFEVDAVNSVQFSNHTGYSHWKGQVLNSDELQELYDGLKLNHVNQYDYVLTGYTRDKSFLAMVVDIVQELKQQNPRLVYVCDPVMGDQRNGEGAMYVPDDLLPVYREKVVPVADIITPNQFEAELLTGRKIHSQEEALEVMDMLHSMGPDTVVITSSNLLSPRGSDYLMALGSQRTRAPDGSVVTQRIRMEMHKVDAVFVGTGDLFAAMLLAWTHKHPNNLKVACEKTVSAMHHVLQRTIKCAKAKSGEGVKPSPAQLELRMVQSKKDIESPEIVVQATVL</moby:String>
	</moby:GenericSequence>";
	my @XML2;
	push(@XML2, '');
	push(@XML2, $AAobject2);
#	push(@XMLinputlist, [@XML2]);
	
	# Case 3: ERROR-The object is built badly
	my $AAobject3 = "
	<moby:AminoAcidSequence namespace='UniProt' id='P82197'>
	<moby:Integer namespace='' id='' articleName='Length'>312</moby:Integer>
	<moby:String namespace='' id='' articleName='SequenceSt'>MEEECRVLSIQSHVVRGYVGNRAATFPLQVLGFEVDAVNSVQFSNHTGYSHWKGQVLNSDELQELYDGLKLNHVNQYDYVLTGYTRDKSFLAMVVDIVQELKQQNPRLVYVCDPVMGDQRNGEGAMYVPDDLLPVYREKVVPVADIITPNQFEAELLTGRKIHSQEEALEVMDMLHSMGPDTVVITSSNLLSPRGSDYLMALGSQRTRAPDGSVVTQRIRMEMHKVDAVFVGTGDLFAAMLLAWTHKHPNNLKVACEKTVSAMHHVLQRTIKCAKAKSGEGVKPSPAQLELRMVQSKKDIESPEIVVQATVL</moby:String>
	</moby:AminoAcidSequence>";
	my @XML3;
	push(@XML3, '');
	push(@XML3, $AAobject3);
#	push(@XMLinputlist, [@XML3]);

	# Case 4: ERROR-The ID is not correct
	my $AAobject4 = "
	<moby:AminoAcidSequence namespace='UniProt' id='AAAAA'>
	<moby:Integer namespace='' id='' articleName='Length'>312</moby:Integer>
	<moby:String namespace='' id='' articleName='SequenceString'>MEEECRVLSIQSHVVRGYVGNRAATFPLQVLGFEVDAVNSVQFSNHTGYSHWKGQVLNSDELQELYDGLKLNHVNQYDYVLTGYTRDKSFLAMVVDIVQELKQQNPRLVYVCDPVMGDQRNGEGAMYVPDDLLPVYREKVVPVADIITPNQFEAELLTGRKIHSQEEALEVMDMLHSMGPDTVVITSSNLLSPRGSDYLMALGSQRTRAPDGSVVTQRIRMEMHKVDAVFVGTGDLFAAMLLAWTHKHPNNLKVACEKTVSAMHHVLQRTIKCAKAKSGEGVKPSPAQLELRMVQSKKDIESPEIVVQATVL</moby:String>
	</moby:AminoAcidSequence>";
	my @XML4;
	push(@XML4, '');
	push(@XML4, $AAobject4);
#	push(@XMLinputlist, [@XML4]);
	
	# Case 5: ERROR-The object is a collection
	my @XML5;
	push(@XML5, '');
	push(@XML5, [$AAobject, $AAobject]);
#	push(@XMLinputlist, [@XML5]);
	
	# Case 6: ERROR- Secondary parameter wrong
	my $AAobject6 = "
	<moby:AminoAcidSequence namespace='UniProt' id='P82197'>
	<moby:Integer namespace='' id='' articleName='Length'>312</moby:Integer>
	<moby:String namespace='' id='' articleName='SequenceString'>MEEECRVLSIQSHVVRGYVGNRAATFPLQVLGFEVDAVNSVQFSNHTGYSHWKGQVLNSDELQELYDGLKLNHVNQYDYVLTGYTRDKSFLAMVVDIVQELKQQNPRLVYVCDPVMGDQRNGEGAMYVPDDLLPVYREKVVPVADIITPNQFEAELLTGRKIHSQEEALEVMDMLHSMGPDTVVITSSNLLSPRGSDYLMALGSQRTRAPDGSVVTQRIRMEMHKVDAVFVGTGDLFAAMLLAWTHKHPNNLKVACEKTVSAMHHVLQRTIKCAKAKSGEGVKPSPAQLELRMVQSKKDIESPEIVVQATVL</moby:String>
	</moby:AminoAcidSequence>";
	my @XML6;
	push(@XML6, '');
	push(@XML6, $AAobject6);
	push(@XML6, 'filtersAAA');
	push(@XML6, "<Value>XNU</Value>");
#	push(@XMLinputlist, [@XML6]);
	
	# Case 7: ERROR- Secondary parameter wrong
	my @XML7;
	push(@XML7, '');
	push(@XML7, $AAobject6);
	push(@XML7, 'filters');
	push(@XML7, "<ValueAA>XNU</ValueAA>");
#	push(@XMLinputlist, [@XML7]);

#----------------------------------------------------------------------------------------------------------
#		MULTIPLE SEQUENCE IN "FASTA" FORMAT USING THE SAME NUMBER OF ROUNDS, FILTERS AND DATABASE
#----------------------------------------------------------------------------------------------------------
# Check if the user wants to test using FASTA sequences.
}elsif (defined($opt_s) && !defined($opt_p) && !defined($opt_i) && defined($opt_a)) {

	# String that stores the AminoacidSequence sequence
	my $AASequence;
	
	# String that stores the AminoacidSequence object	
	my $AAobject;
	
	# Store the id of aminoacid sequence
	my $id;
	
	# Flag that determines when a XML list is created
	my $startXMLInput = 0;

	# Array that contains the filters
	my @filters = undef;
	
	
	#open file which contains the aminoacid sequence.
	open(FILE, $opt_a);
	
	# Declare file line
	my $fileLine;
	
	# Read the file
	while (<FILE>) {
	
		$fileLine = $_;
		
		# Delete the end of line
		chomp($fileLine);
		
		# It is the comment line
		if ($fileLine =~/^>(.*)/) {
		
			# Check if it is the moment to create a new XML input list
			if ($startXMLInput == 1) {
			
				# Declare XML varibles 
				my @XML;
				
				# Create the BioMOBY object aside of sequence
				$AAobject = "<moby:AminoAcidSequence namespace='UniProt' id='".$id."'><moby:Integer namespace='' id='' articleName='Length'>".length($AASequence)."</moby:Integer><moby:String namespace='' id='' articleName='SequenceString'>".$AASequence."</moby:String></moby:AminoAcidSequence>";
				
				# Add the AAsequence object
				push(@XML, '');
				push(@XML, $AAobject);
				
				# Add filters into XML list
				if (defined($opt_f)){
			
					# Delete spaces
					$opt_f =~ s/\s*//;
					
					# Obtain every filters.
					my @filters = split(/,/, $opt_f);
					
					# Add filters into XML list
					push(@XML, 'filters');
					push(@XML, "<Value>@filters</Value>");
				}
				
				# Add rounds into XML list
				if (defined($opt_r)) {
					
					# Add filters into XML list
					push(@XML, 'rounds');
					push(@XML, "<Value>$opt_r</Value>");
				}
				
				# Add database into XML list
				if (defined($opt_d)) {
					
					# Add filters into XML list
					push(@XML, 'database');
					push(@XML, "<Value>$opt_d</Value>");
					
				}
			
				if (defined($opt_e) || defined($opt_m) || defined($opt_l) || defined($opt_c) || defined($opt_w)) {
					
					push(@XML, 'restparameters');
					my $restparameters = "<Value>";
					
					# Add e-value into XML list
					if (defined($opt_e)) {$restparameters .= "evalue=$opt_e ";}
					
					# Add maxsearches into XML list
					if (defined($opt_m)) {$restparameters .= "maxsearches=$opt_m ";}
					
					# Add cutlen into XML list
					if (defined($opt_l)) {$restparameters .= "cutlen=$opt_l ";}
					
					# Add autoarc into XML list
					if (defined($opt_c)) {$restparameters .= "autoarc=$opt_c ";}
					
					# Add framework into XML list
					if (defined($opt_w)) {$restparameters .= "framework=$opt_w ";}
					
					$restparameters .= "</Value>";
					push(@XML, $restparameters);
					
				}

				# Add the input and secondary parameters
				push(@XMLinputlist, [@XML]);
			
				# Undef all global variables
				$AASequence = undef;
				@filters = undef;
				$id = undef;
				
				# Change the flag: Now it is possible to finish the MOBY object
				$startXMLInput = 0;
			
			}
			
			# Split the comment line to get the id
			my @commentLine = split(/ /, $1);
			
			# Get the Id from aminoacid sequence
			$id = $commentLine[0];
			
		
		} else {
		
			# Save the Aminoacid sequence
			$AASequence .= $fileLine;
			
			# Start the flag which determines if the next time pass over a comment line, has to create a new XML input list
			$startXMLInput = 1;
		
		}
			
	}
	
	# The last object of the file has not been stored into XML input list
	if (eof){
	
		# Declare XML varibles 
		my @XML;
		
		# Create the BioMOBY object aside of sequence
		$AAobject = "<moby:AminoAcidSequence namespace='UniProt' id='".$id."'><moby:Integer namespace='' id='' articleName='Length'>".length($AASequence)."</moby:Integer><moby:String namespace='' id='' articleName='SequenceString'>".$AASequence."</moby:String></moby:AminoAcidSequence>";
				
		# Add the AAsequence object
		push(@XML, '');
		push(@XML, $AAobject);
		
		# Add filters into XML list
		if (defined($opt_f)){
			
			# Delete spaces
			$opt_f =~ s/\s*//;
			
			# Obtain every filters.
			my @filters = split(/,/, $opt_f);
		
			# Add filters into XML list
			push(@XML, 'filters');
			push(@XML, "<Value>@filters</Value>");
		}
		
		# Add rounds into XML list
		if (defined($opt_r)) {
			
			# Add filters into XML list
			push(@XML, 'rounds');
			push(@XML, "<Value>$opt_r</Value>");
		}
				
		# Add database into XML list
		if (defined($opt_d)) {
			
			# Add filters into XML list
			push(@XML, 'database');
			push(@XML, "<Value>$opt_d</Value>");
			
		}
				
		if (defined($opt_e) || defined($opt_m) || defined($opt_l) || defined($opt_c) || defined($opt_w)) {
			
			push(@XML, 'restparameters');
			my $restparameters = "<Value>";
			
			# Add e-value into XML list
			if (defined($opt_e)) {$restparameters .= "evalue=$opt_e ";}
			
			# Add maxsearches into XML list
			if (defined($opt_m)) {$restparameters .= "maxsearches=$opt_m ";}
			
			# Add cutlen into XML list
			if (defined($opt_l)) {$restparameters .= "cutlen=$opt_l ";}
			
			# Add autoarc into XML list
			if (defined($opt_c)) {$restparameters .= "autoarc=$opt_c ";}
			
			# Add framework into XML list
			if (defined($opt_w)) {$restparameters .= "framework=$opt_w ";}
			
			$restparameters .= "</Value>";
			push(@XML, $restparameters);
		
		}

		# Add the input and secondary parameters
		push(@XMLinputlist, [@XML]);
		
		# Close the file
		close(FILE);
	
	}
	
#-------------------------------------------------------------------------------------------------------
#		ISS_OUTPUT INPUT WILL BE USING BY PARSING SERVICES
#-------------------------------------------------------------------------------------------------------
# Check if the user wants to test using one sequence.
}elsif (defined($opt_s) && !defined($opt_a) && !defined($opt_i) && defined($opt_p)) {

	my $outputFILE;
	
	# Declare XML varibles 
	my @XML;

	# output file content packaging
	$/=undef;
	local(*INPUTFILE);
	unless (open(INPUTFILE,$opt_p)) {
		
		# Error: we can not find the ISS file.
		print "ERROR: Can't open: $!\n";
		
	} else {
	
		# Store result from service into variable.
		$outputFILE = <INPUTFILE>;
		
		# Add the AAsequence object
		push(@XML, '');
		push(@XML, $outputFILE);
		
		# Add the input
		#push(@XMLinputlist, ['', $outputFILE]);

		# Close files
		close(INPUTFILE);
	}
	$/='\n';
	
	# Add the input and secondary parameters
	push(@XMLinputlist, [@XML]);


#-----------------------------------------------------------------------------------------------------
#		XML FILE WILL BE INSERTING DIRECTLY USED ONLY BY OFUNCUT
#-----------------------------------------------------------------------------------------------------
# Check if the user has inserted multiple inputs by means of a file.
}elsif (defined($opt_s) && !defined($opt_a) && !defined($opt_p) && defined($opt_i)) {
	
	# Declare XML varibles 
	my @XML;
	
	#open file which contains the aminoacid sequence.
	open(INPUTFILE, $opt_i) || die "Can't open: $!\n";
	
	# Store result from service into variable.
	my @outputFILE = <INPUTFILE>;
	
	my ($okAA) = 0;
	my ($aaSeq);
	my ($okBlast) = 0;
	my ($firstBlast);
	my ($okncut) = 0;
	my ($ncutMat);
	
	foreach my $fileLine (@outputFILE) {
	
		if (($fileLine =~ /(.*)moby:CommentedAASequence>(.*)/g) && $okAA) {
			$okAA = 0;
			$aaSeq .= $fileLine;
		} elsif (($fileLine =~ /(.*)moby:NCBI_BLAST_Text>(.*)/g) && $okBlast) {
			$okBlast = 0;
			$firstBlast .= $fileLine;
		} elsif (($fileLine =~ /(.*)(moby:NCut_Matrix>)(.*)/g) && $okncut) {
			$okncut = 0;
			$ncutMat .= $fileLine;
		} elsif (($fileLine =~ /<moby:CommentedAASequence(.*)/) || $okAA) {
			$okAA = 1;
			$aaSeq .= $fileLine;
		} elsif (($fileLine =~ /<moby:NCBI_BLAST_Text(.*)/g) || $okBlast) {
			$okBlast = 1;
			$firstBlast .= $fileLine;
		} elsif (($fileLine =~ /<moby:NCut_Matrix(.*)/g) || $okncut) {
			$okncut = 1;
			$ncutMat .= $fileLine;
		}
		
	}
	
	# Add the AAsequence object
	push(@XML, 'AASequence');
	push(@XML, $aaSeq);
	
	# Add the firstBlast object
	push(@XML, 'firstBlast');
	push(@XML, $firstBlast);
	
	# Add the firstBlast object
	push(@XML, 'NcutMatrix');
	push(@XML, $ncutMat);
	
	# Add the input and secondary parameters
	push(@XMLinputlist, [@XML]);

	# Close the file
	close(INPUTFILE);
	
# Error the user has not configured the options correctly
}else {
	# Show the help
	print help;
	exit 0;

}

#----------------------------------------------
# FIND, RETRIEVE AND EXECUTE THE SERVICE
#----------------------------------------------

# Start time.
my $startTime = time;

# Create MOBY Central to MOBY Server and MOBY URI
my $MOBYCentral = MOBY::Client::Central->new(
			Registries => {mobycentral => {URL => $URL,URI => $URI}});

# Find the service
my ($s, $r) = $MOBYCentral-> findService(
				authURI => "pdg.cnb.uam.es", 
				serviceName => $serviceName);
$s = shift @{$s}; 
 
# Retrieve the service
my $wsdl = $MOBYCentral-> retrieveService($s); 

# Create new instance of MOBY Service
my $MOBYService = MOBY::Client::Service-> new( service => $wsdl); 

# Start time.
my $startTime2 = time;

# Invoke the service twice (in a single message).
my $result = $MOBYService-> execute(XMLinputlist => \@XMLinputlist);

# End time.
my $endTime = time;

# Print difference time.
my $totalTime = $endTime - $startTime;
my $totalTime2 = $endTime - $startTime2;

print "Time used by service execution is $totalTime2 sec.\n";
print "Total time used by process is $totalTime sec.\n\n";

print $result;


1;
