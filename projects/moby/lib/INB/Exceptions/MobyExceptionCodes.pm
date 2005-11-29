package INB::Exceptions::MobyExceptionCodes;

use strict;


###############################################################################
# PROTOTYPES
###############################################################################
sub getExceptionCodeDescription($);

	
##############################################################################
# NAME: getExceptionCodeDescription
#
# DESCRIPTION: An exception code is received by giving input and then, 
#		the method retrieves itself exception code and exception description.
#		Also, if exception message input is defined, then it will be added into 
#		exception description output.
#		The exception error come from INB Exception Codes.
#		See. 0501INB-ExceptionReporting-v1.7 documentation
#
# INPUTS: 	- Exception Code
#		- Dynamic Exception Description
#
# OUTPUTS:	- INB Exception hash: Code, Description
#
# MODIFIED DATE: 02-Sep-2005
#
# AUTHOR: Jose Manuel Rodriguez Carrasco -jmrc@cnb.uam.es- (INB-CNB)
###############################################################################
sub getExceptionCodeDescription($) {
	
	my ($exceptionCode) = @_;
	my ($INB_Exception) = undef;
	
switch: {

# ERROR CODES DEALING WITH ANALYSIS DATA
	if ($exceptionCode == 200) { # UNKNOWN NAME: [200] "Setting input data under a non-existing name, or asking for a result using an unknown name"
		
		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Setting input data under a non-existing name, or asking for a result using an unknown name.",
		};
		
	} elsif ($exceptionCode == 201) { # INPUTS INVALID: [201] "Input data are invalid; not match with its definitions, or with its dependency condition"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Input data are invalid; not match with its definitions, or with its dependency condition.",
		};
		
	} elsif ($exceptionCode == 202) { # INPUT NOT ACCEPTED: [202] "Input data are accepted"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Input data are accepted.",
		};

	} elsif ($exceptionCode == 221) { # INPUT REQUIRED PARAMETER: [221] "Service require parameter X"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Service require parameter.",
		};
	
	} elsif ($exceptionCode == 222) { # INPUT INCORRECT PARAMETER: [222] "Incorrect parameter X"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Incorrect parameter.",
		};

	} elsif ($exceptionCode == 223) { # INPUT INCORRECT SIMPLE: [223] "Incorrect input in simple article"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Incorrect input in simple article.",
		};

	} elsif ($exceptionCode == 224) { # INPUT INCORRECT SIMPLENB: [224] "Service requires two or more simple articles"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Service requires two or more simple articles.",
		};

	} elsif ($exceptionCode == 225) { # INPUT INCORRECT COLLECTION: [225] "Incorrect input in collection article"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Incorrect input in collection article.",
		};

	} elsif ($exceptionCode == 226) { # INPUT EMPTY OBJECT: [226] "Empty input object"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Empty input object.",
		};

	} elsif ($exceptionCode == 231) { # INPUT EMPTY MOBYCONTENT: [231] "Empty MOBYContent"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Empty MOBYContent.",
		};

	} elsif ($exceptionCode == 232) { # INPUT EMPTY MOBYCONTENT: [232] "QueryID does not exists"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "QueryID does not exists.",
		};

	} elsif ($exceptionCode == 233) { # INPUT EMPTY MOBYDATA: [233] "Empty MOBYData"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Empty MOBYData.",
		};

# EXCEPTION CODES DEALING WITH ANALYSIS EXECUTION
	} elsif ($exceptionCode == 300) { # NOT RUNNABLE: [300] "The same job has already been executed, or the data that had been set previously do not exist or are not accessible anymore"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "The same job has already been executed, or the data that had been set previously do not exist or are not accessible anymore.",
		};

	} elsif ($exceptionCode == 301) { # NOT RUNNING: [301] "A job has not yet been started"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "The job has not yet been started.",
		};

	} elsif ($exceptionCode == 302) { # NOT TERMINATED: [302] "A job is not interruptible for some reason"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "The job is not interruptible for some reason.",
		};

# EXCEPTION CODES DEALING WITH ANALYSIS EXECUTION
	} elsif ($exceptionCode == 400) { # NO METADATA AVAILABLE: [400] "There are no metadata available"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "There are no metadata available.",
		};

# EXCEPTION CODES DEALING WITH NOTIFICATION
	} elsif ($exceptionCode == 500) { # PROTOCOLS UNACCEPTED: [500] "Server does not agree on using any of the proposed notification protocols"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Server does not agree on using any of the proposed notification protocols.",
		};

# GENERAL EXCEPTION CODES
	} elsif ($exceptionCode == 600) { # INTERNAL PROCESSING ERROR: [600] "A generic catch-all for errors not specifically mentioned elsewhere in this list"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "A generic error during internal processing.",
		};

	} elsif ($exceptionCode == 601) { # COMMUNICATION FAILURE: [601] "A generic network failure"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "A generic network failure.",
		};

	} elsif ($exceptionCode == 602) { # UNKNOWN STATE: [602] "Used when a network call expects to find an existing state but failed"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Unknown State.",
		};

	} elsif ($exceptionCode == 603) { # NOT IMPLEMENTED: [603] "Not implemented method in question"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Not implemented method in question.",
		};

# NUEVO------------------------------------------------------
# NUEVO------------------------------------------------------
	} elsif ($exceptionCode == 621) { # SERIVCE NOT AVAILABLE: [621] "Service not available"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Service not available.",
		};

# NUEVO------------------------------------------------------
# NUEVO------------------------------------------------------

# SERVICE INTRISIC ERRORS
	} elsif ($exceptionCode == 701) { # SERVICE INTERNAL ERROR: [701] "Specific errors from the BioMOBY service"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Specific errors from the BioMOBY service.",
		};

	} elsif ($exceptionCode == 702) { # OBJECT NOT FOUND: [702] "Object not found with the given input"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Object not found with the given input.",
		};

# NUEVO------------------------------------------------------
# NUEVO------------------------------------------------------

	} elsif ($exceptionCode == 721) { # INCORRECT ARTICLE NAME: [721] "The specified name of MOBYData article is wrong or does not exist"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "The specified name of MOBYData article is wrong or does not exist.",
		};

	} elsif ($exceptionCode == 722) { # INCORRECT OBJECT TYPE: [722] "Incorrect Object type from specified MOBYData article"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Incorrect Object type from specified MOBYData article.",
		};

	} elsif ($exceptionCode == 723) { # INCORRECT ARTICLENAME OBJECT: [723] "The specified article name of BioMOBY Object is wrong or does not exist"
		
		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "The specified article name of BioMOBY Object is wrong or does not exist.",
		};

	} elsif ($exceptionCode == 724) { # INCORRECT NAMESPACE OBJECT: [724] "The namespace of specified BioMOBY Object is invalid"
		
		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "The namespace of specified BioMOBY Object is invalid.",
		};

	} elsif ($exceptionCode == 731) { # INCORRECT ARTICLENAME OF SECONDARY: [731] "The specified name of secondary is wrong or does not exist"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "The specified name of secondary is wrong or does not exist.",
		};

	} elsif ($exceptionCode == 732) { # INCORRECT DATA TYPE OF SECONDARY: [732] "Incorrect data type from specified secondary article"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "Incorrect data type from specified secondary article.",
		};

	} elsif ($exceptionCode == 733) { # INCORRECT VALUE FROM SECONDARY: [733] "The value of secondary article is invalid. It is not inside of correct range"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "The value of secondary article is invalid. It is not inside of correct range.",
		};
		
	} elsif ($exceptionCode == 734) { # INCORRECT VALUE AND DEFAULT VALUE FROM SECONDARY: [734] "There is not SECONDARY value and registered SECONDARY article has not default value"

		($INB_Exception) = {
			'code' => $exceptionCode,
			'message' => "There is not SECONDARY value and registered SECONDARY article has not default value.",
		};
	}
# NUEVO------------------------------------------------------
# NUEVO------------------------------------------------------

} # End Switch

	return ($INB_Exception->{message});
	
} # End getExceptionCodeDescription

1;
