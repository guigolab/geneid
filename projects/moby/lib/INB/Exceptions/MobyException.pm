package INB::Exceptions::MobyException;

# Perl pragma to restrict unsafe constructs
use strict;

# Issue warnings about suspicious programming.
use warnings;

use Carp qw(croak);

use INB::Exceptions::MobyExceptionCodes;

###############
# Constructor #
###############

sub new {
	# Get parameters
	my ($caller, %args) = @_;
	
	my ($caller_is_obj) = ref($caller);
	my ($class) = $caller_is_obj || $caller;
	my ($self) = {};

	# Add attributes related to mobyData where has been generating exception.
	$self->{queryID} = (exists($args{queryID}) && defined($args{queryID})) ? $args{queryID} : '';
	$self->{refElement} = (exists($args{refElement}) && defined($args{refElement})) ? $args{refElement} : '';

	# Add attribute of exception code
	if (exists($args{code}) && defined($args{code})) {
		# Check if input code is correct 
		my ($standardMessage) = MobyExceptionCodes::getExceptionCodeDescription($args{code});
	
		# If standard message is defined that means, the exception code is correct
		if (defined($standardMessage)) {
			$self->{code} = $args{code};
		} else {
			$self->{code} = 0;
		}
	}
	
	# Add attribute of exception "dynamic" message
	$self->{message} = (exists($args{message}) && defined($args{message})) ?  $args{message} : undef;

	# Add attribute of type of exception
	if (exists($args{type}) && defined($args{type})) {
		my ($type) = lc($args{type}); # Converts all characters in the string to lower case
		if ($type eq 'info' || $type eq 'warning' || $type eq 'error') {
			$self->{type} = $type; #  (undef | error | warning | info)
		} else {
			$self->{type} = undef;
		}
	} else {
		$self->{type} = undef;
	}

#	$self->{level} = (exists($args{level}) && defined($args{level})) ?  $args{level} : undef; # (undef | mobySimple | mobyCollection | mobyData)

	# The magic object creation command
	bless ($self, $class);

	# Returns a new 'MobyException' object
	return $self;
}

#################
# Method bodies #
#################

# Return exception code
sub getExceptionCode {
	# Get parameters
	my($self) = shift;
	
	croak("This is an instance method!")  unless(ref($self));
	
	return $self->{code};
}

# Return exception message
sub getExceptionMessage {
	# Get parameters
	my($self) = shift;

	my ($exceptionMessage) = undef;

	croak("This is an instance method!")  unless(ref($self));
	
	# Get standard exception message from given code
	my ($standardMessage) = MobyExceptionCodes::getExceptionCodeDescription($self->{code});

	# If standard message is not defined that means, the exception code is wrong => return undef
	return undef unless(defined($standardMessage));

	# User could add dynamic message into satandard exception message
	($exceptionMessage) = (defined($self->{message})) ? $standardMessage.$self->{message} : $standardMessage;
	
	return $exceptionMessage;
}

# Return type of exception 
sub getExceptionType {
	# Get parameters
	my($self) = shift;
	
	croak("This is an instance method!")  unless(ref($self));
	
	# type of exception is returned, if defined type of exception; otherwise returns undef
	return $self->{type};
}

# Assign exception code
sub setExceptionCode {
	# Get parameters
	my($self, $code) = @_;

	croak("This is an instance method!")  unless(ref($self));

	# Check if input code is correct 
	my ($standardMessage) = MobyExceptionCodes::getExceptionCodeDescription($code);
	
	# If standard message is defined that means, the exception code is correct
	if (defined($standardMessage)) {
		$self->{code} = $code;
	} else {
		croak("input argument not defined"); 
	}
}

# Assign exception message
sub setExceptionMessage {
	# Get parameters
	my($self, $message) = @_;

	croak("This is an instance method!")  unless(ref($self));

	if (defined($message)) {
		$self->{message} = $message; 
	} else {
		croak("input argument not defined"); 
	}
}

# Assign type of exception to attribute of class
sub setExceptionType {
	# Get parameters
	my($self, $type) = @_;

	croak("This is an instance method!")  unless(ref($self));

	# Input type has to be defined and to include within range of values
	if (defined($type)) {
		my ($type) = lc($type); # Converts all characters in the string to lower case
		
		# Input type has to be included within range of values
		if ($type eq 'info' || $type eq 'warning' || $type eq 'error') {
			$self->{type} = $type;
		} else {
			croak("input argument not defined");
		}
	}
}

# Return xml block that will be the exception response (error, warning or information)
sub retrieveExceptionResponse {
	# Get parameters
	my($self)=shift;
	my ($exceptionResponse) = undef;
	
	croak("This is an instance method!")  unless(ref($self));
	
	# Get standard exception message from given code
	my ($standardMessage) = MobyExceptionCodes::getExceptionCodeDescription($self->{code});

	#return undef unless(defined($standardMessage));
	croak("code of exception is wrong or does not exists") unless(defined($standardMessage));

	# User could add dynamic message into satandard exception message
	my ($exceptionMessage) = (defined($self->{message})) ? $standardMessage.$self->{message} : $standardMessage;

	if (defined($self->{type}) && ($self->{type} eq 'info' || $self->{type} eq 'warning' || $self->{type} eq 'error')) {
		if ($self->{type} eq 'error') {
			$exceptionResponse = "<mobyException refQueryID='".$self->{queryID}."' refElement='".$self->{refElement}."' severity='error'>\n\t<exceptionCode>".$self->{code}."</exceptionCode>\n\t<exceptionMessage>$exceptionMessage</exceptionMessage>\n</mobyException>";
		
		} elsif ($self->{type} eq 'warning') {
			$exceptionResponse = "<mobyException refQueryID='".$self->{queryID}."' refElement='".$self->{refElement}."' severity='warning'>\n\t<exceptionCode>".$self->{code}."</exceptionCode>\n\t<exceptionMessage>$exceptionMessage</exceptionMessage>\n</mobyException>";
			
		} elsif ($self->{type} eq 'info') {
			$exceptionResponse = "<mobyException refQueryID='".$self->{queryID}."' refElement='".$self->{refElement}."' severity='information'>\n\t<exceptionCode>".$self->{code}."</exceptionCode>\n\t<exceptionMessage>$exceptionMessage</exceptionMessage>\n</mobyException>"; 
		
		} else {
			croak("type of exception is wrong or does not exists");
		}
	} else {
		croak("type of exception is wrong or does not exists");
	}
	
	return $exceptionResponse;
}

# Return xml block of one empty MobyData
sub retrieveEmptyMobyData {
	# Get parameters
	my($self) = shift;
	
	croak("This is an instance method!")  unless(ref($self));
	
	return "<moby:Data moby:queryID='".$self->{queryID}."' />";
}

# Return xml block of one empty simple MobyArticle
sub retrieveEmptyMobySimple {
	# Get parameters
	my($self, $outputArticle)= @_;
	
	croak("This is an instance method!")  unless(ref($self));
	
	return "<moby:Simple moby:articleName='$outputArticle' />";
}

# Return xml block of one empty collection MobyArticle
sub retrieveEmptyMobyCollection {
	# Get parameters
	my($self, $outputArticle) = @_;
	
	croak("This is an instance method!")  unless(ref($self));
	
	return "<moby:Collection moby:articleName='$outputArticle' />";
}

# Return MOBYData inserting MOBYArticles that has been giving by input
sub embedMOBYArticlesIntoMOBYData {
	# Get parameters
	my($self, $outputMOBYArticles) = @_;

	# Returns MOBYData response
	return "<moby:mobyData moby:queryID='".$self->{queryID}."'>$outputMOBYArticles</moby:mobyData>";

}

#--------------------------------------------------------------------------------------------
#---------------------------------OLD VERSION 2.0--------------------------------------------
#--------------------------------------------------------------------------------------------


# # Return xml block of ERROR exception response
# # If the exception code is wrong, then the method returns undef
# sub retrieveErrorExceptionResponse {
# 	# Get parameters
# 	my($self) = shift;
# 
# 	my ($exceptionResponse) = undef;
# 	
# 	croak("This is an instance method!")  unless(ref($self));
# 	
# 	# Get standard exception message from given code
# 	my ($standardMessage) = MobyExceptionCodes::getExceptionCodeDescription($self->{code});
# 
# 	# If standard message is not defined that means, the exception code is wrong => return undef
# 	return undef unless(defined($standardMessage));
# 
# 	# User could add dynamic message into satandard exception message
# 	my ($exceptionMessage) = (defined($self->{message})) ? $standardMessage.$self->{message} : $standardMessage;
# 
# 	return "<mobyException refQueryID='".$self->{queryID}."' refElement='".$self->{refElement}."' severity='error'>\n\t<exceptionCode>".$self->{code}."</exceptionCode>\n\t<exceptionMessage>".$exceptionMessage."</exceptionMessage>\n</mobyException>";
# }
# 
# # Return xml block of WARNING exception response
# # If the exception code is wrong, then the method returns undef
# sub retrieveWarningExceptionResponse {
# 	# Get parameters
# 	my($self) = shift;
# 
# 	my ($exceptionResponse) = undef;
# 	
# 	croak("This is an instance method!")  unless(ref($self));
# 	
# 	# Get standard exception message from given code
# 	my ($standardMessage) = MobyExceptionCodes::getExceptionCodeDescription($self->{code});
# 
# 	# If standard message is not defined that means, the exception code is wrong => return undef
# 	return undef unless(defined($standardMessage));
# 
# 	# User could add dynamic message into satandard exception message
# 	my ($exceptionMessage) = (defined($self->{message})) ? $standardMessage.$self->{message} : $standardMessage;
# 
# 	return "<mobyException refQueryID='".$self->{queryID}."' refElement='".$self->{refElement}."' severity='warning'>\n\t<exceptionCode>".$self->{code}."</exceptionCode>\n\t<exceptionMessage>".$exceptionMessage."</exceptionMessage>\n</mobyException>";
# }
# 
# # Return xml block of INFORMATION exception response
# # If the exception code is wrong, then the method returns undef
# sub retrieveInfoExceptionResponse {
# 	# Get parameters
# 	my($self) = shift;
# 
# 	my ($exceptionResponse) = undef;
# 	
# 	croak("This is an instance method!")  unless(ref($self));
# 	
# 	# Get standard exception message from given code
# 	my ($standardMessage) = MobyExceptionCodes::getExceptionCodeDescription($self->{code});
# 
# 	# If standard message is not defined that means, the exception code is wrong => return undef
# 	return undef unless(defined($standardMessage));
# 
# 	# User could add dynamic message into satandard exception message
# 	my ($exceptionMessage) = (defined($self->{message})) ? $standardMessage.$self->{message} : $standardMessage;
# 
# 	return "<mobyException refQueryID='".$self->{queryID}."' refElement='".$self->{refElement}."' severity='information'>\n\t<exceptionCode>".$self->{code}."</exceptionCode>\n\t<exceptionMessage>".$exceptionMessage."</exceptionMessage>\n</mobyException>";
# }

#--------------------------------------------------------------------------------------------
#---------------------------------OLD VERSION------------------------------------------------
#--------------------------------------------------------------------------------------------

# # Assign type of exception
# sub setExceptionType {
# 	# Get parameters
# 	my($self, $type) = @_;
# 
# 	croak("This is an instance method!")  unless(ref($self));
# 
# 	if (defined($type) && ($type eq 'info' || $type eq 'warning' || $type eq 'error')) {
# 		$self->{type} = $type; 
# 	} else {
# 		croak("input argument not well-know define"); 
# 	}
# }
# 
# # Assign type of exception
# sub setExceptionLevel {
# 	# Get parameters
# 	my($self, $type) = @_;
# 
# 	croak("This is an instance method!")  unless(ref($self));
# 
# 	if (defined($type) && ($type eq 'info' || $type eq 'warning' || $type eq 'error')) {
# 		$self->{type} = $type; 
# 	} else {
# 		croak("input argument not well-know define"); 
# 	}
# }

# # Return type of exception (info, warning or error)
# sub getExceptionType {
# 	# Get parameters
# 	my($self)=shift;
# 	
# 	croak("This is an instance method!")  unless(ref($self));
# 	
# 	return $self->{type};
# }
# 
# # Return "Moby level" where exception was produced
# sub getExceptionLevel {
# 	# Get parameters
# 	my($self)=shift;
# 	
# 	croak("This is an instance method!")  unless(ref($self));
# 	
# 	return $self->{level};
# }

# # Return xml block that will be the empty Moby response
# sub retrieveEmptyMobyResponse {
# 	# Get parameters
# 	my($self)=shift;
# 	
# 	croak("This is an instance method!")  unless(ref($self));
# 	
# 	my ($mobyResponse) = undef;
# 
# 	if (defined($self->{level}) && ($self->{level} eq 'mobySimple' || $self->{level} eq 'mobyCollection' || $self->{level} eq 'mobyData')) {
# 		if ($self->{level} eq 'mobySimple') {
# 			$mobyResponse =  "\n\t<moby:Simple moby:articleName='".$self->{refElement}."' />\n";
# 		} elsif ($self->{level} eq 'mobyCollection') {
#  			$mobyResponse =  "\n\t<moby:Collection moby:articleName='".$self->{refElement}."' />\n";
# 		} elsif ($self->{level} eq 'mobyData') {
#  			$mobyResponse =  "\n\t<moby:Data moby:queryID='".$self->{queryID}."' />\n";
# 		}
# 	}
# 	
# 	return $mobyResponse;
# }

# # Return xml block that will be the exception response
# sub retrieveExceptionResponse {
# 	# Get parameters
# 	my($self)=shift;
# 	my ($exceptionResponse) = undef;
# 	
# 	croak("This is an instance method!")  unless(ref($self));
# 	
# 	# Get standard exception message from given code
# 	my ($standardMessage) = MobyExceptionCodes::getExceptionCodeDescription($self->{code});
# 
# 	return undef unless(defined($standardMessage));
# 
# 	# User could add dynamic message into satandard exception message
# 	my ($exceptionMessage) = (defined($self->{message})) ? $standardMessage.$self->{message} : $standardMessage;
# 
# 	if (defined($self->{type}) && ($self->{type} eq 'info' || $self->{type} eq 'warning' || $self->{type} eq 'error')) {
# 		if ($self->{type} eq 'error') {
# 			$exceptionResponse = "\n<mobyException refQueryID='".$self->{queryID}."' refElement='".$self->{refElement}."' severity='error'>\n\t<exceptionCode>".$self->{code}."</exceptionCode>\n\t<exceptionMessage>".$exceptionMessage."</exceptionMessage>\n</mobyException>\n";
# 		
# 		} elsif ($self->{type} eq 'warning') {
# 			$exceptionResponse = "\n<mobyException refQueryID='".$self->{queryID}."' refElement='".$self->{refElement}."' severity='warning'>\n\t<exceptionCode>".$self->{code}."</exceptionCode>\n\t<exceptionMessage>".$exceptionMessage."</exceptionMessage>\n</mobyException>\n";
# 			
# 		} elsif ($self->{type} eq 'info') {
# 			$exceptionResponse = "\n<mobyException refQueryID='".$self->{queryID}."' refElement='".$self->{refElement}."' severity='information'>\n\t<exceptionCode>".$self->{code}."</exceptionCode>\n\t<exceptionMessage>".$exceptionMessage."</exceptionMessage>\n</mobyException>\n"; 
# 		
# 		}
# 	}
# 	
# 	return $exceptionResponse;
# }

sub DESTROY {}

1;

