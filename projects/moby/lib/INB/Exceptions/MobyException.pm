#$Id: MobyException.pm ,v 1.2 
# Created: 26-01-2006
# Updated: 01-02-2006

=head1 NAME

MobyException

=head1 DESCRIPTION

Class that contains exception instance and exception methods

=head1 AUTHORS

Jose Manuel Rodriguez Carrasco -jmrc@cnb.uam.es- (INB-CNB)

=head1 METHODS

=cut

# Name of the package
package INB::Exceptions::MobyException;

# Perl pragma to restrict unsafe constructs
use strict;

# Issue warnings about suspicious programming.
use warnings;

use Carp qw(croak);

use INB::Exceptions::MobyExceptionCodes;

=head2 new

 Function  : Create new instance of exception class
 Args      : - querID from one MobyData assign to exception 
 	     - refElement, reference to articleName 
	     - Exception Code
	     - Exception Message
	     - Type of exception: error, information, or warning
 Returns   : - Exception Instace 

=cut

sub new {
	# Get parameters
	my ($caller, %args) = @_;
	
	my ($caller_is_obj) = ref($caller);
	my ($class) = $caller_is_obj || $caller;
	my ($self) = {};

	# Add attributes related to mobyData where has been generating exception.
	$self->{queryID} = (exists($args{queryID}) && defined($args{queryID})) ? $args{queryID} : '';
	$self->{refElement} = (exists($args{refElement}) && defined($args{refElement})) ? $args{refElement} : '';


# UPDATE: 19-01-2006-----------------
# REASON: DUE TO REPORTING INFORMATION MESSAGE. WE DON'T HAVE TO CHECK INPUT CODE AND MESSAGE

# 	# Add attribute of exception code
# 	if (exists($args{code}) && defined($args{code})) {
# 		# Check if input code is correct 
# 		my ($standardMessage) = INBServices::lib::Exceptions::MobyExceptionCodes::getExceptionCodeDescription($args{code});
# 	
# 		# If standard message is defined that means, the exception code is correct
# 		if (defined($standardMessage)) {
# 			$self->{code} = $args{code};
# 		} else {
# 			$self->{code} = 0;
# 		}
# 	} 
	# Add attribute of exception code
	if (exists($args{code}) && defined($args{code})) {
		$self->{code} = $args{code};
	} else {
		$self->{code} = 0;
	}
# UPDATE: 19-01-2006-----------------

	# Add attribute of exception "dynamic" message
	$self->{message} = (exists($args{message}) && defined($args{message})) ?  $args{message} : undef;

	# Add attribute of type of exception
	if (exists($args{type}) && defined($args{type})) {
		my ($type) = lc($args{type}); # Converts all characters in the string to lower case
		if ($type eq 'information' || $type eq 'warning' || $type eq 'error') {
			$self->{type} = $type; #  (undef | error | warning | information)
		} else {
			$self->{type} = undef;
		}
	} else {
		$self->{type} = undef;
	}

	# The magic object creation command
	bless ($self, $class);

	# Returns a new 'MobyException' object
	return $self;
}

#################
# Method bodies #
#################

=head2 getExceptionCode

 Function  : Return exception code
 Args      : 
 Returns   : - Integer: Exception Code

=cut

sub getExceptionCode {
	# Get parameters
	my($self) = shift;
	
	croak("This is an instance method!")  unless(ref($self));
	
	return $self->{code};
} # End getExceptionCode

=head2 getExceptionMessage

 Function  : Return exception message
 Args      : 
 Returns   : - String: Exception message 

=cut

sub getExceptionMessage {
	# Get parameters
	my($self) = shift;

	my ($exceptionMessage) = undef;

	croak("This is an instance method!")  unless(ref($self));
	
	# Get standard exception message from given code
	my ($standardMessage) = INB::Exceptions::MobyExceptionCodes::getExceptionCodeDescription($self->{code});

	# If standard message is not defined that means, the exception code is wrong => return undef
	return undef unless(defined($standardMessage));

	# User could add dynamic message into satandard exception message
	($exceptionMessage) = (defined($self->{message})) ? $standardMessage.$self->{message} : $standardMessage;
	
	return $exceptionMessage;
} # End getExceptionMessage

=head2 getExceptionType

 Function  : Return type of exception
 Args      : 
 Returns   : - String (error, information, warning): Exception type of exception

=cut

sub getExceptionType {
	# Get parameters
	my($self) = shift;
	
	croak("This is an instance method!")  unless(ref($self));
	
	# type of exception is returned, if defined type of exception; otherwise returns undef
	return $self->{type};
} # End getExceptionType

=head2 setExceptionCode

 Function  : Assign exception code
 Args      : - Integer: Exception Code
 Returns   : 

=cut

sub setExceptionCode {
	# Get parameters
	my($self, $code) = @_;

	croak("This is an instance method!")  unless(ref($self));

# UPDATE: 19-01-2006-----------------
# REASON: DUE TO REPORTING INFORMATION MESSAGE. WE DON'T HAVE TO CHECK INPUT CODE AND MESSAGE
# 	# Check if input code is correct 
# 	my ($standardMessage) = INBServices::lib::Exceptions::MobyExceptionCodes::getExceptionCodeDescription($code);
# 	
# 	# If standard message is defined that means, the exception code is correct
# 	if (defined($standardMessage)) {
# 		$self->{code} = $code;
# 	} else {
# 		croak("input argument not defined"); 
# 	}
	if (defined($code)) {
		$self->{code} = $code;
	} else {
		croak("input argument not defined"); 
	}
# UPDATE: 19-01-2006-----------------
} # End setExceptionCode

=head2 setExceptionMessage

 Function  : Assign exception message
 Args      : - String: Exception message
 Returns   : 

=cut

sub setExceptionMessage {
	# Get parameters
	my($self, $message) = @_;

	croak("This is an instance method!")  unless(ref($self));

	if (defined($message)) {
		$self->{message} = $message; 
	} else {
		croak("input argument not defined"); 
	}
} # End setExceptionMessage

=head2 setExceptionType

 Function  : Assign type of exception to attribute of class
 Args      : - String (error, information, warning): type of exception
 Returns   : 

=cut

sub setExceptionType {
	# Get parameters
	my($self, $type) = @_;

	croak("This is an instance method!")  unless(ref($self));

	# Input type has to be defined and to include within range of values
	if (defined($type)) {
		my ($type) = lc($type); # Converts all characters in the string to lower case
		
		# Input type has to be included within range of values
		if ($type eq 'information' || $type eq 'warning' || $type eq 'error') {
			$self->{type} = $type;
		} else {
			croak("input argument not defined");
		}
	}
} # End setExceptionType

=head2 retrieveExceptionResponse

 Function  : Return xml block that will be the exception response (error, warning or information)
 Args      : 
 Returns   : - xml block of exception response

=cut

sub retrieveExceptionResponse {
	# Get parameters
	my($self)=shift;
	my ($exceptionResponse) = undef;
	
	croak("This is an instance method!")  unless(ref($self));

# UPDATE: 19-01-2006-----------------
# REASON: DUE TO REPORTING INFORMATION MESSAGE. WE DON'T HAVE TO CHECK INPUT CODE AND MESSAGE

# 	# Get standard exception message from given code
# 	my ($standardMessage) = INBServices::lib::Exceptions::MobyExceptionCodes::getExceptionCodeDescription($self->{code});
# 
# 	#return undef unless(defined($standardMessage));
# 	croak("code of exception is wrong or does not exists") unless(defined($standardMessage));
# 
# 	# User could add dynamic message into satandard exception message
# 	my ($exceptionMessage) = (defined($self->{message})) ? $standardMessage.$self->{message} : $standardMessage;
# 
# 	if (defined($self->{type}) && ($self->{type} eq 'information' || $self->{type} eq 'warning' || $self->{type} eq 'error')) {
# 		if ($self->{type} eq 'error') {
# 			$exceptionResponse = "<mobyException refQueryID='".$self->{queryID}."' refElement='".$self->{refElement}."' severity='error'>\n\t<exceptionCode>".$self->{code}."</exceptionCode>\n\t<exceptionMessage>$exceptionMessage</exceptionMessage>\n</mobyException>";
# 		
# 		} elsif ($self->{type} eq 'warning') {
# 			$exceptionResponse = "<mobyException refQueryID='".$self->{queryID}."' refElement='".$self->{refElement}."' severity='warning'>\n\t<exceptionCode>".$self->{code}."</exceptionCode>\n\t<exceptionMessage>$exceptionMessage</exceptionMessage>\n</mobyException>";
# 			
# 		} elsif ($self->{type} eq 'information') {
# 			$exceptionResponse = "<mobyException refQueryID='".$self->{queryID}."' refElement='".$self->{refElement}."' severity='information'>\n\t<exceptionCode>".$self->{code}."</exceptionCode>\n\t<exceptionMessage>$exceptionMessage</exceptionMessage>\n</mobyException>"; 
# 		
# 		} else {
# 			croak("type of exception is wrong or does not exists");
# 		}
# 	} else {
# 		croak("type of exception is wrong or does not exists");
# 	}
	
	# If corresponds to free text message, we don't mind the code or message => User is free to insert what he wants
	if (defined($self->{type}) && ($self->{type} eq 'information')) {

		# Chek if is defined information message
		my ($infoMessage) = (defined($self->{message})) ? $self->{message} : '';

		$exceptionResponse = "<mobyException refQueryID='".$self->{queryID}."' refElement='".$self->{refElement}."' severity='information'>\n\t<exceptionCode>".$self->{code}."</exceptionCode>\n\t<exceptionMessage>$infoMessage</exceptionMessage>\n</mobyException>";

	} else {

		# Get standard exception message from given code
		my ($standardMessage) = INB::Exceptions::MobyExceptionCodes::getExceptionCodeDescription($self->{code});

		#return undef unless(defined($standardMessage));
		croak("code of exception is wrong or does not exists") unless(defined($standardMessage));

		# User could add dynamic message into satandard exception message
		my ($exceptionMessage) = (defined($self->{message})) ? $standardMessage.$self->{message} : $standardMessage;

		if (defined($self->{type}) && ($self->{type} eq 'warning' || $self->{type} eq 'error')) {
			if ($self->{type} eq 'error') {
				$exceptionResponse = "<mobyException refQueryID='".$self->{queryID}."' refElement='".$self->{refElement}."' severity='error'>\n\t<exceptionCode>".$self->{code}."</exceptionCode>\n\t<exceptionMessage>$exceptionMessage</exceptionMessage>\n</mobyException>";
		
			} elsif ($self->{type} eq 'warning') {
				$exceptionResponse = "<mobyException refQueryID='".$self->{queryID}."' refElement='".$self->{refElement}."' severity='warning'>\n\t<exceptionCode>".$self->{code}."</exceptionCode>\n\t<exceptionMessage>$exceptionMessage</exceptionMessage>\n</mobyException>";

			} else {
				croak("type of exception is wrong or does not exists");
			}
	
		} else {
			croak("type of exception is wrong or does not exists");
		}
	}
# UPDATE: 19-01-2006-----------------

	return $exceptionResponse;
} # End retrieveExceptionResponse

=head2 retrieveEmptyMobyData

 Function  : Return xml block of one empty MobyData
 Args      : 
 Returns   : - xml block of one empty MobyData

=cut

sub retrieveEmptyMobyData {
	# Get parameters
	my($self) = shift;
	
	croak("This is an instance method!")  unless(ref($self));
	
	return "<moby:mobyData moby:queryID='".$self->{queryID}."' />";
} # End retrieveEmptyMobyData

=head2 retrieveEmptyMobySimple

 Function  : Return xml block of one empty simple MobyArticle
 Args      : - String: name of output article
 Returns   : - xml block of one empty simple MobyArticle

=cut

sub retrieveEmptyMobySimple {
	# Get parameters
	my($self, $outputArticle)= @_;
	
	croak("This is an instance method!")  unless(ref($self));
	
	return "<moby:Simple moby:articleName='$outputArticle' />";
} # End retrieveEmptyMobySimple

=head2 retrieveEmptyMobyCollection

 Function  : Return xml block of one empty collection MobyArticle
 Args      : - String: name of output article
 Returns   : - xml block of one empty collection MobyArticle

=cut

sub retrieveEmptyMobyCollection {
	# Get parameters
	my($self, $outputArticle) = @_;
	
	croak("This is an instance method!")  unless(ref($self));
	
	return "<moby:Collection moby:articleName='$outputArticle' />";
} # End retrieveEmptyMobyCollection

=head2 embedMOBYArticlesIntoMOBYData

 Function  : Return MobyData inserting MobyArticles that has been giving by input
 Args      : - xml block which contains MobyArticles
 Returns   : - xml block of MobyData

=cut

sub embedMOBYArticlesIntoMOBYData {
	# Get parameters
	my($self, $outputMOBYArticles) = @_;

	# Returns MOBYData response
	return "<moby:mobyData moby:queryID='".$self->{queryID}."'>$outputMOBYArticles</moby:mobyData>";

}

=head2 embedExceptionsIntoServiceNotes

 Function  : Return ServiceNotes tag inserting MobyExceptions that has been giving by input
 Args      : - xml block which contains MobyExceptions
 Returns   : - xml block of MobyData

=cut

sub embedExceptionsIntoServiceNotes {
	# Get parameters
	my($self, $outputMOBYExceptions) = @_;

	# Returns MOBYData response
	return "<serviceNotes>$outputMOBYExceptions</serviceNotes>";

}

# UPDATE: 01-02-2006-----------------
# REASON: Method that returns empty mobyStatus during asynchronous callings
# Modified by jmrc

=head2 retrieveEmptyMobyStatus

 Function  : Return xml block of one empty MobyStatus
 Args      : 
 Returns   : - xml block of one empty MobyStatus

=cut

sub retrieveEmptyMobyStatus {
	# Get parameters
	my($self) = shift;
	
	croak("This is an instance method!")  unless(ref($self));
	
	return "<moby:mobyStatus moby:queryID='".$self->{queryID}."' />";
} # End retrieveEmtyMobyStatus
# UPDATE: 01-02-2006-----------------


sub DESTROY {}

1;

