#!/usr/local/bin/perl -w

#--------------------------------------------
# showcolors - show X colors
#--------------------------------------------

@helpText = (
"",
"This Perl/Tk script permits you to look at X colors.",
"You can move the scales to change the sample color based on rgb values",
"or you can select named colors out of the listbox.",
"",
"Of the buttons at the bottom of the window:",
" Closest - finds the named color that is closest to the current sample",
" Undo - undoes the last Closest",
" Put - paints the root with the current sample",
" Filter - brings up a filter dialog for the color names",
""
);

#--------------------------------------------
# Included Library Modules
#--------------------------------------------
use Tk;


#--------------------------------------------
# Global data
#--------------------------------------------
$rgbFile = "/usr/lib/X11/rgb.txt";

%name2value = (); # Sorted by color name
%value2name = (); # Sorted by color value

$value = "#000000";
$name = "Black";

$filter = "";
$caseSensitive = 0;

$closestButtonLabel = "Closest";
$undoButtonLabel = "Undo";
$putButtonLabel = "Put";
$filterButtonLabel = "Filter...";

$applyButtonLabel = 'Apply';
$clearButtonLabel = 'Clear';
$dismissButtonLabel = 'Dismiss';


#--------------------------------------------
# Main
#--------------------------------------------


# If any arguments were given, assume that they're asking for help.
if ($ARGV[0]) {
foreach $line (@helpText) {
print "$line\n";
}
}

# Assemble the list of colors.
getColorList();


# Initialize Tk
$top = MainWindow->new();


# Put a button bar on the bottom, a sample area above that, 
# and two frames side-by-side at the top.
$buttons = $top->Frame( );
$buttons->pack( -side => 'bottom', -fill => 'x' );
$sample = $top->Frame( -height => '2c', -relief => 'ridge' );
$sample->pack( -side => 'bottom', -fill => 'x' );
$left = $top->Frame;
$left->pack( -side => 'left', -fill => 'y');
$right = $top->Frame;
$right->pack( -side => 'right', -fill => 'both', -expand => 'yes' );



# Add the buttons to the button bar.
makeButton( $buttons, \$closestButtonLabel, findClosestColor );
$undo = makeButton( $buttons, \$undoButtonLabel, undoLastClosest );
$undo->configure( -state => 'disabled' );

makeButton( $buttons, \$putButtonLabel, sub{ system "xsetroot -solid $value" });

$filterButton = makeButton( $buttons, \$filterButtonLabel, 
sub { &makeFilterDialog( $filterButton, \$filter, rereadColors ); } );


# Make a scale for each color component.
makeRGBAScales( $left, "Slide colors:" );


# Make a listbox full of color names.
$colorList = makeColorList( $right, "Double-Click on Name:", sort( keys %name2value ));


# Fall into the event loop.
MainLoop();


#--------------------------------------------
# Functions
#--------------------------------------------


#--------------------------------------------
# Assemble the color list.
#
# Sets global varibles %value2name and %name2value
#--------------------------------------------

sub getColorList
{

my ( $name, $value );



# Read the list of colors and intensities from the rgb file.
open(RGB,"< $rgbFile") || die "Couldn't open $rgbFile\n";
while (<RGB>)
{
chomp;
/^\s*(\d+)\s+(\d+)\s+(\d+)\s+(.*?)\s*$/o; # (r) (g) (b) (color name)
$value = sprintf( "#%02x%02x%02x", $1, $2, $3 );
$name = $4;

print STDERR "---> $name($value)\n";

# If a filter has been set, get rid of any color that doesn't match it.
if ( $filter )
{
my $f = ( $caseSensitive ) ? $filter : lc( $filter );
my $n = ( $caseSensitive ) ? $name : lc( $name );
if ( -1 == index( $n, $f ) ) { next }
}


# Get rid of duplicately named colors.
# Also, get rid of Transparent (which doesn't really exist)
# and all grey100's (so that we'll use White instead).
next if ( exists $value2name{$value} );
next if ( substr($name, 1) eq "ransparent" );
# next if ( length($name) > 4 && substr($name, 4) eq "100" );

$value2name{$value} = $name;
$name2value{$name} = $value;
}
close RGB;
}


#--------------------------------------------
# Makes a set of rgb scales
#
# Sets global variables $redScale, $greenScale, $blueScale
# Accesses the global variable $value.
#--------------------------------------------

sub makeRGBAScales
{
my ( $parent, $messageText ) = @_;


# Stick a message at the top.
$label = $parent->Label( -text => $messageText, -relief => "raised", -bd => 2 );
$label->pack( -side => 'top', -fill => 'x' );


# Display the current value at the bottom.
$label = $parent->Label( -textvariable => \$value, -relief => "ridge", -bd => 2 );
$label->pack( -side => 'bottom', -fill => 'x' );


# Make a scale for each color component.
$redScale = makeScale( $parent, 'red' );
$greenScale = makeScale( $parent, 'green' );
$blueScale = makeScale( $parent, 'blue' );

}


#--------------------------------------------
# Make a scale
#--------------------------------------------

sub makeScale
{
my ( $parent, $color ) = @_;
my ( $scale );

$scale = $parent->Scale( -label => substr( $color, 0, 1 ),
-from => 0, -to => 255,
-showvalue => 'no',
-orient => 'vertical',
-command => sub { &scaleCommand; } );
$scale->pack( -side => 'left', -fill => 'y' );
$scale->bind( '<1>' => sub { $scale->focus } );

$scale;
}



#--------------------------------------------
# Make a listbox with message and scrollbar
#--------------------------------------------

sub makeColorList
{
my ( $parent, $messageText, @list ) = @_;
my ( $label, $button, $listbox );


# Stick a message at the top.
$label = $parent->Label( -text => $messageText, 
-relief => "raised", -bd => 2 );
$label->pack( -fill => 'x' );


# Display the current color name at the bottom.
$label = $parent->Label( -textvariable => \$name,
-relief => "ridge", -bd => 2 );
$label->pack( -side => 'bottom', -fill => 'x' );


# Make a listbox/scrollbar pair.
$listbox = &makeList( $parent, sub { $kk = shift; &setColorByName($kk); }, 
'right', @list );

$listbox;
}


#--------------------------------------------
# Reread the list of colors (possibly with a new filter).
#--------------------------------------------

sub rereadColors
{

# Reread the list of colors.
%value2name = ();
%name2value = ();
getColorList();


# Delete the current contents of the listbox and replace it with the new list.
$colorList->delete( 0, 'end' );
$colorList->insert( 'end', sort( keys %name2value ));

}

#--------------------------------------------
# Called if the user changed a scale:
# 
# Sets global variables $value and $name
#--------------------------------------------

sub scaleCommand
{

# Get each scale's setting and recalculate the rgb value.
$value = sprintf( "#%02x%02x%02x", $redScale->get, $greenScale->get, $blueScale->get );

&newColor;

}



#--------------------------------------------
# Find the named color closest to the current value
#--------------------------------------------

sub findClosestColor
{
$top->configure( -cursor => 'watch' );
$top->update;

$prevValue = $value;

my $diff = 3*255;
my $valueBest = "#000000";
my $nameBest = "Black";
my $rBest = hex2dec( substr( $valueBest, 1, 2 ));
my $gBest = hex2dec( substr( $valueBest, 3, 2 ));
my $bBest = hex2dec( substr( $valueBest, 5, 2 ));
my $r = hex2dec( substr( $value, 1, 2 ));
my $g = hex2dec( substr( $value, 3, 2 ));
my $b = hex2dec( substr( $value, 5, 2 ));

while ( ( $v, $n ) = each %value2name )
{
my $rTemp = hex2dec( substr( $v, 1, 2 ));
my $gTemp = hex2dec( substr( $v, 3, 2 ));
my $bTemp = hex2dec( substr( $v, 5, 2 ));

my $tempDiff = dif($r,$rTemp) + dif($g,$gTemp) + dif($b,$bTemp);
if ( $tempDiff < $diff )
{
$diff = $tempDiff;
$rBest = $r;
$gBest = $g;
$bBest = $b;
$valueBest = $v;
$nameBest = $n;
}
}

setColorByName( $nameBest );

$top->configure( -cursor => 'left_ptr' );
$undo->configure( -state => 'normal' );

}


sub undoLastClosest
{
$value = $prevValue;
&newColor;
$undo->configure( -state => 'disabled' );
}


sub dif
{
my ( $p1, $p2 ) = @_;
my $result = ( $p1 < $p2 ) ? ( $p2 - $p1 ) : ( $p1 - $p2 );
}


sub setColorByName
{

# Set the global color name.
my ( $name ) = @_;


# Get its rgb value.
$value = $name2value{$name};

&newColor;
}


sub newColor
{

# If this color has a name then set it.
if ( exists( $value2name{$value} ) ) { $name = $value2name{$value}; }
else { $name = ""; }


# Set each of the scales to the proper setting.
$redScale->set( hex2dec( substr( $value, 1, 2 )));
$greenScale->set( hex2dec( substr( $value, 3, 2 )));
$blueScale->set( hex2dec( substr( $value, 5, 2 )));


# Repaint the sample area.
$sample->configure( -background => $value );
}

#--------------------------------------------
# Turn a hex number into a decimal number
#--------------------------------------------

sub hex2dec
{
my ( $hex ) = @_;
my ( $dec ) = 0;


while ( length($hex) )
{
my ( $digit ) = substr( $hex, 0, 1 );

$hex = substr( $hex, 1 );

if ( $digit eq 'a' ) { $digit = '10'; }
elsif ( $digit eq 'b' ) { $digit = '11'; }
elsif ( $digit eq 'c' ) { $digit = '12'; }
elsif ( $digit eq 'd' ) { $digit = '13'; }
elsif ( $digit eq 'e' ) { $digit = '14'; }
elsif ( $digit eq 'f' ) { $digit = '15'; }

$dec = ( $dec * 16 ) + int($digit);
}

# Return the result.
$dec;
}


#--------------------------------------------
# makeList - Create a listbox and associated scrollbar.
#--------------------------------------------

sub makeList
{
my ( $parent, $command, $side, @list ) = @_;
my ( $listbox, $scrollbar );


# Stick a scroll bar on the requested side.
$scrollbar = $parent->Scrollbar();
$scrollbar->pack( -side => $side, -fill => 'y' );


# Create a listbox and tie it to the scrollbar.
# (If the you move the listbox with button2, the scrollbar will also move.)
$listbox = $parent->Listbox( -setgrid => 'yes',
-yscrollcommand => [ 'set', $scrollbar ] );
$listbox->insert( 'end', @list );
$listbox->pack( -fill => 'both', -expand => 'yes' );


# Change the cursor when button2 is depressed to remind people that 
# they can move the listbox. 
$listbox->bind( '' => sub {
$listbox->configure( -cursor => 'spider' )}); 
$listbox->bind( '' => sub { 
$listbox->configure( -cursor => 'left_ptr' )});



# Tie the scrollbar to the listbox
# (If you move the scrollbar, the listbox will also move.)
$scrollbar->configure( -command => [ 'yview', $listbox ] );


# If an element is selected with either the mouse or the keyboard,
# do what needs to be done.
$listbox->bind( '<Double-1>' => sub {
my $selection = $listbox->get('active');
&$command( $selection );
} );
$listbox->bind( '<KeyPress-Up>' => sub {
my $selection = $listbox->get('active');
&$command( $selection );
} );
$listbox->bind( '<KeyPress-Down>' => sub {
my $selection = $listbox->get('active');
&$command( $selection );
} );
$listbox->bind( '<KeyPress-Return>' => sub {
my $selection = $listbox->get('active');
&$command( $selection );
} );



# Shift input focus to the listbox's components when the mouse touches them.
$listbox->bind( '<1>' => sub { $listbox->focus } );
$scrollbar->bind( '<1>' => sub { $scrollbar->focus } );


$listbox;
}



#--------------------------------------------
# Called when the Filter... button is pressed
#--------------------------------------------

sub makeFilterDialog
{

my ( $frame, $button, $apply, $label );

$fButton = shift;
$filterP = shift;
my $command = shift;


# Disable the filter button while the filter dialog is up.
$fButton->configure( -state => 'disabled' );


# Find the right edge of the parent window.
$geom = getRightEdge( $fButton );


# Create the dialog box
$dialog = $fButton->Toplevel( -class => 'Filter' );
$dialog->geometry( $geom );
$dialog->wm( 'protocol', 'WM_DELETE_WINDOW', sub{ &destroyDialog; } );


# Put the filter text widget at the top.
$frame = $dialog->Frame(); $frame->pack( -fill => 'x');

$label = $frame->Label( -text => "Filter:", -relief => 'flat' );
$label->pack( -side => 'left' );

$filterText = $frame->Entry( -width => 10, -relief => 'sunken', -bd => 2 );
$filterText->insert( 'end', $$filterP );
$filterText->pack( -side => 'left', -fill => 'x', -expand => 'yes' );
$filterText->bind( '' => sub { $apply->flash; 
$$filterP = $filterText->get; 
&$command } );
$filterText->focus;


# Put a case-sensitive check button under that.
$button = $dialog->Checkbutton( -text => "Case Sensitive",
-variable => \$caseSensitive );
$button->pack( -expand => 'yes', -pady => '1m' );


# Put some action buttons in the space at the bottom.
$frame = $dialog->Frame( -relief => 'groove', -bd => 5 );
$frame->pack( -side => 'left', -padx => '1m', -pady => '1m', -expand => 'yes' );
$apply = makeButton( $frame, \$applyButtonLabel, 
sub { $$filterP = $filterText->get; &$command } );


makeButton( $dialog, \$clearButtonLabel,
sub { $$filterP = ""; $filterText->delete( 0, 'end' );} );

makeButton( $dialog, \$dismissButtonLabel, destroyDialog );

}


#--------------------------------------------
# Destroy the dialog window and re-enable the filter button.
#--------------------------------------------

sub destroyDialog
{
destroy $dialog; 
$fButton->configure( -state => 'normal' );
}



#--------------------------------------------
# Find the right edge of the parent window.
#--------------------------------------------

sub getRightEdge
{
my ( $w ) = @_;
my $geom;

# Get the geometry of the eldest parent.
while ( $w )
{
$geom = $w->winfo( 'geometry' );
$w = $w->winfo( 'parent' ); 
}

# Split it into its component parts.
my ( $xCoord, $yCoord, $xSize, $ySize ) = parseGeometry( $geom );


# Add the x size and coordinate to get the right edge.
my $x = $xCoord + $xSize;


# Find the midpoint of that edge for the y coordinate.
my $y = $yCoord + int($ySize/2);


# Assemble x and y into a geometry statement (coordinates only, no size).
$geom = "+$x+$y";
}


sub parseGeometry
{
my ( $geom ) = @_;
my ( $size, $xCoord, $yCoord, $xSize, $ySize );


# Split the geometry string into its component parts.
my ( $size, $xCoord, $yCoord ) = split '\+', $geom ;
my ( $xSize, $ySize ) = split 'x', $size;


# Return the components.
return ( $xCoord, $yCoord, $xSize, $ySize );
}





#--------------------------------------------
# Make a button
#--------------------------------------------

sub makeButton
{
my ( $parent, $text, $command ) = @_;
my ( $button );

$button = $parent->Button( -textvariable => $text, -command => sub { &$command } );
$button->bind( '' => sub { &$command } );
$button->pack( -side => 'left', -padx => '1m', -pady => '1m', -expand => 'yes' );

$button;
}
