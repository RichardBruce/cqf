#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use LWP::Simple;
use Date::Parse;
use Getopt::Long;

# RIC extension to google code
my %ric_extension_to_google_code =
(
    'HK' => 'HKG'
);

my $help            = 0;
my $verbose         = 0;
my $ric             = "";
my $bbg             = "";
my $date_from_str   = "";
my $date_to_str     = "";
my $days_back       = -1;

# Argument handling
GetOptions (
    "from:s"    => \$date_from_str,
    "to:s"      => \$date_to_str,
    "days:i"    => \$days_back,
    "ric:s"     => \$ric,
    "bbg:s"     => \$bbg,
    "v"         => \$verbose,
    "help"      => \$help
);

my $secs_now    = time();
my $secs_from   = str2time( $date_from_str );
my $secs_to     = 0;
if ( $days_back > -1 )
{
    print "Info: Creating to by subtracting $days_back days.\n";
    $secs_to = $secs_from - ( $days_back * 24 * 60 * 60 );
}
else
{
    $secs_to = str2time( $date_to_str );
}

# Checks
if ( not defined $secs_from )
{
    print "Error: From date ($date_from_str) not correctly specified\n";
    exit 1;
}

if ( not defined $secs_to )
{
    print "Error: To date ($date_to_str) not correctly specified\n";
    exit 1;
}

if ( $secs_from < $secs_to )
{
    print "Error: From time ($date_from_str) less than time to ($date_to_str).\n";
    exit 1;
}

if ( $secs_now < $secs_from )
{
    print "Error: From time ($secs_from) greater than now ($secs_now).\n";
    exit 1;
}

# Range of dates to fetch
my $secs_offset = $secs_now - $secs_from;
my $secs_period = $secs_from - $secs_to;
my $days_offset = $secs_offset / ( 24 * 60 * 60 );
my $days_period = $secs_period / ( 24 * 60 * 60 );

print "$secs_from $secs_to\n";
exit;

# Connect to DB
my $username    = 'bruce';
my $pass        = 'BU**4rit';
my $db          = 'market_data';
my $dbh         = DBI->connect( "dbi:mysql:$db", $username, $pass, { 'PrintError' => 1, 'RaiseError' => 1 } );

# Prepare SQL
my $insert          = 'INSERT INTO stocks(id, price, day, volume, high, low) VALUE (?, ?, STR_TO_DATE(?, \'%M %d,%Y\'), ?, ?, ?)';
my $insert_handle   = $dbh->prepare($insert);

# Check the codes specified exists and are consistent
my $exists      = 0;
my $consistent  = 0;
if ( $ric and $bbg )
{
    print "Info: Checking code consistency.\n" unless !$verbose;
    $consistent = check_id_consistency($dbh, $ric, $bbg);
}

if ( $ric )
{
    print "Info: Checking ric code exists.\n" unless !$verbose;
    $exists = check_id_exists($dbh, "ric", $ric);
}
elsif ( $bbg )
{
    print "Info: Checking bbg code exists.\n" unless !$verbose;
    $exists = check_id_exists($dbh, "bbg", $bbg);
}
else
{
    print "Error: A BBG or RIC code must be specified.\n";
    exit;
}

# Tidy up any mess possible
if (( $exists == 0 ) and $ric and $bbg )
{
    # Add id
    print "Info: Adding id for bbg: $bbg, ric: $ric\n" unless !$verbose;
    create_new_id_for_codes($dbh, $ric, $bbg);
}
elsif ( $consistent == 0 )
{
    print "Error: Inconsistent codes.\n";
    exit 1;
}

# Get all the codes
my $id = get_id_from_ric_or_bbg($dbh, $ric, $bbg);
$ric = get_ric_from_id($dbh, $id);
$bbg = get_bbg_from_id($dbh, $id);

# Get the web page to parse
my ( $ric_num, $ric_ext ) = split(/\./, $ric);
my $google_code = $ric_extension_to_google_code{$ric_ext};
my $page_base = "http://www.google.com/finance/historical?q=$google_code%3A$ric_num&start=0&num=200";
my $webpage = get($page_base);
$webpage =~ s/\n/ /g;

# Find the historic data table
if ($webpage =~ m/<table class=\"gf-table historical_price\">(.*)<\/table>/i)
{
    # Split the table into rows
    my @rows = split(/<tr[\s\w\=\"]*>/, $1);
    
    # Split the header into columns
    shift( @rows );
    my $header = shift( @rows );
    my @headers = split(/<th[\s\w\=\"]*>/, $header);
    
    # Get the index of the required column names
    my ( $day_idx, $price_idx, $volume_idx, $high_idx, $low_idx );
    my $idx = 0;
    foreach ( @headers )
    {
        if ( $_ =~ m/date/i )
        {
            $day_idx = $idx;
        }
        elsif ( $_ =~ m/high/i )
        {
            $high_idx = $idx;
        }
        elsif ( $_ =~ m/low/i )
        {
            $low_idx = $idx;
        }
        elsif ( $_ =~ m/close/i )
        {
            $price_idx = $idx;
        }
        elsif ( $_ =~ m/volume/i )
        {
            $volume_idx = $idx;
        }
        $idx++;
    }
    
    # Get the date from each row
    foreach ( @rows )
    {
        # Split the row
        my @data = split(/<td[\s\w\=\"]*>/, $_);
        
        # Get the data by index
        my $day     = clean_string($data[$day_idx    ]);
        my $price   = clean_string($data[$price_idx  ]);
        my $volume  = clean_string($data[$volume_idx ]);
        my $high    = clean_string($data[$high_idx   ]);
        my $low     = clean_string($data[$low_idx    ]);
        
        # Add to data base
        insert_stock_data_into_db($insert_handle, $id, $price, $day, $volume, $high, $low);
    }
}
else
{
    print $webpage;
    print "Error: Couldn't find historic data table on page.\n";
}

# Check the DB
my $sql         = 'SELECT * FROM stocks WHERE id = ? ORDER BY day';
my $sql_handle  = $dbh->prepare($sql);
$sql_handle->execute($id);

# Print output
my @data;
while (@data = $sql_handle->fetchrow_array())
{
    print join(", ",@data);
    print "\n";
}
    
exit 0;

# Runs SQL on DB to insert
sub insert_stock_data_into_db
{
    my ( $ih, $id, $price, $day, $volume, $high, $low ) = @_;
    $ih->execute($id, $price, $day, $volume, $high, $low);
}

# Check an id exists for a ric or bbg code
sub check_id_exists
{
    my ( $dbh, $code, $id ) = @_;
    my $sql         = "SELECT * FROM indentities WHERE $code = ?";
    my $sql_handle  = $dbh->prepare($sql);
    $sql_handle->execute($id);
    return $sql_handle->rows();
}

# Check a ric and bbg code are consistent for one id
sub check_id_consistency
{
    my ( $dbh, $ric, $bbg ) = @_;
    my $sql         = 'SELECT id FROM indentities WHERE ric = ? AND bbg = ?';
    my $sql_handle  = $dbh->prepare($sql);
    $sql_handle->execute($ric, $bbg);
    return $sql_handle->rows();
}

sub create_new_id_for_codes
{
    my ( $dbh, $ric, $bbg ) = @_;
    my $sql         = 'INSERT INTO indentities (ric, bbg) VALUES (?, ?)';
    my $sql_handle  = $dbh->prepare($sql);
    $sql_handle->execute($ric, $bbg);
}

sub get_id_from_ric_or_bbg
{
    my ( $dbh, $ric, $bbg ) = @_;
    my $sql         = 'SELECT id FROM indentities WHERE ric = ? OR bbg = ?';
    my $sql_handle  = $dbh->prepare($sql);
    $sql_handle->execute($ric, $bbg);
    
    # Check
    if ( $sql_handle->rows() == 0 )
    {
        print "Error: No mapping found.\n";
        exit 1;
    }
    if ( $sql_handle->rows() > 1 )
    {
        print "Error: Code mapped to multiple ids.\n";
        exit 1;
    }
    
    my @data = $sql_handle->fetchrow_array();
    return $data[0];
}

sub get_ric_from_id
{
    my ( $dbh, $id ) = @_;
    my $sql         = 'SELECT ric FROM indentities WHERE id = ?';
    my $sql_handle  = $dbh->prepare($sql);
    $sql_handle->execute($id);
    
    # Check
    if ( $sql_handle->rows() == 0 )
    {
        print "Error: No mapping found.\n";
        exit 1;
    }
    if ( $sql_handle->rows() > 1 )
    {
        print "Error: Code mapped to multiple ids.\n";
        exit 1;
    }
    
    my @data = $sql_handle->fetchrow_array();
    return $data[0];
}

sub get_bbg_from_id
{
    my ( $dbh, $id ) = @_;
    my $sql         = 'SELECT bbg FROM indentities WHERE id = ?';
    my $sql_handle  = $dbh->prepare($sql);
    $sql_handle->execute($id);
    
    # Check
    if ( $sql_handle->rows() == 0 )
    {
        print "Error: No mapping found.\n";
        exit 1;
    }
    if ( $sql_handle->rows() > 1 )
    {
        print "Error: Code mapped to multiple ids.\n";
        exit 1;
    }
    
    my @data = $sql_handle->fetchrow_array();
    return $data[0];
}

# Removes leading and trailing white space
sub clean_string
{
    my ( $str ) = @_;
    $str =~ s/^\s+//; # Remove leading spaces
    $str =~ s/\s+$//; # Remove trailing spaces
    return $str;
}
