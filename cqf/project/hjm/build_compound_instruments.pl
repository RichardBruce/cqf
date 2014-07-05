#!/usr/bin/env perl

use strict;
use Getopt::Long;


my $help                = 0;
my $verbose             = 0;

my $build_coupon_bond   = 0;
my $build_cap           = 0;
my $build_floor         = 0;
my $build_cap_collar    = 0;
my $build_floor_collar  = 0;
my $build_option        = 0;

my $cur_date        = 0;
my $period          = 1/12;
my $end_date        = 1;
my $cap_strike      = 0.035;
my $floor_strike    = 0.035;
my $rate            = 1;
my $floor_barrier   = 0.0;
my $cap_barrier     = 0.0;
my $coupon_rate     = 0.03;
my $option_strike   = 0.0;

my $nr_of_paths = 100000;#0.0770169

# Argument handling
GetOptions (
    "cap"               => \$build_cap,
    "floor"             => \$build_floor,
    "cap_collar"        => \$build_cap_collar,
    "floor_collar"      => \$build_floor_collar,
    "coupon_bond"       => \$build_coupon_bond,
    "option"            => \$build_option,
    "nr_of_paths:i"     => \$nr_of_paths,
    "end_date:i"        => \$end_date,
    "period:f"          => \$period,
    "cap_strike:f"      => \$cap_strike,
    "floor_strike:f"    => \$floor_strike,
    "cap_barrier:f"     => \$cap_barrier,
    "floor_barrier:f"   => \$floor_barrier,
    "coupon_rate:f"     => \$coupon_rate,
    "option_strike:f"   => \$option_strike,
    "rate:f"            => \$rate,
    "v"                 => \$verbose,
    "help"              => \$help
    );

if ($help)
{
    print "Usage: build_compound_instruments [options]\n";
    print "     -cap                Build a cap\n";
    print "     -floor              Build a floor\n";          
    print "     -cap_collar         Build a collar that is long the cap\n";
    print "     -floor_collar       Build a collar that is long the floor\n";
    print "     -coupon_bond        Build a coupon paying bond\n";
    print "     -option             Build an option on the other payoffs\n";
    print "     -end_date:i         Time to maturity (years)\n";
    print "     -period:i           Length of the period between coupons or evaluations (years)\n";
    print "     -cap_strike:f       Strike of the cap\n";
    print "     -floor_strike:f     Strike of the floor\n";
    print "     -cap_barrier:f      Knock out barrier for the cap\n";
    print "     -floor_barrier:f    Knock out barrier for the floor\n";
    print "     -coupon_rate:f      The size of the coupon paid on the coupon bond\n";
    print "     -rate:f             The rate the cap/floor is struck on\n";
    print "     -v                  Run in verbose mode\n";
    print "     -help               Print this message\n";
}

# Sanity checks    
if ($cap_strike > 0.1)
{
    print "Warning: Thats a very high cap strike.\n";
}

if ($floor_strike > 0.1)
{
    print "Warning: Thats a very high floor strike.\n";
}

if ($cap_barrier > 0.1)
{
    print "Warning: Thats a very high cap barrier.\n";
}

if ($floor_barrier > 0.1)
{
    print "Warning: Thats a very high floor barrier.\n";
}

if ($rate > 25)
{
    print "Warning: Thats a very high rate.\n";
}

if (($build_cap or $build_floor or $build_cap_collar or $build_floor_collar or $build_coupon_bond) == 0)
{
    print "Warning: No options to build.\n";
}


# Build cap and/or floor options
my $cap             = "";
my $floor           = "";
my $cap_collar      = "";
my $floor_collar    = "";
my $coupon_bond     = "";

my $caplet = " -caplet ";
if ($cap_barrier > 0.0)
{
    $caplet = " -ko_caplet ";
}

my $floorlet = " -floorlet ";
if ($floor_barrier > 0.0)
{
    $floorlet = " -ko_floorlet ";
}

while (($cur_date + $period) <= $end_date)
{
    if ($build_coupon_bond)
    {
        $coupon_bond = $coupon_bond . " -zcb 0 " . ($cur_date + $period) . " " . $coupon_rate . " ";
    }
    
    if ($build_cap)
    {
        $cap = $cap . $caplet . " 0 " . ($cur_date + $period) . " " . ($cur_date + $period + $period) . " " . $period . " 1 " . $cap_strike . " " . $rate . " ";
        if ($cap_barrier > 0.0)
        {
            $cap = $cap . " " . $cap_barrier;
        }
    }
    
    if ($build_floor)
    {
        $floor = $floor . $floorlet . " 0 " . ($cur_date + $period) . " " . ($cur_date + $period + $period) . " " . $period . " 1 " . $floor_strike . " " . $rate . " ";
        if ($floor_barrier > 0.0)
        {
            $floor = $floor . " " . $floor_barrier;
        }        
    }
    
    if ($build_cap_collar)
    {
        # Long cap
        $cap_collar = $cap_collar . $caplet . " 0 " . ($cur_date + $period) . " " . ($cur_date + $period + $period) . " " . $period . " 1 " . $cap_strike . " " . $rate . " ";
        if ($cap_barrier > 0.0)
        {
            $cap_collar = $cap_collar . " " . $cap_barrier;
        }

        # Short floor
        $cap_collar = $cap_collar . $floorlet . " 0 " . ($cur_date + $period) . " " . ($cur_date + $period + $period) . " " . $period . " -1 " . $floor_strike . " " . $rate . " ";
        if ($floor_barrier > 0.0)
        {
            $cap_collar = $cap_collar . " " . $floor_barrier;
        }        
    }
    
    if ($build_floor_collar)
    {
        # Short cap
        $floor_collar = $floor_collar . $caplet . " 0 " . ($cur_date + $period) . " " . ($cur_date + $period + $period) . " " . $period . " -1 " . $cap_strike . " " . $rate . " ";
        if ($cap_barrier > 0.0)
        {
            $floor_collar = $floor_collar . " " . $cap_barrier;
        }

        # Long floor
        $floor_collar = $floor_collar . $floorlet . " 0 " . ($cur_date + $period) . " " . ($cur_date + $period + $period) . " " . $period . " 1 " . $floor_strike . " " . $rate . " ";
        if ($floor_barrier > 0.0)
        {
            $floor_collar = $floor_collar . " " . $floor_barrier;
        }        
    }
    
    # Move to next period
    $cur_date += $period;
}

# Print option tags
if ($build_option)
{
    print $option . " " . $option_stike . "\n";
}

# Print payoffs
if ($build_coupon_bond)
{
    $coupon_bond = $coupon_bond . " -zcb 0 " . $end_date . " 1 ";
    print $coupon_bond . "\n";
}

if ($build_cap)
{
    print $cap . "\n";
}

if ($build_floor)
{
    print $floor . "\n";
}

if ($build_cap_collar)
{
    print $cap_collar . "\n";
}

if ($build_floor_collar)
{
    print $floor_collar . "\n";
}

# Print end of option tags
if ($build_option)
{
    print "-end_option\n";
}

exit 0;
