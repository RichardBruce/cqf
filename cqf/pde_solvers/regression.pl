#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $help            = 0;
my $verbose         = 0;

my @exclude_solver;

my %adi_2d_params = (
    low_p =>
    {
        nas1 => 600,
        nas2 => 600,
        nts  => 100
    },
    med_p =>
    {
        nas1 => 900,
        nas2 => 900,
        nts  => 100
    },
    hi_p =>
    {
        nas1 => 1200,
        nas2 => 1200,
        nts  => 100
    }
);

my %adi_2d_tests = (
    wo_call =>
    {
        cashflows => [ [ "wo_call", 100.0, 100.0 ] ],
        economics =>
        {
            s1 => [ 100.0, 0.05, 0.2 ],
            s2 => [ 100.0, 0.05, 0.2 ]
        },
        exact   => 6.0
    },
    bo_call =>
    {
        cashflows => [ [ "bo_call", 100.0, 100.0 ] ],
        economics =>
        {
            s1 => [ 100.0, 0.05, 0.2 ],
            s2 => [ 100.0, 0.05, 0.2 ]
        },
        exact   => 6.0
    },
    wo_binary_call =>
    {
        cashflows => [ [ "wo_binary_call", 100.0, 100.0 ] ],
        economics =>
        {
            s1 => [ 100.0, 0.05, 0.2 ],
            s2 => [ 100.0, 0.05, 0.2 ]
        },
        exact   => 6.0
    }
);

my %pde_1d_params = (
    param_s06 =>
    {
        nas1 => 600,
        nts  => 500
    },
    param_s09 =>
    {
        nas1 => 900,
        nts  => 500
    },
    param_s12 =>
    {
        nas1 => 1200,
        nts  => 500
    },
    param_s15 =>
    {
        nas1 => 1500,
        nts  => 500
    },
    param_s18 =>
    {
        nas1 => 1800,
        nts  => 500
    },
    param_s21 =>
    {
        nas1 => 2100,
        nts  => 500
    },
    param_s24 =>
    {
        nas1 => 2400,
        nts  => 500
    },
    param_s27 =>
    {
        nas1 => 2700,
        nts  => 500
    },
    param_s30 =>
    {
        nas1 => 3000,
        nts  => 500
    },
    param_t05 =>
    {
        nas1 => 3000,
        nts  => 50
    },
    param_t1 =>
    {
        nas1 => 3000,
        nts  => 100
    },
    param_t2 =>
    {
        nas1 => 3000,
        nts  => 200
    },
    param_t3 =>
    {
        nas1 => 3000,
        nts  => 300
    },
    param_t4 =>
    {
        nas1 => 3000,
        nts  => 400
    },
    param_t5 =>
    {
        nas1 => 3000,
        nts  => 500
    },
    param_t6 =>
    {
        nas1 => 3000,
        nts  => 600
    },
    param_t7 =>
    {
        nas1 => 3000,
        nts  => 700
    },
    param_t8 =>
    {
        nas1 => 3000,
        nts  => 800
    },
    param_t9 =>
    {
        nas1 => 3000,
        nts  => 900
    }
);

my %grid_params = (
    uniform =>
    {
        nas1    => 90,
        nts     => 100,
        grid    => "uniform 1.0 1.0"        
    },
    sinh_1 =>
    {
        nas1    => 90,
        nts     => 100,
        grid    => "sinh 100.0 1.0"        
    },
    sinh_5 =>
    {
        nas1    => 90,
        nts     => 100,
        grid    => "sinh 100.0 5.0"        
    },
    sinh_10 =>
    {
        nas1    => 90,
        nts     => 100,
        grid    => "sinh 100.0 10.0"        
    },
    sinh_15 =>
    {
        nas1    => 90,
        nts     => 100,
        grid    => "sinh 100.0 15.0"        
    },
    sinh_20 =>
    {
        nas1    => 90,
        nts     => 100,
        grid    => "sinh 100.0 20.0"        
    },
    sinh_25 =>
    {
        nas1    => 90,
        nts     => 100,
        grid    => "sinh 100.0 25.0"        
    },
    sinh_30 =>
    {
        nas1    => 90,
        nts     => 100,
        grid    => "sinh 100.0 30.0"        
    },
    sinh_35 =>
    {
        nas1    => 90,
        nts     => 100,
        grid    => "sinh 100.0 35.0"        
    },
    sinh_40 =>
    {
        nas1    => 90,
        nts     => 100,
        grid    => "sinh 100.0 40.0"        
    },
    sinh_45 =>
    {
        nas1    => 90,
        nts     => 100,
        grid    => "sinh 100.0 45.0"        
    },
    sinh_50 =>
    {
        nas1    => 90,
        nts     => 100,
        grid    => "sinh 100.0 50.0"        
    }
);

my %pde_1d_tests = (
    call =>
    {
        cashflows => [ [ "call", 100.0 ] ],
        economics =>
        {
            s1 => [ 100.0, 0.05, 0.2 ]
        },
        exact   => 10.4506
    },
    put =>
    {
        cashflows => [ [ "put", 100.0 ] ],
        economics =>
        {
            s1 => [ 100.0, 0.05, 0.2 ]
        },
        exact   => 5.57352
    },
    binary_call =>
    {
        cashflows => [ [ "binary_call", 100.0 ] ],
        economics =>
        {
            s1 => [ 100.0, 0.05, 0.2 ]
        },
        exact   => 0.532325
    },
    binary_put =>
    {
        cashflows => [ [ "binary_put", 100.0 ] ],
        economics =>
        {
            s1 => [ 100.0, 0.05, 0.2 ]
        },
        exact   => 0.418905
    }
);

my %test_cases = (
    explicit =>
    {
        tests => { %pde_1d_tests },
        params => { %pde_1d_params },
    },
    implicit =>
    {
        tests => { %pde_1d_tests },
        params => { %pde_1d_params },
    },
    cn =>
    {
        tests => { %pde_1d_tests },
        #params => { %pde_1d_params },
        params => { %grid_params }
    },
    doug =>
    {
        tests => { %pde_1d_tests },
        params => { %pde_1d_params },
    },
    cs =>
    {
        tests => { %adi_2d_tests },
        params => { %adi_2d_params },
    },
    mcs =>
    {
        tests => { %adi_2d_tests },
        params => { %adi_2d_params },
    },
    hv =>
    {
        tests => { %adi_2d_tests },
        params => { %adi_2d_params },
    }
);


# Argument handling
GetOptions (
    "exclude:s"         => \@exclude_solver,
    "v"                 => \$verbose,
    "help"              => \$help
    );

# Help message
if ($help)
{
    print "Usage: regression [--exclude <solver>] [--v] [--help]\n.";
    exit 0;
}


# Convert excluded solvers to a hash
my %exclude_solvers = map { $_ => 1 } @exclude_solver;

# Call the PDE solver
# Foreach PDE solver under test
my %results;
for my $solver (sort keys %test_cases)
{
    # Fetch the definition of the test cases
    print $solver . "\n";
    next if (exists($exclude_solvers{$solver}));
    
    my %tests = %{ $test_cases{$solver}{tests} };
    my %params = %{ $test_cases{$solver}{params} };
    
    # Foreach test on this solver
    for my $test (sort keys %tests)
    {
        # Build the payoff
        my $cashflow = "";
        for (@{$tests{$test}{cashflows}})
        {
            $cashflow = $cashflow . " -cashflow " . join(' ', @$_);
        }
        
        # Build the economics arguement to the PDE solver
        my $economics = "";
        for (keys %{$tests{$test}{economics}})
        {
            $economics = $economics . " -economics $_ " . join(' ', @{$tests{$test}{economics}{$_}});
        }
            
        # Foreach set of parameters to test
        for my $param (sort keys %params)
        {
            # Build the parameters to the PDE solver
            my $pde_params = "";
            for (keys %{$params{$param}})
            {
                $pde_params = $pde_params . " -$_ $params{$param}{$_}";
            }
           
            # Call the pde solver
            my $call = "./pde_solver -solver $solver $pde_params $cashflow $economics";
            print "$call\n" unless !$verbose;
            my @result = `$call`;
            for (@result)
            {
                if (m/result: (-?[\d\.\-e]+)/i)
                {
                    $results{$solver}{$test}{$param}{result} = $1;
                    my $exact = $test_cases{$solver}{tests}{$test}{exact};
                    $results{$solver}{$test}{$param}{error} = (($1 - $exact) / $exact) * 100.0;
                }
            }
        }
    }
}

# Dump as csv
print "Solver,Test,Nas1,Nas2,Nts,Value,Error\n";
for my $solver (sort keys %results)
{
    for my $test (sort keys %{$results{$solver}})
    {
        for my $param (sort keys %{$results{$solver}{$test}})
        {
            my $nas1 = $test_cases{$solver}{params}{$param}{nas1};
            my $nas2 = $test_cases{$solver}{params}{$param}{nas2};
            my $nts  = $test_cases{$solver}{params}{$param}{nts};
            my $value = $results{$solver}{$test}{$param}{result};
            my $error = $results{$solver}{$test}{$param}{error};
            print "$solver,$test,$nas1,$nas2,$nts,$value,$error\n";
        }
    }
}
