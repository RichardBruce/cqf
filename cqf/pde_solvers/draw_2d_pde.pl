#!/usr/bin/env perl

use strict;
use Getopt::Long;

my $help            = 0;
my $verbose         = 0;
my $no_pde          = 0;
my $no_plot         = 0;
my $no_video        = 0;
my $nts             = 100;
my $nas1            = 600;
my $nas2            = 600;
my $solver          = "doug";
my $payoff          = "wo_binary_call";
my $k_s1            = 100.0;
my $k_s2            = 100.0;
my $t               = 1.0;
my $pde_log         = "./logs/pde_log_";
my $plot_log        = "./logs/pde_plot_";
my $plot_max_value  = 300.0;
my $video_fps       = 30;
my $video_x_res     = 1280;
my $video_y_res     = 720;
my $video_file      = "./logs/pde_video.avi";


# Argument handling
GetOptions (
    "nts:i"             => \$nts,
    "nas1:i"            => \$nas1,
    "nas2:i"            => \$nas2,
    "solver:s"          => \$solver,
    "payoff:s"          => \$payoff,
    "k_s1:f"            => \$k_s1,
    "k_s2:f"            => \$k_s2,
    "t:f"               => \$t,
    "pde_log:s"         => \$pde_log,
    "plot_log:s"        => \$plot_log,
    "plot_max_value:f"  => \$plot_max_value,
    "video_fps:i"       => \$video_fps,
    "video_file:s"      => \$video_file,
    "no_pde"            => \$no_pde,
    "no_plot"           => \$no_plot,
    "no_video"          => \$no_video,
    "v"                 => \$verbose,
    "help"              => \$help
    );

#if ($nas1 == $nas2)
#{
#    print "Error: Must have different number of steps in s1 and s2\n";
#}

    
# Call the PDE solver
if (!$no_pde)
{
    `./pde_solver -nts $nts -nas1 $nas1 -nas2 $nas2 -log_file $pde_log -t $t -cashflow $payoff $k_s1 $k_s2 -solver $solver`;
}

# Plot 3D graphs
if (!$no_plot)
{
    my $i;
    for ($i = 0; $i < $nts; $i++)
    {
        my $time = $i * ($t / $nts);
        my $title = "Time: " . sprintf("%.3f", $time);
        my $plot_cmd = "gnuplot -e \"set ticslevel 0;";
        $plot_cmd = $plot_cmd . "set view 60,310;";
        $plot_cmd = $plot_cmd . "set zrange [0:$plot_max_value];";
        $plot_cmd = $plot_cmd . "set xlabel 'S1';";
        $plot_cmd = $plot_cmd . "set ylabel 'S2';";
        $plot_cmd = $plot_cmd . "set zlabel 'V';";
        $plot_cmd = $plot_cmd . "set terminal jpeg enhanced size $video_x_res,$video_y_res;";
        $plot_cmd = $plot_cmd . "set output '$plot_log$i.jpg';";
        $plot_cmd = $plot_cmd . "set title '$title';";
        $plot_cmd = $plot_cmd . "unset key;";
        $plot_cmd = $plot_cmd . "splot '$pde_log$i' u 1:2:3 every 25:25 with lines\"";
        `$plot_cmd`;
    }
}

# Compress graphs into video
if (!$no_video)
{
    `../../lib/ffmpeg/ffmpeg -g $video_fps -flags +bitexact -flags2 +wpred -qscale 1 -i $plot_log%d.jpg -y $video_file`;
}
