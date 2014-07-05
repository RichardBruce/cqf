set pm3d
set ticslevel 0
splot "./fv_datafile.dat" u 2:1:3 every 20:1 with lines

set zrange [-0.1:0.1]
set autoscale
splot "./delta_datafile.dat" u 2:1:3 every 50:1 with lines

splot "./datafile.dat" u 2:1:3 every 1:::0::0 with lines
splot "./datafile.dat" u 2:1:3 every 1::950::1050 with lines

splot "output" u 1:2:3 every 5:5 with lines
