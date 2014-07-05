/usr/bin/time ./pde_solver -economics 50 0.05 0.2 -cashflow call asian_in 1 0.0,0.1,0.2,0.3,1 . . . . -solver -grid uniform 3.0 0.0 0.5 uniform 3.0 0.0 0.5 -risk dspot1 spot.log grid 0.5 . . . -risk value value.log grid 0.5 . . .
/usr/bin/time ./pde_solver -economics 50 0.05 0.2 -cashflow call asian 50 0.94,0.96,0.98,1 . . . . -solver -grid uniform 3.0 0.0 0.5 uniform 3.0 0.0 0.5 -risk dspot1 spot.log grid 0.0 . . . -risk value value.log grid . . . .
valgrind ./pde_solver -cashflow put asian 100 9.0,9.5,10 . . . . -solver -grid uniform 3.0 0.0 0.5 uniform 3.0 0.0 0.5 -economics 100 0.05 0.2 -risk dspot1 spot.log grid . . . .

set pm3d
set ticslevel 0
splot "./fv_datafile.dat" u 2:1:3 every 20:1 with lines

set zrange [-0.1:0.1]
set autoscale
splot "./delta_datafile.dat" u 2:1:3 every 50:1 with lines

splot "./datafile.dat" u 2:1:3 every 1:::0::0 with lines
splot "./datafile.dat" u 2:1:3 every 1::950::1050 with lines

splot "output" u 1:2:3 every 5:5 with lines
