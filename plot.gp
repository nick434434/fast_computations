set term png
set datafile separator ','
set key autotitle columnhead
set output "monte_time_dim3_10.png"
set key outside
set title "Npoints = 1000000"
plot for [col=2:3] "out_monte_time.csv" using 1:col with lines

set output "monte_time_dim10.png"
set title "Npoints = 100000000"
plot for [col=4:4] "out_monte_time.csv" using 1:col with lines
