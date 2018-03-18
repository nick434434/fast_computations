set term png
set datafile separator ','
set key autotitle columnhead
set output "graphs_1.png"
set key outside
plot for [col=2:6] "10e7_5thrds" using 1:col with lines
