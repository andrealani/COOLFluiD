#!/bin/bash
#
# first parameter:  update time of the pdf file
# second parameter: number of last steps

# default values
UPDATETIME=10;
NUMSTEPS=100




#### no changes below ####

if [ $# -eq 1 ]; then
	UPDATETIME=$1;
fi
if [ $# -eq 2 ]; then
	UPDATETIME=$1;
	NUMSTEPS=$2
fi

if [ ! -f convergence-P0.plt ]; then
	echo "convergence file not found"
	exit 1
fi

while [ 1 ]; do

	gnuplot << EOF
#set terminal png
set terminal postscript eps color
set output "convergence.eps"
set datafile separator ' '
set multiplot           # engage multiplot mode
set size 1.0,1.0 
set grid
plot "convergence-P0.plt" every ::2 using 1:2  notitle with linespoints
set object 1 rectangle from screen 0.0,0.0 to screen 1.0,1.0 fillcolor rgb"#ffffff" behind
#set style line 1 lt 2 lc rgb "red" lw 1 pt 1 ps 2
#plot "convergence-P0.plt" every ::2 using 1:2  notitle with linespoints
#set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]
plot "convergence-P0.plt" every ::2 using 1:2  notitle with linespoints
set size 0.4,0.4        # set the size of the second plot in plot units
set origin 0.55,0.55     # set the origin for the second plot in plot units
set xrange [GPVAL_DATA_X_MAX-$NUMSTEPS:GPVAL_DATA_X_MAX]
set object 1 rectangle from screen 0.58,0.55 to screen 0.95,0.95 fillcolor rgb"#ffffff" behind
plot "convergence-P0.plt" every ::2 using 1:2  notitle with linespoints
EOF
	epstopdf convergence.eps
	rm convergence.eps
	sleep $UPDATETIME
done
