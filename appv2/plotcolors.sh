#!/usr/bin/gnuplot --persist
RGB(R,G,B) =  int(255.*R) * 2**16 + int(255.*G) * 2**8  + int(255.*B)

set isosamples 1000,1000
set xyplane at 0.0
set view map
#set cbrange [0:2]

#set xrange [-1:1]
#set yrange [-1:1]
#set zrange [-1:1]
set title "RGB coloring of pm3d surface"
unset key
do for[i=1:11]{
	set term pngcairo
	outfile = sprintf("chapter_6_%d", i)
	set output outfile

	infile = sprintf("FCT_%d.dat", i)
	stat infile using 3 prefix 'R' nooutput
	stat infile using 4 prefix 'G' nooutput
	stat infile using 5 prefix 'B' nooutput
	splot infile using 1:2:(RGB($3/R_max,$4/G_max,$5/B_max)) \
			   with pm3d lc rgb variable
}

