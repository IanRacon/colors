#!/usr/bin/gnuplot
set term png
set pm3d map
set view map
set cbrange [0:2]
#filenames = system "ls *dat"
do for [i=1:10]{
	outfile = sprintf("chapter_1_2_%d_1", i)
	set output outfile
	infile = sprintf("bicgstab_2_%d.dat", i)
	splot infile with pm3d
}


