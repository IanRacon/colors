#!/usr/bin/gnuplot
reset
n=100 #n frames
set term gif animate #size 1024,1024
set output "animation.gif"
set pm3d
set view map
set size ratio -1
set cbrange [0:2]
i=3
load "animation.plt"
set output
