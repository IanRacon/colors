# ------------------------------------------------------------------
# Rysowanie gifow
# ------------------------------------------------------------------

reset
set term gif animate size 800,300
set output "animation.gif"
n=50    #n frames
set pm3d
set view map
set size ratio -1
set cbrange [0:1]

i=3
load "animation.plt"
set output

