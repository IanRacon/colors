splot "density.dat" u 1:2:i w pm3d  title sprintf("t=%i",i)
i=i+1
if (i < n) reread
