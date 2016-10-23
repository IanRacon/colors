#!/usr/bin/env python3
it = 0
tab = []
tab2 = []
for i in open("result.dat"):
    tab.append(i.split()[0])
    it+=1
    if not it%10:
        tab2.append(tab)
        tab = []
with open("formatted.dat", "w") as f:
    for i in tab2:
        f.write("".join(i) + "\n")
    
