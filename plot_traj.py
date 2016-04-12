#!/usr/bin/python

import sys
from matplotlib import pyplot as plt
if len(sys.argv) < 3:
	print './plot_traj src.txt dst.txt'
	sys.exit(1)

x1=[]
y1=[]
x2=[]
y2=[]
f = open(sys.argv[1])
for l in f:
	t = l.split()
	x1.append(float(t[1]))
	y1.append(float(t[2]))
f.close()
f = open(sys.argv[2])
for l in f:
	t = l.split()
	x2.append(float(t[1]))
	y2.append(float(t[2]))
f.close()

plt.plot(x1,y1,label=sys.argv[1].split('.')[0])
plt.plot(x2,y2,label=sys.argv[2].split('.')[0])
plt.legend()
plt.show()
