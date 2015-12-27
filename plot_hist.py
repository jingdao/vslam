#!/usr/bin/python

import sys
import numpy
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
	print "./plot_hist map_point.txt"
	sys.exit(1)

f = open(sys.argv[1],"r")
RMSE = []
for l in f:
	RMSE.append(float(l.split()[4]))

plt.hist(numpy.log(RMSE),bins=10,normed=True)
plt.xlabel("log RMSE")
plt.ylabel("Probability")
plt.show()
