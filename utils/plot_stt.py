#!/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt

fd = open(sys.argv[1], 'r')
fd = fd.readlines()
fd = fd[6:]
fd = [x.strip() for x in fd]
fd = [x.split() for x in fd]
fd = [[float(y) for y in x] for x in fd]
fd = np.array(fd)


plt.imshow(fd, cmap='hot', interpolation='nearest')
plt.colorbar()
plt.show()
