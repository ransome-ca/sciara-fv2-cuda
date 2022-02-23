#!/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

fd = open(sys.argv[1], 'r')
fd = fd.readlines()
fd = fd[6:]
fd = [x.strip() for x in fd]
fd = [x.split() for x in fd]
fd = [[float(y) for y in x] for x in fd]
fd[-1][-1] = 0.0
fd = np.array(fd)

# if sys.argv[1] == 'data/2006/output_2006_000000016000_Temperature.stt':

#     unique, counts = np.unique(fd, return_counts=True)

#     print("Unique values count in Temperature:")
#     for i in range(len(unique)):
#         print(f'Value {unique[i]} has count {counts[i]}')
    
#     plt.imshow(fd, cmap='autumn_r', interpolation='nearest')

# else:

plt.imshow(fd, cmap='hot_r', interpolation='nearest')
plt.colorbar()
plt.show()
