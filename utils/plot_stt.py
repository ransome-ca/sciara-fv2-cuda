#!/bin/python
#                                                                      
# GPL3 License                                                                         
#                                                                      
#                                                                      
# Copyright (c) 2022 Ransome CA                              
#                                                                      
# This file is part of SCIARA-fv2-CUDA.                                          
#                                                                      
# SCIARA-fv2-CUDA is free software: you can redistribute it and/or modify        
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or    
# (at your option) any later version.                                  
#                                                                      
# SCIARA-fv2-CUDA is distributed in the hope that it will be useful,             
# but WITHOUT ANY WARRANTY; without even the implied warranty of       
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        
# GNU General Public License for more details.                         
#                                                                      
# You should have received a copy of the GNU General Public License    
# along with SCIARA-fv2-CUDA.  If not, see <http://www.gnu.org/licenses/>.       
#     
#

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

plt.imshow(fd, cmap='terrain', interpolation='nearest')
plt.colorbar()
plt.show()
