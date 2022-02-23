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

fd1 = open(sys.argv[1], 'r')
fd1 = fd.readlines()
fd1 = fd[6:]
fd1 = [x.strip() for x in fd]
fd1 = [x.split() for x in fd]
fd1 = [[float(y) for y in x] for x in fd]
fd1 = np.array(fd)

fd2 = open(sys.argv[2], 'r')
fd2 = fd.readlines()
fd2 = fd[6:]
fd2 = [x.strip() for x in fd]
fd2 = [x.split() for x in fd]
fd2 = [[float(y) for y in x] for x in fd]
fd2 = np.array(fd)


diff = np.subtract(fd1, fd2)
indices = np.where(diff != 0)

print(f'Difference between {sys.argv[1]} and {sys.argv[2]}:')

for i in range(len(indices[0])):
    print(f'{indices[0][i]}, {indices[1][i]}: {'cuda' if diff[indices[0][i]][indices[1][i]] else 'serial'}')
