'''
Calculates the kappa squared of the the local file 1MBD
imports the coordinates from the txt file

calculates the kappa squared for both the normal and disordered transition dipoles 
    estimate: uses only the coordinates of the txt file
    exact: generates coordinates of the file for specific angles
'''

import numpy
import scipy
from scipy import *
from numpy import *
from scipy.optimize import fsolve

#need to change the directory to the file you want to open    
datafile = open('1MBD.pdb', 'r')
data = []


#I think this gets rid of all spaces in "data" variable?
# maybe it generates 2D array?
for row in datafile:
    data.append(row.strip().split())

