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


# grabs XYZ coordinates from a given row and returns them as a list of 3 dimensions
def getXYZfromRow(row):
    return [float(row[6]),float(row[7]),float(row[8])]


#this grabs the tryptophan coordinates, only for the specific case. 
#should change this so that the range is more flexible, and also so that it prompts for which tryptophan is of interest
for row in data[372:386]:
    if row[2]=='CD1':
        trpCD1=getXYZfromRow(row)
    if row[2]=='CD2':
        trpCD2=getXYZfromRow(row)
    if row[2]=='CE2':
        trpCE2=getXYZfromRow(row)
    if row[2]=='NE1':
        trpNE1=getXYZfromRow(row)
