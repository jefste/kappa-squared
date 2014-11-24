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


#calculates the trpOrigin, needed for calculating the tryptophan coordinates
#does this need to be global? may want to change to local if dealing with more tryptophans??
trpOrigin = numpy.add(trpCD2,trpCE2)/2

#this grabs the heme coordinates
for row in data:
    if row[0]=='HETATM':
        if row[2]=='FE':
            hemeFe=getXYZfromRow(row)
        if row[2]=='CHC':
            hemeCHC=getXYZfromRow(row)
        if row[2]=='CHD':
            hemeCHD=getXYZfromRow(row)
        if row[2]=='C3C':
            hemeC3C=getXYZfromRow(row)
        if row[2]=='C2B':
            hemeC2B=getXYZfromRow(row)

#returns the kappa squared value given the angles between the dipoles
def kappaSquared(angle_DA,angle_DT,angle_AT):
    return numpy.square((numpy.cos(angle_DA)-3*numpy.cos(angle_DT)*numpy.cos(angle_AT)))

#gives the angle between 2 vectors, in radians
def dotProductAngle(vector1,vector2):
    numerator1=numpy.dot(vector1,vector2)
    denominator1=numpy.power(numpy.dot(vector1,vector1)*numpy.dot(vector2,vector2),0.5)
    return numpy.arccos(numerator1/denominator1)

def estimationOfKappaSquared():
    #shorten the function name to make code more readable
    dPA = dotProductAngle
    #this is the transition vector from the trp to the heme
    Trans_V=numpy.subtract(hemeFe,trpOrigin)
    #use this as a rough estimate of the heme normal and disorder dipoles
    hemeD_Norm_Est = numpy.subtract(hemeC3C,hemeFe)
    hemeD_Dis_Est = numpy.subtract(hemeC2B,hemeFe)
    #use this as a rough estimate of tryptophan dipole
    trpD_Est = numpy.subtract(trpNE1,trpOrigin)
    kappa_est_normal=kappaSquared(dPA(trpD_Est,hemeD_Norm_Est),dPA(trpD_Est,Trans_V),dPA(hemeD_Norm_Est,Trans_V))
    kappa_est_dis=kappaSquared(dPA(trpD_Est,hemeD_Dis_Est),dPA(trpD_Est,Trans_V),dPA(hemeD_Dis_Est,Trans_V))
    return kappa_est_normal,kappa_est_dis
    
print 'normal: %.3f, disordered:, %.3f' % estimationOfKappaSquared()
