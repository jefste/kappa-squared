'''
Calculates the kappa squared of the the local file 1MBD
imports the coordinates from the PDB file

calculates the kappa squared for both the normal and disordered transition dipoles 
    estimate: uses only the coordinates of the PDB file
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

# create a blank dictionary to write to and read from for XYZ coordinates
atomXYZ={}
#lists of tryptophan and heme coordinates to write into atomXYZ dictionary
trpList=['CD1','CD2','CE2','NE1']
hemeList=['FE','CHC','CHD','C3C','C2B']

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


#gets the XYZ coordinates for relevant tryptophan atoms and writes them to a dictionary
def getTrpCoordfromPDB(trpNumber):
    #searches entire rows of file
    for row in data:
        if row[3]=="TRP" and str(row[5])==str(trpNumber):
            #checks to see if atom type is TRP as well as the trpNumber (for example 7)
            for name in trpList:
                #checks to see if the name of the atom is in the trpList
                if row[2]==name:
                    #writes XYZ coordinates to atomXYZ dictionary
                    atomXYZ[name]=getXYZfromRow(row)
    #writes trpOrigin to atomXYZ dictionary
    atomXYZ['trpOrigin']= numpy.add(atomXYZ["CD2"],atomXYZ["CE2"])/2

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


#this grabs the heme coordinates
def getHemeCoordfromPDB():
    #searches entire rows of file
    for i in range(0,len(data)):
        #checks to see if atom type is HETATM
        if data[i][0]=="HETATM":
            for name in hemeList:
                #checks to see if the name of the atom is in the hemeList
                if data[i][2]==name:
                    #writes XYZ coordinates to atomXYZ dictionary
                    atomXYZ[name]=[float(data[i][6]),float(data[i][7]),float(data[i][8])]



#returns the kappa squared value given the angles between the dipoles
def kappaSquared(angle_DA,angle_DT,angle_AT):
    return numpy.square((numpy.cos(angle_DA)-3*numpy.cos(angle_DT)*numpy.cos(angle_AT)))

#gives the angle between 2 vectors, in radians
def dotProductAngle(vector1,vector2):
    numerator1=numpy.dot(vector1,vector2)
    denominator1=numpy.power(numpy.dot(vector1,vector1)*numpy.dot(vector2,vector2),0.5)
    return numpy.arccos(numerator1/denominator1)


#estimates kappa squared without doing any calcuations, only uses coordinates of file
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
    
# used for debugging output of estimated k^2
#print 'normal: %.3f, disordered:, %.3f' % estimationOfKappaSquared()

# generates the parameter 'x' to be used in the function hemeDipoleCoordinates(x,orientation)
# an initial guess for 'x' is needed, typically use .5 as the inital guess. The program works such that it creates a line
# using the points CHC and CHD (for the normal case) for which a coordinate generated on the this line is will 
# be a dipole for the desired angle. By using 0.5, the midpoint of the line is used as a starting guesss.
# Typically an angle of 55 degrees is used for the generated heme dipole, but if it is desired to generate a 'fish' plot
# flexibility for this coordinate will be desireable. 

def generateHemeDipoleParameter(x,angle,orientation):
    
    if orientation=="normal":
        return dotProductAngle(numpy.subtract(hemeCHC,hemeFe),numpy.subtract((x*hemeCHD+(1-x)*hemeCHC),hemeFe))-angle*numpy.pi/180
    if orientation=="disordered":
        return dotProductAngle(numpy.subtract(hemeCHC,hemeFe),numpy.subtract((x*hemeC2B+(1-x)*hemeCHC),hemeFe))-angle*numpy.pi/180

# generates heme diople coordinates. as stated in the comments above, the 'x' parameter needs to be found using the 
# generateHemeDipoleParameter for a specific angle and orientation.
# this function returns a list of [X,Y,Z] coordinates for the generated heme dipole.
def hemeDipoleCoordinates(x,orientation):
    if orientation=="normal":
        return x*hemeCHD+(1-x)*hemeCHC
    if orientation=="disordered":
        return x*hemeC2B+(1-x)*hemeCHC

# works in a very similar manner to generateHemeDipoleParameter function, but don't need flexibility of normal and disordered orientation
# realistically, don't need the angle parameter, because it should always be 38 degrees.
def generateTrpDipoleParameter(x,angle):
    return dotProductAngle(numpy.subtract(trpCD1,trpOrigin),numpy.subtract((x*trpCD1+(1-x)*trpNE1),trpOrigin))-angle*numpy.pi/180

# works in similar manner to hemeDipoleCoordinates.  Needs an 'x' parameter that should be generated by generateTrpDipoleParameter fuction.
def trpDipoleCoordinates(x):
    return x*trpCD1+(1-x)*trpNE1


def kappaSquaredRoutine():
    #this is the transition vector from the trp to the heme
    Trans_V=numpy.subtract(hemeFe,trpOrigin)
    #use this as a rough estimate of the heme normal and disorder dipoles
    #shorten function names just to make it a bit 'easier' to read the lines following code
    hDC= hemeDipoleCoordinates
    gHDP=generateHemeDipoleParameter
    dPA = dotProductAngle
    hemeD_Norm =numpy.subtract(hDC(fsolve(gHDP,.5,(55,"normal")),"normal"),hemeFe)
    hemeD_Dis =numpy.subtract(hDC(fsolve(gHDP,.5,(55,"disordered")),"disordered"),hemeFe)
    #use this as a rough estimate of tryptophan dipole
    trpD = numpy.subtract(trpDipoleCoordinates(fsolve(generateTrpDipoleParameter,.5,38)),trpOrigin)
    #calculates the normal and disordered
    normKap=kappaSquared(dPA(trpD,hemeD_Norm),dPA(trpD,Trans_V),dPA(hemeD_Norm,Trans_V))
    disorderedKap=kappaSquared(dPA(trpD,hemeD_Dis),dPA(trpD,Trans_V),dPA(hemeD_Dis,Trans_V))
    return normKap,disorderedKap


#prints output for user to see difference between estimated and generated kappa squared values    
print "kappa squared for generated coordinates"
print "normal: %.3f \t  disordered: %.3f" % kappaSquaredRoutine()
print 'kappa squared for estimated coordinates'
print 'normal: %.3f \t disordered: %.3f' % estimationOfKappaSquared()

