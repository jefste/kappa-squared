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


#gets the XYZ coordinates for relevant tryptophan atoms and writes them to a dictionary
def getTrpCoordfromPDB(trpNumber):
    #searches entire rows of file
    for row in data: 
        if row[0]=='ATOM':
            if row[3]=='TRP' and str(row[5])==str(trpNumber):
                #checks to see if atom type is TRP as well as the trpNumber (for example 7)
                for name in trpList:
                    #checks to see if the name of the atom is in the trpList
                    if row[2]==name:
                        #writes XYZ coordinates to atomXYZ dictionary
                        atomXYZ[name]=getXYZfromRow(row)
        #writes trpOrigin to atomXYZ dictionary
    atomXYZ['trpOrigin']= numpy.add(atomXYZ['CD2'],atomXYZ['CE2'])/2


#this grabs the heme coordinates
def getHemeCoordfromPDB():
    #searches entire rows of file
    for row in data:
        #checks to see if atom type is HETATM
        if row[0]=="HETATM":
            for name in hemeList:
                #checks to see if the name of the atom is in the hemeList
                if row[2]==name:
                    #writes XYZ coordinates to atomXYZ dictionary
                    atomXYZ[name]=getXYZfromRow(row)



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
    Trans_V=numpy.subtract(atomXYZ['FE'],atomXYZ['trpOrigin'])
    #use this as a rough estimate of the heme normal and disorder dipoles
    hemeD_Norm_Est = numpy.subtract(atomXYZ['C3C'],atomXYZ['FE'])
    hemeD_Dis_Est = numpy.subtract(atomXYZ['C2B'],atomXYZ['FE'])
    #use this as a rough estimate of tryptophan dipole
    trpD_Est = numpy.subtract(atomXYZ['NE1'],atomXYZ['trpOrigin'])
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
    #add to make code a bit easier to read?
    aXYZ=atomXYZ
    
    if orientation=="normal":
        return dotProductAngle(numpy.subtract(aXYZ['CHC'],aXYZ['FE']),numpy.subtract((x*aXYZ['CHD']+(1-x)*aXYZ['CHC']),aXYZ['FE']))-angle*numpy.pi/180
    if orientation=="disordered":
        return dotProductAngle(numpy.subtract(aXYZ['CHC'],aXYZ['FE']),numpy.subtract((x*aXYZ['C2B']+(1-x)*aXYZ['CHC']),aXYZ['FE']))-angle*numpy.pi/180

# generates heme diople coordinates. as stated in the comments above, the 'x' parameter needs to be found using the 
# generateHemeDipoleParameter for a specific angle and orientation.
# this function returns a list of [X,Y,Z] coordinates for the generated heme dipole.
def hemeDipoleCoordinates(x,orientation):
    if orientation=="normal":
        return x*atomXYZ['CHD']+(1-x)*atomXYZ['CHC']
    if orientation=="disordered":
        return x*atomXYZ['C2B']+(1-x)*atomXYZ['CHC']

# works in a very similar manner to generateHemeDipoleParameter function, but don't need flexibility of normal and disordered orientation
# realistically, don't need the angle parameter, because it should always be 38 degrees.
def generateTrpDipoleParameter(x,angle):
    return dotProductAngle(numpy.subtract(atomXYZ['CD1'],atomXYZ['trpOrigin']),numpy.subtract((x*atomXYZ['CD1']+(1-x)*atomXYZ['NE1']),atomXYZ['trpOrigin']))-angle*numpy.pi/180

# works in similar manner to hemeDipoleCoordinates.  Needs an 'x' parameter that should be generated by generateTrpDipoleParameter fuction.
def trpDipoleCoordinates(x):
    return x*atomXYZ['CD1']+(1-x)*atomXYZ['NE1']


def distanceCentertoCenter(trpNumber):
    getTrpCoordfromPDB(trpNumber)
    getHemeCoordfromPDB()
    return numpy.linalg.norm(subtract(atomXYZ["FE"],atomXYZ["trpOrigin"]))

def kappaSquaredRoutine(trpNumber):
    #shorten function names just to make it a bit 'easier' to read the lines following code
    hDC= hemeDipoleCoordinates
    gHDP=generateHemeDipoleParameter
    dPA = dotProductAngle
    aXYZ = atomXYZ
    #grabs coordinates for trp(specify number) and heme
    getTrpCoordfromPDB(trpNumber)
    getHemeCoordfromPDB()
    #this is the transition vector from the trp to the heme
    Trans_V=numpy.subtract(aXYZ['FE'],aXYZ['trpOrigin'])
    #use this as a rough estimate of the heme normal and disorder dipoles
    hemeD_Norm =numpy.subtract(hDC(fsolve(gHDP,.5,(55,"normal")),"normal"),aXYZ['FE'])
    hemeD_Dis =numpy.subtract(hDC(fsolve(gHDP,.5,(55,"disordered")),"disordered"),aXYZ['FE'])
    #use this as a rough estimate of tryptophan dipole
    trpD = numpy.subtract(trpDipoleCoordinates(fsolve(generateTrpDipoleParameter,.5,38)),aXYZ['trpOrigin'])
    #calculates the normal and disordered
    normKap=kappaSquared(dPA(trpD,hemeD_Norm),dPA(trpD,Trans_V),dPA(hemeD_Norm,Trans_V))
    disorderedKap=kappaSquared(dPA(trpD,hemeD_Dis),dPA(trpD,Trans_V),dPA(hemeD_Dis,Trans_V))
    return normKap,disorderedKap


#prints output for user to see difference between estimated and generated kappa squared values    
print 'center to center distance between dipoles (in Angstroms):',distanceCentertoCenter()
print "kappa squared for generated coordinates"
print "normal: %.3f \t disordered: %.3f" % kappaSquaredRoutine(7)
print 'kappa squared for estimated coordinates'
print 'normal: %.3f \t disordered: %.3f' % estimationOfKappaSquared()

