import matplotlib.pyplot
from mpl_toolkits.mplot3d import Axes3D

import sys

import urllib
import requests
import os
import string


import numpy
import scipy
from scipy import *
from numpy import *
from scipy.optimize import fsolve

#use this now for trouble shooting    
datafile = open('PDB/1MBD.pdb', 'r')

data = []
# create a blank dictionary to write to and read from for XYZ coordinates
atomXYZ={}
#lists of tryptophan and heme coordinates to write into atomXYZ dictionary
trpList=['CD1','CD2','CE2','NE1']
hemeList=['FE','CHC','CHD','C3C','C2B']

#I think this gets rid of all spaces in "data" variable?
# maybe it generates 2D array?
def readFromDatafile(pdbFileName):
    datafile = open('PDB/'+pdbFileName+'.pdb', 'r')
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


def grabPDB():
    while True:
        #prompts user for input of PDB structure. PDB is not case sensitive
        # checks to see if an value was given for the argument
        if len(sys.argv)>1:
            # reads 1st argument after k-s.py and converts to uppercase
            pdbID = sys.argv[1].upper()
        else:
            # prompts user for input
            pdbID = str(raw_input("enter the PDB Structure: ")).upper()
        
        url = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId="+pdbID
        pathnamePDB="PDB/"+pdbID+".pdb"
        #checks to see if the path and the file name exists. Is this redunant? Maybe only need to check for file exisiting?
        #Note that when checking exists or isfile, the function does not appear to be case dependant
        if os.path.exists(pathnamePDB) and os.path.isfile(pathnamePDB):
            print  "reading file locally"
            return pdbID.upper()
        else:
            #checks to see if PDB exists
            if requests.get(url).status_code!=200:
                print "PDB structure name: "+ pdbID +" not found"    
                # the program gets stuck in this loop if the PDB is not found, need to address this issue
                # pdbID = str(raw_input("enter another PDB Structure: ")).upper()
                # url = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId="+pdbID
                # pathnamePDB="PDB/"+pdbID+".pdb"
            else:
                print "pdb found"
                urllib.urlretrieve (url, pathnamePDB)
                print "Did the PDB download and save (does file exist)?",os.path.isfile(pathnamePDB)
                return pdbID.upper()

             
grabPDB()
readFromDatafile('1MBD')
getTrpCoordfromPDB(7)
getHemeCoordfromPDB()
#print atomXYZ



fig = matplotlib.pyplot.figure()
ax  = fig.add_subplot(111, projection = '3d')

#plot reads coordinates as [x1,x2,x3...],[y...][z...] so need to split up the function
# there is probably a more elegant way to write this, but this works for now
# make the return work for x, y, z not just one coordinate. I tried 3 but it didnt work
def mergeCoordinatesForPlotting(listOfkeys,index):
	x=[]
	for key in listOfkeys:
		x.append(atomXYZ[key][index])
	return x

#simplify for readibility
mCFP=mergeCoordinatesForPlotting
#plots sets of coordinates
#
ax.plot(mCFP(trpList,0),mCFP(trpList,1),mCFP(trpList,2),color = 'r')
ax.plot(mCFP(hemeList,0),mCFP(hemeList,1),mCFP(hemeList,2),color = 'b')


matplotlib.pyplot.show()