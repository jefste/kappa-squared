import matplotlib.pyplot
from mpl_toolkits.mplot3d import Axes3D

#by adding this line, i can get rid of all my other code that was copied from the other python script
from kappaSquared import *
             
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
def addListtoPlot(listOfkeys,plotname,linecolor):
	plotname.plot(mCFP(listOfkeys,0),mCFP(listOfkeys,1),mCFP(listOfkeys,2),color = linecolor)


addListtoPlot(['CD1','NE1','CE2','CD2'],ax,'r')
addListtoPlot(['NE1','trpOrigin'],ax,'b')
addListtoPlot(['FE','trpOrigin'],ax,'g')

addListtoPlot(['CHC','FE','CHD'],ax,'r')
#ax.plot(mCFP(hemeList,0),mCFP(hemeList,1),mCFP(hemeList,2),color = 'b')



matplotlib.pyplot.show()