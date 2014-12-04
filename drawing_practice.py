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
ax.plot(mCFP(hemeList,0),mCFP(hemeList,1),mCFP(hemeList,2),color = 'b')


#shows plot
#matplotlib.pyplot.show()

# add function to read the 'CONECT' lines in 


def parseConnectList(listtoparse,startCoordinate):
    #change this function to only add the entry before and after the atom of interest. 
    #don't need every single connection, only the ones that are connected to the atom of interest.    
    mylist=[]    
    print startCoordinate
    for row in listtoparse:
        for item in row:
            #adds line, but omits 'Conect' command.        
            if str(item) != 'CONECT':
                mylist.append(item)  
        #add newline as place holder for now. need to decide best way forward on combining the atom numbers.q
        mylist.append('newline')
    #creates list of atoms
    return mylist



def listToDraw():
    connections=[]
    hemeListFull=[]
    for row in data:
        if row[0]=='CONECT':
            connections.append(row)

    for row in data:
        if row[0]=='HETATM':
            hemeListFull.append(row)
    
#    for i in range(5,6):
#        print hemeListFull[i][1],hemeListFull[i]
    # prints the row if the index is found. just use FE to see if logic works. (seems to)
    connectList= [row for row in connections if hemeListFull[5][1] in row]
    print parseConnectList(connectList,hemeListFull[5][1])

    

listToDraw()
#print hemeListFull
#print connections
# add function to grab all heme coordinates and all trp coordinates