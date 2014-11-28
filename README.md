kappa-squared
=============
This python script computes the orientation factor or ''kappa squared'' (k^2) between 2 dipoles in 
Forster Resonance Energy Transfer (FRET). A PDB file is read and the user inputs the number of 
tryptopahn that the orientation factor between the tryptophan and the heme. The orientation factor 
is calculated with 2 different methods, for both the normal and disordered orientations of hemes. 
For further information, read: 
"Intrinsic fluorescence of hemoglobins and myoglobins." 
Gryczynski, Z.; Lubkowski, J.; Bucci, E. Methods Enzymol. 1997, 278, 538–569. 
Specifically, start reading at page 548.
DOI: 10.1016/S0076-6879(97)78030-3

Currently the program:
*calculates the kappa squared (also known as the orientation factor) between the transition dipoles for heme and tryptophan in myoglobin.
*calculates the kappa squared for both the normal and disordered transition dipoles:
    *estimate: uses only the coordinates of the pdb file
    *'exact': generates coordinates of the file for specific angles
*calculates distance between tryptophan and heme dipole origins 
*allows user to select which PDB to determine and which tryptophan to measure.
	*if the file does not exist locally, program will download the PDB from the PDB db
	*if the file does exist locally, the program checks a database file (csv file) to see if your data has already been calculated
	*if no tryptophan is specified, the program will determine at what locations the tryptophans are stored and allow you to select one of them
*has some error handling

TO RUN THIS PROGRAM:
The files/directories needed are 'kappa-squared.py', 'k_sq_db.csv' and a directory 'PDB' in the same directory as the .py and .csv file.
Example to execute this script:
$ python kappa-squared.py 1MBD 7
This will calculate kappa squared for the pdb file 1MBD for tryptophan 7.
Example of other valid usage:
$ python kappa-squared.py 1MBD 
This will prompt user to input which tryptophan to use to calculate kappa squared.
Example of other valid usage:
$ python kappa-squared.py
This will prompt user to input a PDB and tryptophan to use to calculate kappa squared.

  

Future things to implement:
*error handling:
	*handle invalid pdb name
	*handle invalid pdb file (abort if no trp or heme are found)
		*currently aborts script if there are no trps found
*generate fish plot for each pdb
*use mplot3d/matplotlib to draw the tryptophan and heme structures, draw the transition dipoles for each as well as the transition vector
*read from large MD file to map distance and kappa squared over time (need old MD files for this)
*determine if hemes are normal or disordered depending on structure?
*create relation database (SQL?) that contains trp-heme kappa squared for all heme proteins that stores the image for each
	*compare with other literature (FRET book?)
