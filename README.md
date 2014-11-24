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

*calculates the kappa squared (also known as the orientation factor) between 
  the transition dipoles for heme and tryptophan in myoglobin.
  *calculates the kappa squared for both the normal and disordered transition dipoles:
    *estimate: uses only the coordinates of the pdb file
    *'exact': generates coordinates of the file for specific angles


Possible things to implement:
*calculate distance between tryptophan and heme dipole origins (easy to add)
*generalize program: 
  *use any pdb, not just 1MBD
  *tryptophan can be selected (for instance 7 or 14 in wild type Mb)
*get PDB from web/local data base
*generate fish plot for heme
*use mplot3d/matplotlib to draw the tryptophan and heme structures, draw the transition dipoles for each as well as the transition vector
*read from large MD file to map distance and kappa squared over time (need old MD files for this)
*determine if hemes are normal or disordered depending on structure
*create CSV (eventually database?) that contains trp-heme kappa squared for other heme proteins, compare with other literature (FRET book?)
*create local file where information gets stored
