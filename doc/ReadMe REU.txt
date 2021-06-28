FractalModel v1.0
Jeremy Riousset*
Samuel T. Poulos**
CSSL Pennsylvania State University, College Station

* Student at the Pennsylvania State University at State College, PA
** Undergraduate Student at the University of Texas at Austin, TX

Throughout the course of the REU program, I (Sam) concentrated on several key aspects of this program.  Changes to parts of the program relating to these aspects are described below:

New Method of Input

In addition to allowing the user to specify charge configuration via cylinders of uniform charge and other predefined objects, charge can now be specified via an input file generated from the program 'buck' (written by Sam).  The process by which this is accomplished is described below:

XLMA* plain text output --(buck)--> file of discretized points --(new method in Sources.cpp)--> charge density matrix and grid attributes initialized

The format of the file of discretized points is as follows:

%d %d %d			grid dimensions (xn, yn, zn)
%lf %lf %lf			grid side lengths (xl, yl, zl) (i.e., locations [1, 1, 2] and [1, 1, 3] are separated by zl meters) 
%d				altitude (above sealevel) of ground
(%d %d %d %lf)_n		source point data (xi, yi, zi, charge)

The program 'buck', as well as XLMA, should be freely available.

*program written to analyze data from the L-ightning M-apping A-rray


Formatted Output

Previously, several important interim messages were displayed to the console as follows:


/******************************
/* Initiating Tree		 *
/******************************/

Other messages were displayed with no special formatting at all.  In an attempt to make the output more meaningful and easy to understand, the following descriptions were appended to output messages such as the following:

ii: Maximum charge density: 9.84 nC/m^3.					(ii = information)
ii:	 Maximum charge density occurred at [4.5 km, 3.9 km, 9.8 km].	(tab denotes subpoint)
..: Initiating the tree.									(.. = beginning task)
++: Tree initiated!									(++ = success)
ss: Saving intermediate results.							(ss = saving)
ee: Boundary not specified.							(ee = error)
xx: Program terminating.								(xx = termination)


File IO Subfolder Modification

Previously, every IO::write function wrote to the file specified in its argument.  Now, the IO class owns a static member pointing to a string containing the path name, and every IO::write function first appends its argument to that path name in order to calculate the complete file name.  The path name is set through two function calls: SetPathName() and SortPathName().  SetPathName() creates a new directory based on the current time of day and sets the IO path name accordingly.  SortPathName creates a new directory based on the time of day at the beginning of the simulation run and on the number of established links.  It too sets the IO path name accordingly.

In addition, the IO class now has functions which allow the previous subfolder to be copied into the new directory, and the previous subfolder to be deleted.  It is by this method that the interim results file is copied to its final location, based on the number of established links in the current simulation run.


Miscellaneous

Various sections of code were timed, with the running time output to the screen via function calls such as:

printf("Time required to write results: %lf s.\n", running_time/10000);

However, the time (as returned by calls to the clock() function) were given in milliseconds, thereby overcalculating the running time by a factor of 10.  This problem was remedied with the CLOCKS_PER_SEC macro as follows:

printf("Time required to write results: %lf s.\n", running_time/CLOCKS_PER_SEC);

