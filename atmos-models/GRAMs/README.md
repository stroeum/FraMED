# Folder Description
This folder contains more sample input files for FraMED simulations. The .dat files are derived by coupling profile outputs from the planetary suite of NASA GRAMs with the breakdown threshold calculator TREBEC. The provided values can be interpolated for different domain resolutions and sizes with the use of the resizeDomain.m MATLAB file.

# Usage of resizeDomain.m 
Script uses the .dat files derived from coupling the NASA GRAMs with the TREBEC model to create additional, customized sets of input files for the FraMED model through the use of interpolation between critical values.

## Steps for Run-Time:
1. Run resizeDomain.m in MATLAB.
2. Input desired changes into terminal window.
3. (Optional) Define the horizontal domain size based on the spacings derived from Step 2.

### Three options are available to modify the provided .dat files:
- Change the spacings between the nodes (D).
  - Increasing the spacing (e.g. from 25m to 100m) decreases the runtime and resolution for lightning propagation.
  - Decreasing the spacing (e.g. from 25m to 10m) increases the runtime and resolution for lightning propagation.
- Change the maximum altitude of the domain (L).
  - Only allows for a maximum altitude less than the available limit (i.e. script interpolates but does not extrapolate).
- Change both the spacings and the maximum altitude (B).
   
## Resulting Output
### On-Screen:
- List of unavailable files; script checks for the existence of the following potential files (filename suffixes): 
  - Particle density ( '_ng.dat' ).
  - Positive propagation threshold ( '_Eth_positive_Vm.dat' )
  - Negative propagation threshold ( '_Eth_negative_Vm.dat' )
  - Initiation threshold ( '_E_initiation_Vm.dat' )
  - Temperature ( '_T_g.dat' )
  - Reduced electric field ( '_Ek.dat' )
- Path and prefix for associated output files.
- (Optional) Summary table; only displayed if the horizontal domain size is defined.

### Files:
The associated output files will be found within the same object's folder. The files will have the prefix output to the screen, with the same suffixes as the existing potential files.

## Sample MATLAB Command Window:

````{verbatim}
>>>> resizeDomain

What is the planetary body that the simulation is focused on?
(No quotation marks needed for string input)
-->?

	Not an acceptable input.
	Default options: Earth, Mars, Titan, Venus.

What is the planetary body that the simulation is focused on?
(No quotation marks needed for string input)
-->Mars

What type of discharge should it be? (leader/streamer)
-->streamer

*** Domain for Mars has an altitude domain between 0m --> 2500m with 25m spacings in between (total of 101 grid points). ***

Would you like to change the spacings (D), maximum domain altitude (L), or both (B)?
-->B

*** (B) CHANGING MAXIMUM ALTITUDE AND SPACINGS IN-BETWEEN ***

Input the new spacing value in meters (currently 25 meters).
-->50

Input the new maximum altitude in meters (currently 2500 meters).
-->2000

*** New domain for Mars has an altitude domain between 0m --> 2000m with 50m spacings in between (total of 41 grid points). ***
	Temperature file not found for Mars.
	Reduced electric field file not found for Mars.

Would you like to define the number of grid points for the x and y dimensions too?
This is not required but simplifies the use of main.cpp to initialize the simulation. (Y/N)
-->Y

There are currently 41 nodes for the z-dimension.
How many grid points shall the x-dimension span?
-->31

How many grid points shall the y-dimension span?
-->31

	******* SUMMARY OF INITIALIZED SIMULATION DOMAIN *******
				(x)	(y)	(z)
	Domain Size 	(L):	1500	1500	2000	(meters)
	Spacings 	(D):	50	50	50	(meters)
	Grid Points 	(N):	31	31	41	(nodes)

	Total nodes = 39401	(9.8% of previous 400869, should run 10.2x faster)

	Path and prefix of associated files: Mars/streamer-N41-D50m_
>> 
````