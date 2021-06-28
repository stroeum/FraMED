/*
 *  ResGrid.h
 *  Created by Jeremy Riousset on 10/25/07.
 */

#ifndef RESGRID_H
#define RESGRID_H

#include <cmath>
#include "SizeDomain.h"
#include "SizeGrid.h"

/**************************************************************************************/
class ResGrid
{
public:
	double x,y,z,xy,yz,xz,xyz;
	ResGrid();									// Default Constructor
	ResGrid(SizeDomain L, SizeGrid N);			// Surcharged constructor
	/* SAM function.  Similar to 'init(SizeDomain, SizeGrid)'.  Rationale:
	 *
	 * The input file generated from the 'buck' program outputs grid dimensions
	 * and grid lengths.  Hence it is easier to initialize the ResGrid directly from
	 * the grid lengths, and not from the SizeDomain and SizeGrid (as was done
	 * previously).
	 */
	bool init(double xi, double yi, double zi);
	bool init(SizeDomain L, SizeGrid N);		// Initiation after declaration
	~ResGrid();									// Destructor
};
/**************************************************************************************/

#endif
