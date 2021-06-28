/*
 *  ResGrid.cpp
 *  Created by Jeremy Riousset on 10/25/07.
 */

#include "ResGrid.h"

/**************************************************************************************/
ResGrid::ResGrid() {x = 1; y = 1; z = 1;};			// Default Constructor

ResGrid::ResGrid(SizeDomain L, SizeGrid N)
{ResGrid::init(L,N);};								// Surcharged constructor

bool ResGrid::init(SizeDomain L, SizeGrid N)
{
	x	= L.x/(N.x-1);
	y	= L.y/(N.y-1);
	z	= L.z/(N.z-1);
	xy	= sqrt(pow(x,2) + pow(y,2));
	yz	= sqrt(pow(y,2) + pow(z,2));
	xz	= sqrt(pow(x,2) + pow(z,2));
	xyz	= sqrt(pow(x,2) + pow(y,2) + pow(z,2));
	return true;
}												// Initiation after declaration

/* SAM function.  Similar to 'init(SizeDomain, SizeGrid)'.  Rationale:
 *
 * The input file generated from the 'buck' program outputs grid dimensions
 * and grid lengths.  Hence it is easier to initialize the ResGrid directly from
 * the grid lengths, and not from the SizeDomain and SizeGrid (as was done
 * previously).
 */
bool ResGrid::init(double xi, double yi, double zi)
{
	x = xi;
	y = yi;
	z = zi;
	xy = sqrt(pow(xi, 2) + pow(yi, 2));
	yz = sqrt(pow(yi, 2) + pow(zi, 2));
	xz = sqrt(pow(xi, 2) + pow(zi, 2));
	xyz = sqrt(pow(xi, 2) + pow(yi, 2) + pow(zi, 2));
	return true;
}

ResGrid::~ResGrid(){};
/**************************************************************************************/