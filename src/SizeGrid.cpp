/*
 *  SizeGrid.cpp
 *  Created by Jeremy Riousset on 10/25/07.
 */

#include "SizeGrid.h"

/**************************************************************************************/
SizeGrid::SizeGrid(){x = 1; y = 1; z = 1;};

SizeGrid::SizeGrid(int xx, int yy, int zz)
{SizeGrid::init(xx,yy,zz);};

double SizeGrid::max()
{
	double mmax;
	double ttmp = (x >= y)*x + (x < y)*y;
	mmax = (ttmp>= z)*ttmp + (ttmp < z)*z;
	return mmax;
};

bool SizeGrid::init(int xx, int yy, int zz)
{
	x = xx; y = yy; z = zz;
	return true;
};

bool SizeGrid::IsOnBoundary(Point &p)
{
	return(	(p.i == 0)	||	(p.j == 0)	||	(p.k == 0)	||
			(p.i == (x - 1))	||	(p.j == (y - 1))	||	(p.k == (z - 1))	);
}

SizeGrid::~SizeGrid(){};
/**************************************************************************************/
