/*
 *  Point.cpp
 *  Created by Jeremy Riousset on 10/25/07.
 */

#include "Point.h"

void Point::init(int ii, int jj, int kk)
{
	i = ii;
	j = jj;
	k = kk;
}

bool Point::operator==(const Point& Pt) const
{
	if (i == Pt.i && j == Pt.j && k == Pt.k)
		return true;
	else
		return false;
}  // Operator ==

bool Point::operator!=(const Point& Pt) const
{
	if (i != Pt.i || j != Pt.j || k != Pt.k)
		return true;
	else
		return false;
}  // Operator !=

Point& Point::operator=(const Point& Pt)
{
	i = Pt.i;
	j = Pt.j;
	k = Pt.k;
	return *this;
} // Operator =

double Point::di;
double Point::dj;
double Point::dk;
double Point::oi;
double Point::oj;
double Point::ok;

void Point::initDist(double idi, double idj, double idk, double ioi, double ioj, double iok)
{
	di = idi;
	dj = idj;
	dk = idk;
	oi = ioi;
	oj = ioj;
	ok = iok;
}
