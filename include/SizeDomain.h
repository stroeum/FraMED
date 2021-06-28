/*
 *  SizeDomain.h
 *  Created by Jeremy Riousset on 10/25/07.
 */

#ifndef SIZEDOMAIN_H
#define SIZEDOMAIN_H
#include <iostream>
using namespace std;

/**************************************************************************************/
class SizeDomain
{
public:
	double x,y,z;
	SizeDomain();
	SizeDomain(double xx,double yy,double zz);
	bool init(double xx, double yy, double zz);
	~SizeDomain();
}; // Dimensions of the box
/**************************************************************************************/

#endif
