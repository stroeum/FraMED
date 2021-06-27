/* Links.h */

#ifndef LINKS_H
#define LINKS_H

#include <iostream>
#include <list>
#include <time.h>

#include "Matrix.h"
#include "SimParameters.h"
using namespace std;

/**************************************************************************************/
enum LinkType {x, y, z, xy, yz, xz, xyz};
/**************************************************************************************/

/**************************************************************************************/
/* Point																			  */
/**************************************************************************************/
class Point
{
public:
	int i,j,k;										// Coordinates of the point
	Point(int ii=0, int jj=0, int kk=0)
	{
		i = ii;
		j = jj;
		k = kk;
	}												// Default constructor
	bool operator==(const Point&) const;			// Compare two points, return true if identical
	bool operator!=(const Point&) const;			// Compare two points, return false if identical
	Point& operator=(const Point& Pt);				// Overloading affectation operator
	void init(int ii, int jj, int kk) {i = ii; j = jj; k = kk;};
													// Set point coordinates
	~Point(){};										// Destructor
};
/**************************************************************************************/

/**************************************************************************************/
/* Link																				  */
/**************************************************************************************/
class Link
{
private:
public:
	LinkType type;									// Type of link
													//	* linear along x, y, z axis
													//  * planar in xy, yz or xz
													//  * tridimensional in xyz
	double	l;										// Length of the link
	
	
	Point	start;									// Starting point declared public to
													//		simplify access
	Point	end;									// Ending point declared public to
													//		simplify access
	double	efield;									// Electric field along the link
	double	deltaV;									// Potential drop along the channel
	double	proba;									// Probability of establishment of the link (bounding)
	
	Link(){};										// default constructor
	Link(const Link& LL);							// Constructor surcharge
	Link(Point sstart, Point eend, double ddeltaV, StepsSizes dd, double pphiStart, double pphiEnd);
													// Constructor surcharge
	Link(int SStarti,int SStartj,int SStartk, int EEndi,int EEndj,int EEndk, double ddeltaV, StepsSizes dd, double pphiStart, double pphiEnd);
													// Constructor surcharge	
	Link(Point sstart, Point eend, double ddeltaV, StepsSizes dd, CMatrix3D pphi);
													// Constructor surcharge	
	Link(int iiStart, int jjStart, int kkStart, int iiEnd, int jjEnd, int kkEnd, double ddeltaV, StepsSizes dd, CMatrix3D pphi);
													// Constructor surcharge
	
	void init(const Link& LL);						// Define a link after its declaration
	void init(Point sstart, Point eend, double ddeltaV, StepsSizes dd, double pphiStart, double pphiEnd);
													// Define a link after its declaration
	void init(int SStarti,int SStartj,int SStartk, int EEndi,int EEndj,int EEndk, double ddeltaV, StepsSizes dd, double pphiStart, double pphiEnd);
													// Define a link after its declaration
	void init(Point sstart, Point eend, double ddeltaV, StepsSizes dd, CMatrix3D pphi);	
													// Define a link after its declaration
	void init(int iiStart, int jjStart, int kkStart, int iiEnd, int jjEnd, int kkEnd, double ddeltaV, StepsSizes dd, CMatrix3D pphi);
													// Define a link after its declaration

	void set_deltaV(double ddeltaV){deltaV = ddeltaV;};
													// Set voltage drop
	void fixLink(CMatrix3D&);						// Set value 1 at each extremities of the link in matrix Un
	bool operator==(const Link&) const;				// Compare two Link, return true if identical
	bool operator!=(const Link&) const;				// Compare two Link, return false if identical
	bool isCrossing(const Link&) const;				// Check 2 Links crossing
	bool isCrossing(list<Link>&) const;				// Check if crossing any element of a list.
	Link& operator=(const Link& matrix);			// Assign value of an link to another one.
	friend ostream & operator<< (ostream &, const Link &);
													// Display the link
	~Link(){};										// Destructor
};
/**************************************************************************************/

/**************************************************************************************/
typedef list<Link>   ListLink;
/**************************************************************************************/

#endif LINKS_H