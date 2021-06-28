/* Links.h */

#ifndef LINKS_H
#define LINKS_H

#include <list>
#include "Point.h"
#include "Matrix.h"
#include "ResGrid.h"
using namespace std;

/**************************************************************************************/
enum LinkType {x, y, z, xy, yz, xz, xyz};
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


	Point	start;									// Starting point declared public to simplify access
	Point	end;									// Ending point declared public to simplify access
	double	efield;									// Electric field along the link
	double	deltaV;									// Potential drop along the channel
	double	proba;									// Probability of establishment of the link (bounding)

	Link(){};										// default constructor
	Link(const Link& lk);							// Constructor surcharge
	Link(Point A, Point B, double dV, ResGrid d, double phiStart, double phiEnd);
													// Constructor surcharge
	Link(int Ai,int Aj,int Ak, int Bi,int Bj,int Bk, double dV, ResGrid d, double phiStart, double phiEnd);
													// Constructor surcharge
	Link(Point A, Point B, double dV, ResGrid d, CMatrix3D phi);
													// Constructor surcharge
	Link(int iA, int jA, int kA, int iB, int jB, int kB, double dV, ResGrid d, CMatrix3D phi);
													// Constructor surcharge

	void init(const Link& lk);						// Define a link after its declaration
	void init(Point A, Point B, double dV, ResGrid d, double phiStart, double phiEnd);
													// Define a link after its declaration
	void init(int Ai,int Aj,int Ak, int Bi,int Bj,int Bk, double dV, ResGrid d, double phiStart, double phiEnd);
													// Define a link after its declaration
	void init(Point A, Point B, double dV, ResGrid d, CMatrix3D phi);
													// Define a link after its declaration
	void init(int iA, int jA, int kA, int iB, int jB, int kB, double dV, ResGrid d, CMatrix3D phi);
													// Define a link after its declaration

	void set_deltaV(double dV){deltaV = dV;};       // Set voltage drop
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

#endif
