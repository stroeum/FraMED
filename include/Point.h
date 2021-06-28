/*
 *  Point.h
 *  Created by Jeremy Riousset on 10/25/07.
 */

#ifndef POINT_H
#define POINT_H

class Point
{
public:
	int i,j,k;										// Coordinates of the point
	Point(int ii=0, int jj=0, int kk=0)				// Default constructor
	{
		init(ii, jj, kk);
	}
	void init(int ii, int jj, int kk);
	bool operator==(const Point&) const;			// Compare two points, return true if identical
	bool operator!=(const Point&) const;			// Compare two points, return false if identical
	Point& operator=(const Point& Pt);				// Overloading affectation operator

	/* SAM functions */

	/* Each point can now return its exact position in units of length.  However,
	 * this requires the class variables, which specify the distance between grid points
	 * and offset of the grid, to be set.
	 */
	double GetX(){return i*di+oi;}
	double GetY(){return j*dj+oj;}
	double GetZ(){return k*dk+ok;}

	double GetXk(){return GetX()/1e3;} // k = (in) kilometers
	double GetYk(){return GetY()/1e3;}
	double GetZk(){return GetZ()/1e3;}

	// Set point coordinates
	~Point(){};										// Destructor

	// Distances between discretized locations
	static double di, dj, dk;
	// Offsets
	static double oi, oj, ok;
	// Set the aforementioned class variables.
	static void initDist(double idi, double idj, double idk, double ioi, double ioj, double iok);
};

#endif
