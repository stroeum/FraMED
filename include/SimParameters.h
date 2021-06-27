/* SimParameters.h */
/* SIMulation PARAMETERS */

#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H

#include <cmath>

/**************************************************************************************/
class BoxSteps
{
public:
	int x,y,z;
	BoxSteps();
	BoxSteps(int xx, int yy, int zz);	
	double max();
	bool init(int xx, int yy, int zz);
	~BoxSteps();
}; // Number of steps in the box
/**************************************************************************************/

/**************************************************************************************/
class BoxLengths
{
public:
	double x,y,z;
	BoxLengths();
	BoxLengths(double xx,double yy,double zz);
	bool init(double xx, double yy, double zz);
	~BoxLengths();
}; // Dimenstions of the box
/**************************************************************************************/

/**************************************************************************************/
class StepsSizes
{
public:
	double x,y,z,xy,yz,xz,xyz;
	StepsSizes();									// Default Constructor
	StepsSizes(BoxLengths L, BoxSteps N);			// Surcharged constructor
	bool init(BoxLengths L, BoxSteps N);			// Initiation after declaration
	~StepsSizes();									// Destructor
};
/**************************************************************************************/

#endif SIMPARAMETERS_H