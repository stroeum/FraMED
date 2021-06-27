#include "SimParameters.h"
#include <cmath>

/**************************************************************************************/
BoxSteps::BoxSteps(){x = 1; y = 1; z = 1;};

BoxSteps::BoxSteps(int xx, int yy, int zz)
{BoxSteps::init(xx,yy,zz);};

double BoxSteps::max()
{
	double mmax;
	double ttmp = (x >= y)*x + (x < y)*y;
	mmax = (ttmp>= z)*ttmp + (ttmp < z)*z;
	return mmax;
};

bool BoxSteps::init(int xx, int yy, int zz)
{
	x = xx; y = yy; z = zz;
	return true;
};

BoxSteps::~BoxSteps(){};
/**************************************************************************************/

/**************************************************************************************/
BoxLengths::BoxLengths(){x = 1; y = 1; z = 1;};
BoxLengths::BoxLengths(double xx,double yy,double zz)
{BoxLengths::init(xx,yy,zz);};

bool BoxLengths::init(double xx, double yy, double zz)
{
	x = xx; y = yy; z = zz;
	return true;
};

BoxLengths::~BoxLengths(){};
/**************************************************************************************/

/**************************************************************************************/
StepsSizes::StepsSizes() {x = 1; y = 1; z = 1;};			// Default Constructor

StepsSizes::StepsSizes(BoxLengths L, BoxSteps N)
{StepsSizes::init(L,N);};									// Surcharged constructor

bool StepsSizes::init(BoxLengths L, BoxSteps N)
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

StepsSizes::~StepsSizes(){};
/**************************************************************************************/
