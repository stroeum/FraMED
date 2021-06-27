/* CriticalFields.h */

#ifndef CRITICALFIELDS_H
#define CRITICALFIELDS_H

#include "Matrix.h"
#include "SimParameters.h"
using namespace std;

/**************************************************************************************/
class ScaledFields
{
protected:
	double z_gnd;   								// Altitude of the ground plane
	double positive0;								// Value of the propagation positive threshold at z = z_gnd
	double negative0;								// Value of the propagation negative threshold at z = z_gnd
public:
	CMatrix1D positive;								// Positive propagation threshold
	CMatrix1D negative;								// Negative propagation threshold
	
	ScaledFields(){};								// Default constructor
	ScaledFields(double ppositive0, double nnegative0,
				   double zz_gnd, StepsSizes dd, BoxSteps NN);
													// Constructor surcharge
	CMatrix1D getParams();							// Retrieve z_gnd, initiation0,positive0,Negative0
	bool init(double ppositive0, double nnegative0,
			  double zz_gnd, StepsSizes dd, BoxSteps NN);
													// Initiate after declaration
	~ScaledFields(){};								// Destructor
};
/**************************************************************************************/

/**************************************************************************************/
class VoltageDrops : public ScaledFields
{
public:
	VoltageDrops(){};
	VoltageDrops(double ppositive0, double nnegative0, double zz_gnd, StepsSizes dd, BoxSteps NN) :
		ScaledFields(ppositive0,nnegative0,zz_gnd,dd,NN){};
													// Constructor
	~VoltageDrops(){};								// Destructor
};
/**************************************************************************************/


/**************************************************************************************/
class CriticalFields : public ScaledFields
{
protected:
	double initiation0;								// Value of the initiation field at z = z_gnd
public:
	CMatrix1D initiation;							// Initiation Threshold

	CriticalFields(){};
	CriticalFields(double iinitiation0, double ppositive0, double nnegative0,
				   double zz_gnd, StepsSizes dd, BoxSteps NN);
													// Constructor
	CMatrix1D getParams();							// Retrieve z_gnd, initiation0,positive0,Negative0
	bool init(double iinitiation0, double ppositive0, double nnegative0,
			  double zz_gnd, StepsSizes dd, BoxSteps NN);
													// Initiate after declaration
	~CriticalFields(){};							// Destructor
};
/**************************************************************************************/

/**************************************************************************************/
double StdAtmosphere(double);
/**************************************************************************************/
#endif CRITICALFIELDS_H
