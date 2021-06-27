/* SorSolution.h */

#ifndef SORSOLUTION_H
#define SORSOLUTION_H

#include <string>
#include <time.h>
#include "Matrix.h"
#include "SimParameters.h"
#include "Constants.h"
#include "Sources.h"

/**************************************************************************************/
enum SourceType {ChargeDistribution, PotentialDistribution};
													// Allowed type of sources to use
													//  SOR solver
/**************************************************************************************/

/**************************************************************************************/
class SorSolution
{
private:
	double epsilon;									// Maximum tolerable error
	int MaxStep;									// Maximum number of iterations
	double a,b,c,d,e,f,g;							// Laplace's Equation coefficients
	CMatrix3D h;									// Laplace's Equation coefficients (cont.)
	SourceType type;								// Type of source
	double eErrDen;									// Normalisation of the error

public:	
	SorSolution(){};								// Default constructor
	SorSolution(CMatrix3D&, double,int,StepsSizes,BoxSteps, Charge&, const CMatrix3D&);
													// Constructor surcharge
	SorSolution(CMatrix3D&, double,int,StepsSizes,BoxSteps, Potential&, CMatrix3D&);
													// Constructor surcharge
	~SorSolution(){};								// Destructor
	void init(CMatrix3D&, double,int,StepsSizes,BoxSteps, Charge&, const CMatrix3D&);
	void init(CMatrix3D&, double,int,StepsSizes,BoxSteps, Potential&, CMatrix3D&);
	void Solve(StepsSizes dd, BoxSteps NN, const CMatrix3D& UUUn, CMatrix3D& ppphi);
													// Solve solution using SOR method
};
/**************************************************************************************/

#endif SORSOLUTION_H

