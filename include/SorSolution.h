/* SorSolution.h */

#ifndef SORSOLUTION_H
#define SORSOLUTION_H

#include "Sources.h"

/**************************************************************************************/
enum SourceType {ChargeDistribution, PotentialDistribution};
													// Allowed type of sources to use SOR solver
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
	SorSolution(CMatrix3D&, double, int, ResGrid, SizeGrid, Charge&, const CMatrix3D&);
													// Constructor surcharge
	SorSolution(CMatrix3D&, double, int, ResGrid, SizeGrid, Potential&, CMatrix3D&);
													// Constructor surcharge
	~SorSolution(){};								// Destructor
	void init(CMatrix3D&, double, int, ResGrid, SizeGrid, Charge&, const CMatrix3D&);
	void init(CMatrix3D&, double, int, ResGrid, SizeGrid, Potential&, CMatrix3D&);
	void Solve(ResGrid _d, SizeGrid N, const CMatrix3D& Un, CMatrix3D& phi);
													// Solve solution using SOR method
};
/**************************************************************************************/

#endif
