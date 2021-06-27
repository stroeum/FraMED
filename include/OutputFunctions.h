/* OutputFunctions.h */
#ifndef OUTPUTFUNCTIONS_H
#define OUTPUTFUNCTIONS_H

#include "Matrix.h"
#include "SimParameters.h"
#include "Constants.h"
#include "CriticalFields.h"
#include "Trees.h"
#include "Links.h"
#include <fstream>
#include <math.h>

struct Vector {double x,y,z;};
typedef list<Vector> ListVector;
typedef list<double> ListDouble;
typedef list<CMatrix1D> ListCMatrix1D;

/**************************************************************************************/
/* For Electric Field																  */
/**************************************************************************************/
CMatrix1D Eijk(int, int, int,const CMatrix3D&, const StepsSizes&, const BoxSteps&);
													// Derive local Ex, Ey, Ez at point
													//	 i,j,k as well as its magnitude.
CMatrix3D GlobalE(const CMatrix3D&,const StepsSizes&, const BoxSteps&);				
													// Derive magnitude of E everywhere
CMatrix1D eFieldFlux(const CMatrix3D& pphi, const StepsSizes& dd, const BoxSteps& NN);
													// Derive flux of E-field through the boundaries in C.
/**************************************************************************************/

/**************************************************************************************/
/* For Charge																		  */
/**************************************************************************************/
double rhoijk(int, int, int,const CMatrix3D&, const StepsSizes&, const BoxSteps&);
													// Derive the charge distribution
													// at point i,j,k
CMatrix3D Globalrho(const CMatrix3D&,const StepsSizes&, const BoxSteps&);
													// Derive the charge distribution
double ChannelCharge(const CMatrix3D& rrho, const CMatrix3D& UUn, const StepsSizes& dd, const BoxSteps& NN);
double ChannelChargePositive(const CMatrix3D& rrho, const CMatrix3D& UUn, const StepsSizes& dd, const BoxSteps& NN);
double ChannelChargeNegative(const CMatrix3D& rrho, const CMatrix3D& UUn, const StepsSizes& dd, const BoxSteps& NN);
double TotalCharge(const CMatrix3D& rrho, const StepsSizes& dd, const BoxSteps& NN);
CMatrix1D ChannelLinearDensity(const CMatrix3D& rrho, const CMatrix3D& UUn, const StepsSizes& dd, const BoxSteps& NN);
/**************************************************************************************/

/**************************************************************************************/
/* For Dipole Moment																  */
/**************************************************************************************/
Vector DipoleMoment(double& CCarriedCharge, const CMatrix3D& pphi, const CMatrix3D& UUn, const BoxLengths& LL, const StepsSizes& dd, const BoxSteps& NN);
													// Also derive charge carried by the channel
/**************************************************************************************/

/**************************************************************************************/
/* Misc.																			  */
/**************************************************************************************/
void SwitchValues(double& x1, double& x2);
void SwitchValues(int& n1, int& n2);
bool isfinite(const CMatrix3D& MM, const BoxSteps& NN);

ListDouble read(string);
void write(ListDouble&, string);
void write(ListVector&, string);
void write(ListCMatrix1D&, string);
/**************************************************************************************/


#endif OUTPUTFUNCTIONS_H
