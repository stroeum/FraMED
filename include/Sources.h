/* Sources.h */
#ifndef SOURCESBIS_H
#define SOURCESBIS_H

#include <iostream>
#include <iomanip>
#include <list>
#include <string>
#include "Constants.h"
#include "Matrix.h"
#include "SimParameters.h"
using namespace std;

/**************************************************************************************/
class Charge
{
protected:
	/* Charge Type */
	string Type;
	/* Charge Content */
	double Q;
	/* Charge Position */
	double Xq,Yq,Zq;
	/* Charge Geometrical Parameters */
	double Rq1,Rq2,Rq3;
	
public:
	CMatrix3D Un;
	CMatrix3D rho;
	Charge();											// Default constructor
	Charge(StepsSizes dd, BoxSteps NN);					// Initialize Un and rho size
	Charge(double QQ, double XXq, double YYq, double ZZq, double RRq1, double RRq2, double RRq3);
														// Initialize geometrical parameters
	bool		init(double QQ, double XXq, double YYq, double ZZq, double RRq1, double RRq2, double RRq3);
														// Initiate geometrical parameters
	bool		reset(StepsSizes dd, BoxSteps NN);
														// Reset geometrical parameters
	bool		gaussian(double QQ, double XXq, double YYq, double ZZq, double aaq, double bbq, double ccq, StepsSizes dd, BoxSteps NN);
														// Distribute charge assuming Gaussian distribution
	bool        gaussian(double QQ, double XXq, double YYq, double ZZq, double llambda,double mmu,double nnu, double aa,double bb, StepsSizes dd, BoxSteps NN);
														// Distribute charge assuming INCLINED Gaussian distribution
	bool		disk(double QQ, double XXq,double YYq,double ZZq, double RR, double HH, StepsSizes dd, BoxSteps NN);
														// Distribute charge assuming disk geometry
	bool		ellipse(double QQ, double XXq,double YYq,double ZZq, double aa, double bb, double hh, StepsSizes dd, BoxSteps NN);
														// Distribute charge assuming elliptical geometry
	bool		ellipsoid(double QQ, double XXq,double YYq,double ZZq, double aa, double bb, double cc, StepsSizes dd, BoxSteps NN);
														// Distribute charge assuming ellipsoidal geometry
	bool		sphere(double QQ, double XXq, double YYq, double ZZq, double RR, StepsSizes dd, BoxSteps NN);
														// Distribute charge assuming spherical geometry
	bool		rectangle(double QQ, double XXq,double YYq,double ZZq, double llx,double lly,double llz, StepsSizes dd, BoxSteps NN);
														// Distribute charge assuming parallelepipedic geometry
	
	// Calculate analytical solution for a SPHERE of radius Rq1 and charge Q at Xq,Yq,Zq
	CMatrix1D	MonopoleAnalyticalSolution(	StepsSizes dd, BoxSteps NN);
														// ... in free space
	CMatrix1D	DipoleAnalyticalSolution(	StepsSizes dd, BoxSteps NN);
														// ... above a PEC ground plane
	CMatrix1D	MultipoleAnalyticalSolution(StepsSizes dd, BoxSteps NN);
														// ... between two PEC ground planes
	CMatrix1D	getParams();							// Get Charge parameters Q, Xq,Yq,Zq, Rq1,Rq2,Rq3
	string		getType();								// Get Charge type
	
	/*Note: All the surcharged operators below does not consider the position of the 
			centers of the charges, it is basically a shortcut to sum the total charge as
			well as the Status of the lattices. After the sum, the center is affected to
			be 0,0,0 but does not really matter since it won't be used. Our objective while
			defining those operator was to ease the use of several charges.*/
	Charge		operator+=(const Charge&);
	Charge&		operator=(const Charge&);
	Charge		operator+(const Charge&) const;
	friend ostream & operator<< (ostream &, const Charge &);
	~Charge();
};
/**************************************************************************************/

/**************************************************************************************/
class Potential
{
private:
	bool	EquiPotential;							// Type of source
	double	Vo;										// Potential value
	double	Xc,Yc,Zc;								// Center Position
	double	L,W,H;									// Dimensions of the electrode,
													//  * L: Length/Radius
													//  * W: Width
													//  * H: Heigth
public:
	CMatrix3D rho;									// Potential Distribution
	CMatrix3D Un;									// Status of the lattice
													// 0 is modifiable point, 1 if not.
													// Insofar as no point other than the
													// boundary points, Un is no required when
													// source is a charge distribution.
	Potential();									// Defaut constructor
	Potential(CMatrix3D& pphi, CMatrix3D& UUn);		// Constructor by definition of the potential 
													// distribution and Un
	bool init(CMatrix3D& pphi, CMatrix3D& UUn);		// update potential as if constructed in the 
													// previous way
	Potential(double,double,double,double,double,double,double,StepsSizes,BoxSteps);
                                                    // Constructor surcharge - Box
	Potential(double,double,double,double,double,double,StepsSizes,BoxSteps);
                                                    // Constructor surcharge - Vertical cylinder
	Potential(double,double,double,double,double,StepsSizes,BoxSteps);
                                                    // Constructor surcharge - Sphere
	double	getVo()	{return Vo;};					// Return the potential Vo
	bool	getEquiPotential() {return EquiPotential;};	
                                                    // Return EquiPotential
	void	updateUn(const CMatrix3D&);				// Update Un when a new link is added.
	Potential& operator=(const Potential&);			// Equal 2 Potentials distributions.
	friend ostream & operator<< (ostream &, const Potential &);
	~Potential();
};
/**************************************************************************************/

/**************************************************************************************/
typedef list<Charge>   ListCharge;
/**************************************************************************************/
#endif //SOURCESBIS_H
