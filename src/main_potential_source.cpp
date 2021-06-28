/* main.cpp Sept. 5, 2006, 3.45PM  */

/**************************************************************************************/
/* Include libraries and header files												  */
/**************************************************************************************/
#include <time.h>
#include <list>
#include <iostream>
#include <fstream>

#include "BoundaryConditions.h"
#include "Constants.h"
#include "CriticalFields.h"
#include "Links.h"
#include "Matrix.h"
#include "OutputFunctions.h"
#include "SimParameters.h"
#include "SorSolution.h"
#include "Sources.h"
#include "Trees.h"
#include <random>

using namespace std;
/**************************************************************************************/

/**************************************************************************************/
/* Declare functions																  */
/**************************************************************************************/
double EoverEk(const double tt,
ListCharge& CChargeCfg,
const double II1, const double II2, const double IIscreen,
SorSolution SSOR, const double eepsilon, const int MMaxStep,
const int BBCtype, CriticalFields OOverShotEc,
CMatrix3D& pphi, CMatrix3D& UUn,
const StepsSizes dd, const BoxSteps NN);

void LoadTripoleModel(ListCharge& CChargeCfg,
const double II1, const double II2, const double IIscreen,
SorSolution SSOR, const double eepsilon, const int MMaxStep,
const int BBCtype, CriticalFields EEc,
const double TThresholdOvershoot,
CMatrix3D& pphi, CMatrix3D& UUn,
const StepsSizes dd, const BoxSteps NN);

bool InitiateTree(const int IInitiationType,
Point& IInitiationPoint,
double& pphi0, const double TThresholdOvershoot,
CriticalFields EEc, CMatrix1D& EEzNum, CMatrix1D& pphiNum,
SorSolution SSOR,
CMatrix3D& pphi, CMatrix3D& UUn,
const StepsSizes dd, const BoxSteps NN);

void GrowTree(bool AAddNew, const bool iisBndXingPossible, const bool iisFlashAccoutedInBC,
const bool iisRsDeveloped,
const bool iisLinkXingPossible,  const bool iisChannelEquipotential,
const bool iisBCerrorCalculated, const bool iisEsEnergyCalculated,
const int BBCtype, const int rrho3dCalculationStep,
double& pphi0, double VVmin, double VVmax,
const double rrhoAmbMin, const double rrhoAmbMax,
double& QQchannelPlus, double& QQchannelMinus, double& QQtot_af,
CriticalFields& EEc, VoltageDrops& VVd,
const double TThresholdOvershoot, const double IInitR,
const Point& IInitiationPoint, ListLink& EEstablishedLinks,
ListVector& BBndUpdateErrors, ListVector& DDischargeDipoleMoment,
ListDouble& CCarriedCharge, ListDouble& CChannelPotential, ListDouble& EEsEnergy,
ListCMatrix1D& TTotalPotential, ListCMatrix1D& TTotalEfield,
SorSolution SSOR, const double eepsilon, const int MMaxStep,
CMatrix3D& pphi_cha, CMatrix3D& pphi_amb,
CMatrix3D& pphi, CMatrix3D& UUn,
const BoxLengths& LL, const StepsSizes dd, const BoxSteps NN, ListCharge& CChargeCfg);

void StoreData(const bool iisBCerrorCalculated, const bool iisEsEnergyCalculated,
const Point& IInitiationPoint, const double IInitR, ListCharge& CChargeCfg,
ListLink EEstablishedLinks,
ListVector BBndUpdateErrors, ListVector DDischargeDipoleMoment,
ListDouble CCarriedCharge, ListDouble CChannelPotential, ListDouble EEsEnergy,
ListCMatrix1D& TTotalPotential, ListCMatrix1D& TTotalEfield,
CriticalFields EEc, CMatrix1D& EEzNum, CMatrix1D& pphiNum,
CMatrix3D& EE, CMatrix3D& pphi, CMatrix3D& pphi_amb,
const StepsSizes dd, const BoxSteps NN);

void StoreIntermData(const Point& IInitiationPoint, const double IInitR, ListCharge& CChargeCfg,
ListLink EEstablishedLinks,
CriticalFields EEc, CMatrix3D& pphi, CMatrix3D& pphi_amb,
const StepsSizes dd, const BoxSteps NN, int cccptLinks);
/**************************************************************************************/

/**************************************************************************************/
/* Main																				  */
/**************************************************************************************/
int main()
{
time_t seconds;

time(&seconds);

srand((unsigned int) seconds);


/* The next lines test the random number generator */
//    double rr;
//    for (int n=0; n<10; n++) {
//        rr = rand()/(double)RAND_MAX;
//        printf("r = %f \n",rr);
//    };
//    exit(121);

clock_t startTime = clock();

/**********************************************************************************/
/* Variables Declaration														  */
/**********************************************************************************/
cout<<"\nstart <Equipotential lightning flash>\n\n";
BoxSteps		N;						// Number of discretization points
BoxLengths		L;						// Dimensions of the simulation domain
StepsSizes		d;						// Lengths of the discretization-grid

double			z_gnd;					// Altitude of ground plane
double			z_shift;				// Vertical displacement of the cloud
double			y_shift;				// Horizontal displacement of the upper part of the cloud

const double	epsilon	= 1e-10;		// SOR precision
const int		MaxStep	= 7500;			// Allowed maximum number of point per SOR iteration

double			InitX;					// X-coordinate of the initiation point
double			InitY;					// Y-coordinate of the initiation point
double			InitZ;					// Z-coordinate of the initiation point
Point			InitiationPoint;		// i,j,k coordinates of the initiation point
double			InitR = 0.0;			// Radius of the initiation region (if applicable)

// The following values are only use to avoid warning at compilation //
double			phi0			= 0;	// Channel potential after minization
double			Vmin			= 0;	// Minimum in cloud potential
double			Vmax			= 0;	// Maximum in cloud potential
double			rhoAmbMin		= 0;	// Minimum in cloud charge density
double			rhoAmbMax		= 0;	// Maximum in cloud charge density
double			QchannelPlus	= 0;	// Positive charge content of the channel (== charge transfer)
double			QchannelMinus	= 0;	// Negative charge content of the channel
double			Eps_bf			= 0;	// Electrostatic energy before the discharge
double			Eps_af			= 0;	// Electrostatic energy after the discharge
double			Qtot_bf			= 0;	// Total charge in the simulation domain before the discharge
double			Qtot_af			= 0;	// Total charge in the simulation domain after the discharge
// The previous  values are only use to avoid warning at compilation //

double			I1, I2, Iscreen;		// Loading currents I1 and I2, and screening current Iscreen
double			Vo, Xp,Yp,Zp, Rp1,Rp2,Rp3; // Potential center parameters:
// * Potential (V),
// * X-,Y-,Z-coordinates of a potential center (Xp,Yp,Zp),
// * 1st, 2nd and 3rd geometrical parameters (Rp1, Rp2, Rp3)
double			Q,  Xq,Yq,Zq, Rq1,Rq2,Rq3; // Charge center parameters:
// * Charge (Q),
// * X-,Y-,Z-coordinates of a charge center (Xq,Yq,Zq),
// * 1st, 2nd and 3rd geometrical parameters (Rq1, Rq2, Rq3)
double			ThresholdOvershoot;		// %-age by which threshold must be exceeded to initiate discharge
int				BCtype = 0;				// Boundary conditions:
// = 0 ``Tin Can"
// = 1 Open BC
// = 2 Moving Capacitor Plates
// = 3 Free Space
int				InitiationType;			// Initiation type:
// = 1 randomly in available region
// = 2 at point of maximum E-field magnitude
// = 3 at implemented position
int				rho3dCalculationStep    = 0;    // Charge density is calculated and store every rho3dCalculationStep steps. (0: 3-D charge density is never calculated)

bool			isChargeContentAssumed = false;	// Choose either charge content or loading currents
bool			isInitiationPrevented  = false;	// Only simulate cloud electrical structure
bool			AddNew                 = false;	// Channel is allowed to propagate: Y/N
bool			isBCerrorCalculated    = false;	// Error at the boundary is calculated at each step: Y/N
bool			isBndXingPossible      = false;	// Channel is allowed to cross boundaries: Y/N
bool			isChannelEquipotential = false;	// Channel potential is adapted to ensure charge neutrality: Y/N
bool			isEsEnergyCalculated   = false;	// Electrostatic energy is calculated: Y/N
bool			isFlashAccoutedInBC    = false;	// Channel charge is accounted for in derivation of BC: Y/N
bool			isInitiationPossible   = false;	// Initiation is possible in the simulation domain after charge load: Y/N
bool			isLinkXingPossible     = false;	// Channels crosses are allowed: Y/N
bool			isRsDeveloped          = false;	// Return stroke development: Y/N

Charge			C;						// Charge center
SorSolution		SOR;					// Parameters of SOR algorithms

ListDouble		ChannelPotential;		// Channel potential at each stage of development
ListDouble		CarriedCharge;			// Channel transfer at each of stage development
ListDouble		EsEnergy;				// Electrostatic energy in the simulation domain
ListLink		EstablishedLinks;		// List of links constituting the discharge tree
ListVector		BndUpdateErrors;		// List of errors at the boundaries
ListVector		DischargeDipoleMoment;	// Channel dipole moment at each of stage development
ListCMatrix1D	TotalPotential;			// Potential on a vertical axis at the center of the simulation domain
ListCMatrix1D	TotalEfield;			// E-field on a vertical axis at the center of the simulation domain
ListCharge		ChargeCfg;				// Table with all parameters of the charge configuation

ListCharge::iterator it;				// Table with all parameters of the charge configuation

CriticalFields	Ec;						// Initiation, propagation of positive and negative channels
VoltageDrops	Vd;						// Voltage drops in positive and negative channels

CMatrix3D		E;						// _V/m	Total electric field
CMatrix3D		phi;					// _V	Total electric potential
CMatrix3D		phi_cha;				// _V	Channel electric potential
CMatrix3D		phi_amb;				// _V	Cloud electric potential
CMatrix3D		rho;					// _C/m3	total charge density
CMatrix3D		Un;						// Map of occupied grid points
CMatrix1D		phiNum;					// _V	Total electric potential on a vertical axis in the center of simulation domain
CMatrix1D		EzNum;					// _V/m	Total electric field on a vertical axis in the center of simulation domain
CMatrix1D		eFlux(7);				// _C Flux of electrostatic field through the boundaries
/**********************************************************************************/
/* Endof Variables Declaration													  */
/**********************************************************************************/

N.init(51,51,51);					// Number of discretization points
L.init(10e+3,10e+3,10e+3);			// _m	Size of the simulation domain
d.init(L,N);						// _m	Sizes of grid-steps

// Initiate matrices dimensions //
E.init(N.x,N.y,N.z);				// _V/m	Total electric field
phi.init(N.x,N.y,N.z);				// _V	Total electric potential
phi_cha.init(N.x,N.y,N.z);			// _V	Channel electric potential
phi_amb.init(N.x,N.y,N.z);			// _V	Cloud electric potential
rho.init(N.x,N.y,N.z);				// _C/m3	total charge density
Un.init(N.x,N.y,N.z);				// Map of occupied grid points
phiNum.init(N.z);					// _V	Total electric potential on a vertical axis in the center of simulation domain
EzNum.init(N.z);					// _V/m	Total electric field on a vertical axis in the center of simulation domain
// Endof Initiate matrices dimensions //

z_gnd   = 0e3;						// _m Altitude of ground level
z_shift = 0e3;						// _m Vertical displacement of the cloud
y_shift = 0e3;						// _m Horizontal influence of a windshear

//	We assume that the first link is somehow established, then only the propagation threshold needs to be exceeded to develop the flash.
Ec.init(2.16e+5,2.16e+5,-2.16e+5, z_gnd,d,N,1);
// _V/m Initiation-, Positive channel propagation-, Negative channel propagation-threshold
Vd.init(0.21e+5,-0.21e+5, z_gnd,d,N,1);
// _V/m Voltage drop in positive, negative channels
// _V/m Voltage drop in positive, negative channels
//cout << Vd.positive << endl;
//exit(112);

InitiationPoint.init((int)round(InitX/d.x), (int)round(InitY/d.y),(int)round(InitZ/d.z));
// Initiation point with coordinates expressed as i,j,k and not x,y,z

// TEST WITH EQUIPOTENTIAL SURFACE //
Vo = 1e6;	Xp = L.x/2;	Yp = L.y/2+y_shift;	Zp = L.z;	Rp1 = 2.5e+3; Rp2=Rp1; Rp3 = Rp1;
Potential V(Vo, Xp,Yp,Zp, Rp1, d,N);

// GJ //
InitX	= L.x/2;					// _m X-coordinate of initiation point
InitY	= L.y/2;					// _m Y-coordinate of initiation point
InitZ	= L.z-Rp1;					// _m Z-coordinate of initiation point


ThresholdOvershoot		= -100;		// % by which the E-field initiation threshold must be exceeded
// arbitrarily fixed to 1 in this case.
BCtype					= 0;		// Boundary conditions:
// = 0 ``Tin Can"
// = 1 Open BC
// = 2 Moving Capacitor Plates - Greifinger and Greifinger boundary
// = 3 Free Space
InitiationType			= 2;		// Initiation type:
// = 1 randomly in available region
// = 2 at point of maximum E-field magnitude
// = 3 at implemented position
rho3dCalculationStep	= 10;		// Charge density is calculated and store every rho3dCalculationStep steps.
// = 0 3-D charge density is never calculated
AddNew					= true;		// Channel is allowed to propagate: Y/N
isBCerrorCalculated		= true;		// Error at the boundary is calculated at each step: Y/N
isBndXingPossible		= false;	// Channel is allowed to cross boundaries: Y/N
isChannelEquipotential	= true;		// Channel potential is adapted to ensure charge neutrality: Y/N
isEsEnergyCalculated    = true;		// Electrostatic energy is calculated at each step: Y/N
isFlashAccoutedInBC 	= false;		// Channel charge is accounted for in derivation of BC: Y/N
isInitiationPossible	= false;	// Initiation is possible in the simulation domain after charge load: Y/N
isLinkXingPossible		= false;	// Channels crosses are allowed: Y/N
isRsDeveloped			= false;	// Return stroke development: Y/N

isInitiationPrevented	= false; 	// Only simulate cloud electrical structure

cout<<"/*******************************************************************/"<<endl;
cout<<"/* Loading charge layers                                           */"<<endl;
cout<<"/*******************************************************************/"<<endl;

ThresholdOvershoot /= 100;			// Convert % into decimal

//		C.reset(d,N);
//		for(it=ChargeCfg.begin() ; it!=ChargeCfg.end() ; it++)
//			C += *it;

//SOR.init(phi,epsilon, MaxStep, d, N, C, Un);
SOR.init(phi,epsilon, MaxStep, d, N, V, Un);
//ApplyBC(BCtype,phi,C.rho,d,N);
SOR.Solve(d,N,Un,phi);
/*
for(it=ChargeCfg.begin() ; it!=ChargeCfg.end() ; it++)
{
Q  = it->getParams()[0];
Zq = it->getParams()[3];
Rq1= it->getParams()[4];
Rq2= it->getParams()[5];
Rq3= it->getParams()[6];
Q /= (Rq1*Rq2*Rq3);				// refers to the charge density and no longer the charge content //
if (Q!=0)
for(int ii=0 ; ii<N.x ; ii++) for(int jj=0 ; jj<N.y ; jj++) for(int kk=0 ; kk<N.z ; kk++)
{
if(    (kk*d.z-Zq) > Rq3/2) phi[ii][jj][kk]+= Q/PMC.eps0*Rq3/2*   (kk*d.z - Zq - Rq3/4);
if(fabs(kk*d.z-Zq) <=Rq3/2) phi[ii][jj][kk]+= Q/PMC.eps0    /2*pow(kk*d.z - Zq      ,2);
if(    (kk*d.z-Zq) <-Rq3/2) phi[ii][jj][kk]+=-Q/PMC.eps0*Rq3/2*   (kk*d.z - Zq + Rq3/4);
};
}

// Modify propagation thresholds //
for(int kk=1 ; kk<N.z ; kk++)
{
Ec.initiation[kk] =	Ec.initiation[kk-1];
Ec.positive[kk]   =	Ec.positive[kk-1];
Ec.negative[kk]   =	Ec.negative[kk-1];
}
*/
// Derive Electrostatic energy before the discharge //
Eps_bf=0;
for(int ii=0 ; ii<N.x ; ii++) for(int jj=0 ; jj<N.y ; jj++) for(int kk=0 ; kk<N.z ; kk++)
Eps_bf += PMC.eps0*pow(Eijk(ii,jj,kk,phi,d,N)[0],2)/2*d.x*d.y*d.z;
cout<<"\ndone\n"<<endl;

cout<<"/*******************************************************************/"<<endl;
cout<<"/* Searching for Potential Extrema                                 */"<<endl;
cout<<"/*******************************************************************/"<<endl;
rho		= Globalrho(phi,d,N);
Qtot_bf = TotalCharge(rho,d,N);
phi.MinMax(Vmin,Vmax);
rho.MinMax(rhoAmbMin,rhoAmbMax);
phi_amb=phi;
cout<<"\nExtrema of Ambient Potential\n";
cout<<"phi_min = "<<setw(12)<<Vmin<<" ; phi_max = "<<setw(12)<<Vmax<<endl;
cout<<"\ndone\n"<<endl;

cout<<"/*******************************************************************/"<<endl;
cout<<"/* Initializing the tree                                           */"<<endl;
cout<<"/*******************************************************************/"<<endl;
isInitiationPossible = InitiateTree(InitiationType,
InitiationPoint,
phi0, ThresholdOvershoot,
Ec, EzNum, phiNum,
SOR, phi, Un,
d,N);
if(isInitiationPrevented == true)
{
isInitiationPossible = false;
cout<<"Initiation not allowed."<<endl;
}
cout<<"\ndone\n"<<endl;

cout<<"/*******************************************************************/"<<endl;
cout<<"/* Growing the tree                                                */"<<endl;
cout<<"/*******************************************************************/"<<endl;
if(isInitiationPossible)
GrowTree(AddNew, isBndXingPossible, isFlashAccoutedInBC, isRsDeveloped,
isLinkXingPossible,  isChannelEquipotential,
isBCerrorCalculated, isEsEnergyCalculated,
BCtype, rho3dCalculationStep,
phi0, Vmin, Vmax,
rhoAmbMin, rhoAmbMax,
QchannelPlus, QchannelMinus, Qtot_af,
Ec, Vd, ThresholdOvershoot, InitR,
InitiationPoint, EstablishedLinks,
BndUpdateErrors, DischargeDipoleMoment,
CarriedCharge, ChannelPotential, EsEnergy,
TotalPotential, TotalEfield,
SOR, epsilon, MaxStep,
phi_cha, phi_amb, phi, Un,
L, d, N, ChargeCfg);

// Calculate eletrostatic energy after the discharge //
Eps_af=0;
for(int ii=0 ; ii<N.x ; ii++) for(int jj=0 ; jj<N.y ; jj++) for(int kk=0 ; kk<N.z ; kk++)
Eps_af += PMC.eps0*pow(Eijk(ii,jj,kk,phi,d,N)[0],2)/2*d.x*d.y*d.z;
eFlux = eFieldFlux(phi,d,N);
cout<<"\ndone\n"<<endl;

cout<<"/*******************************************************************/"<<endl;
cout<<"/* Estimating fields & Storing the results                         */"<<endl;
cout<<"/*******************************************************************/"<<endl;

StoreData(isBCerrorCalculated, isEsEnergyCalculated,
InitiationPoint, InitR, ChargeCfg,
EstablishedLinks, BndUpdateErrors, DischargeDipoleMoment,
CarriedCharge, ChannelPotential, EsEnergy,
TotalPotential, TotalEfield,
Ec, EzNum, phiNum, E, phi, phi_amb,
d, N);
cout<<"\ndone\n"<<endl;

cout<<"/*******************************************************************/"<<endl;
cout<<"/* Results Summary                                                 */"<<endl;
cout<<"/*******************************************************************/"<<endl;

if(BCtype == 0)
cout<<"Closed Boundary Conditions"<<endl;
else if(BCtype == 1)
cout<<"Open Boundary Conditions"<<endl;
else if(BCtype == 2)
cout<<"Moving Capacitor Plate Boundary Conditions"<<endl;
else if(BCtype == 3)
cout<<"Free Space Boundary Conditions"<<endl;
if(isFlashAccoutedInBC==true)
cout<<"Updated Boundary Conditions"<<endl;
else
cout<<"Fast Boundary Conditions (no update)"<<endl;
cout<<"Threshold Overshoot = "<<ThresholdOvershoot*100<<" %"<<endl;

cout<<"N = ["<<N.x<<" "<<N.y<<" "<<N.z<<"]\n";
cout<<"d = ["<<d.x<<" "<<d.y<<" "<<d.z<<"]\n"<<endl;
cout.precision(10);
cout<<"El.Stat. Energy before channel propagation    : "<<Eps_bf<<" J\n";
cout<<"El.Stat. Energy after  channel propagation    : "<<Eps_af<<" J\n"<<endl;
cout<<"Total    Charge before channel propagation    : "<<Qtot_bf<<" C\n";
cout<<"Total    Charge after  channel propagation    : "<<Qtot_af<<" C\n"<<endl;
cout<<"Positive Charge in the channel                : "<<QchannelPlus<<" C\n";
cout<<"Negative Charge in the channel                : "<<QchannelMinus<<" C\n"<<endl;
cout<<"Flux of electric field through the boundaries : "<<eFlux;//[0]<<"\n";
cout<<"Equivalent charge                             : "<<Qtot_af<<" C\n";
//	cout<<"Relative error                                : "<<fabs((eFlux[0]-Qtot_af)/Qtot_af)*100<<"%\n"<<endl;
cout<<"Charge layers: "<<endl;
for(it=ChargeCfg.begin() ; it!=ChargeCfg.end() ; it++) cout<<*it;

/**************************************************************************
double		tmp;
CMatrix3D	rho_amb;
CMatrix3D	rho_cha;

QchannelPlus  = 0;
QchannelMinus = 0;
tmp			  = 0;
rho_cha = Globalrho(phi_cha,d,N);
rho_amb = Globalrho(phi_amb,d,N);
rho     = Globalrho(phi    ,d,N);
for(int ii=1 ; ii<N.x-1 ; ii++) for(int jj=1 ; jj<N.y-1 ; jj++) for(int kk=1 ; kk<N.z-1 ; kk++)
{
if(Un[ii][jj][kk] == 1)
{
tmp = rho(ii,jj,kk)-rho_amb(ii,jj,kk);
//		tmp = rho_cha(ii,jj,kk);
if(tmp>=0)
QchannelPlus  += tmp*d.x*d.y*d.z;
if(tmp<=0)
QchannelMinus += tmp*d.x*d.y*d.z;
}
};
cout<<"Positive Charge in the channel (tmp)          : "<<QchannelPlus<<" C\n";
cout<<"Negative Charge in the channel (tmp)          : "<<QchannelMinus<<" C\n"<<endl;
**************************************************************************/


/*
double Qtot_diff = TotalCharge(Globalrho(phi_cha,d,N),d,N)-ChannelCharge(Globalrho(phi_cha,d,N),Un,d,N);
double deltaQ2 = fabs(Qtot_diff);
double deltaQ1 = fabs(Qtot_bf-Qtot_af);
cout<<"Qdiff = "<<Qtot_diff<<endl;
cout<<"dQ1 = "<<setw(12)<<deltaQ1<<"; dQ2 = "<<setw(12)<<deltaQ2<<endl;
cout<<"Noise in rho_cha represents "<<(1-fabs(deltaQ2-deltaQ1)/fabs(deltaQ1))*100<<"% of the difference between Q_af and Q_bf"<<endl<<endl;
*/

clock_t endTime = clock();
clock_t runTime = endTime - startTime;
printf("Run time: %fs\n",(double)runTime/CLOCKS_PER_SEC);
//    printf("Run time: %fs\n",(double)runTime/100);
cout<<"\ndone\n";

return 0;
};
/**************************************************************************************/

/**************************************************************************************/
/* Load a Tripolar Structure for the ThunderCloud									  */
/**************************************************************************************/
void LoadTripoleModel(ListCharge& CChargeCfg,
const double II1, const double II2, const double IIscreen,
SorSolution SSOR, const double eepsilon, const int MMaxStep,
const int BBCtype, CriticalFields EEc,
const double TThresholdOvershoot,
CMatrix3D& pphi, CMatrix3D& UUn,
const StepsSizes dd, const BoxSteps NN)
{
clock_t	sstartTime	= clock();
clock_t	eendTime	= clock();
clock_t	rrunTime	= clock();
ListCharge::iterator  it;	// Table with all parameters of the charge configuation

CriticalFields OOverShotEc((1+TThresholdOvershoot)*EEc.getParams()[0],
EEc.getParams()[1],EEc.getParams()[2],EEc.getParams()[3],
dd,NN, 1);
int		cchoice		= 2;	// Define choice of method to derive bd time
double	TT			= 0;	// Bd Time
double	rr			= 0;	// E/[(1+ThresholdOvershoot)*Einit]
int		kk			= 0;	// Nb of iterations

if(cchoice == 0) // Bisection Method
{
double	TTl		= 0;				// _s
double	TTr		= 180;				// _s
double	TTav	= (TTl+TTr)/2;		// _s
double	rr_tmp	= 0;				// return value of E/[(1+ThresholdOvershoot)*Einit]
int		rr_l	= 0;				// =1 if rr_tmp>1 (ie E>Einit) ; =-1 else
int		rr_r	= 0;				// =1 if rr_tmp>1 (ie E>Einit) ; =-1 else
int		rr_av	= 0;				// =1 if rr_tmp>1 (ie E>Einit) ; =-1 else

rr_tmp	= EoverEk( TTl, CChargeCfg, II1,II2,IIscreen, SSOR,eepsilon,MMaxStep, BBCtype, OOverShotEc, pphi,UUn, dd,NN);
rr_l	= (rr_tmp>=1)*1+(rr_tmp<1)*-1;
rr_tmp	= EoverEk( TTr, CChargeCfg, II1,II2,IIscreen, SSOR,eepsilon,MMaxStep, BBCtype, OOverShotEc, pphi,UUn, dd,NN);
rr_r	= (rr_tmp>=1)*1+(rr_tmp<1)*-1;

if(TTl > TTr) {SwitchValues(TTl,TTr); SwitchValues(rr_l,rr_l);}
while(fabs(2*(TTr-TTl)/(TTr+TTl))>eepsilon)
//	while(rr_r>1e-5)
{
kk++;
TTav	= (TTl+TTr)/2;
rr_tmp	= EoverEk( TTav, CChargeCfg, II1,II2,IIscreen, SSOR,eepsilon,MMaxStep, BBCtype, OOverShotEc, pphi,UUn, dd,NN);
rr_av	= (rr_tmp>=1)*1+(rr_tmp<1)*-1;
if(rr_l*rr_av > 0) { TTl = TTav; rr_l = rr_av; }
else {TTr = TTav; rr_r = rr_av;}
}
TT = (TTl+TTr)/2;
}
else if(cchoice == 1)	//Adaptive timestep
{
double	dt		= 60;
double	rr_tmp	= 0;

while(dt>eepsilon)
{
kk++;
rr_tmp	=  EoverEk( TT+dt, CChargeCfg, II1,II2,IIscreen, SSOR,eepsilon,MMaxStep, BBCtype, OOverShotEc, pphi,UUn, dd,NN);
rr		= (rr_tmp>=1)*1+(rr_tmp<1)*-1;
while(rr>=0 && dt>eepsilon)
{
dt /=2;
rr_tmp	=  EoverEk( TT+dt, CChargeCfg, II1,II2,IIscreen, SSOR,eepsilon,MMaxStep, BBCtype, OOverShotEc, pphi,UUn, dd,NN);
rr		= (rr_tmp>=1)*1+(rr_tmp<1)*-1;
}
TT +=dt;
//	cout<<"[TT dt] = ["<<setw(8)<<TT<<setw(13)<<dt<<"]"<<endl;
}
}
else if(cchoice == 2)	// Linearity
{
double rr = 1;
TT	= 1; // _s
rr = EoverEk(TT, CChargeCfg, II1,II2,IIscreen, SSOR,eepsilon,MMaxStep, BBCtype, OOverShotEc, pphi,UUn, dd,NN);
TT *= 1/rr;
}
/**********************************************************************************/
/* Estimate how much the initiation threshold is exceeded						  */
/**********************************************************************************/
rr = EoverEk(TT, CChargeCfg, II1,II2,IIscreen, SSOR,eepsilon,MMaxStep, BBCtype, OOverShotEc, pphi,UUn, dd,NN);

/**********************************************************************************/
/* display results																  */
/**********************************************************************************/
cout<<"T = "<<setw(12)<<TT<<" ; max(E/((1+alpha)*Ek)) = "<<setw(13)<<rr<<endl;
cout<<"   Nb of iterations    = "<<kk+1<<endl;;
cout<<"Initiation Threshold is exceeded by "<<(rr*(1+TThresholdOvershoot)-1)*100<<"%\n";
cout<<"Final values of the charge domains.\n";
cout<<"   [Q1 Q2 Q3 Q4]       = [ ";
for (it=CChargeCfg.begin() ; it!=CChargeCfg.end() ; it++)
cout<<it->getParams()[0]<<" ";
cout<<"] C\n";
eendTime = clock();
rrunTime = eendTime - sstartTime;
printf("\nRun time for layers load   : %fs\n",(double)rrunTime/CLOCKS_PER_SEC);
//    printf("\nRun time for layers load   : %fs\n",(double)rrunTime/100);
};

double EoverEk(const double tt,
ListCharge& CChargeCfg,
const double II1, const double II2, const double IIscreen,
SorSolution SSOR, const double eepsilon, const int MMaxStep,
const int BBCtype, CriticalFields OOverShotEc,
CMatrix3D& pphi, CMatrix3D& UUn,
const StepsSizes dd, const BoxSteps NN)
{
double					MMaxRatio(0);
double					rratio(0);
int						nn;				// Counter of iteration
double					QQ(0.), XXq(0.),YYq(0.),ZZq(0.), RRq1(0.),RRq2(0.),RRq3(0.);
Charge					CC(dd,NN);		// Temporary charge structure for SOR derivation
ListCharge::iterator	it;				// Table with all parameters of the charge configuation

nn=0;
for (it=CChargeCfg.begin() ; it!=CChargeCfg.end() ; it++)
{
/* retrieve tripole geometrical parameters */
XXq = it->getParams()[1];
YYq = it->getParams()[2];
ZZq = it->getParams()[3];
RRq1= it->getParams()[4];
RRq2= it->getParams()[5];
RRq3= it->getParams()[6];

/* evaluate layers charge */
if(nn==0) QQ = -II2*tt;
if(nn==1) QQ = (-II1+II2)*tt;
if(nn==2) QQ = (II1*tt);
if(nn==3) QQ = -IIscreen*tt - ( -II2*tt + (-II1+II2)*tt + (II1*tt));
if(nn>=4) QQ = 0;

/* update layers */
it->disk(QQ, XXq,YYq,ZZq, RRq1,RRq3, dd,NN);
CC+=*it;
nn++;
}

/* solve fields */
SSOR.init(pphi,eepsilon, MMaxStep, dd, NN, CC, UUn);
ApplyBC(BBCtype,pphi,CC.rho,dd,NN);
SSOR.Solve(dd,NN,UUn,pphi);

/* Find ddeltaE */
for(int kk=0; kk<NN.z; kk++) for(int jj=1; jj<NN.y-1; jj++) for(int ii=0; ii<NN.x; ii++)
{
rratio = Eijk(ii,jj,kk,pphi,dd,NN)[0]/OOverShotEc.initiation[kk];
if(rratio>=MMaxRatio)	MMaxRatio	= rratio;
};
return MMaxRatio;
}
/**************************************************************************************/

/**************************************************************************************/
/* Initiate Tree (Find initiation point)											  */
/**************************************************************************************/
bool InitiateTree(const int IInitiationType,
Point& IInitiationPoint,
double& pphi0, const double TThresholdOvershoot,
CriticalFields EEc, CMatrix1D& EEzNum, CMatrix1D& pphiNum,
SorSolution SSOR,
CMatrix3D& pphi, CMatrix3D& UUn,
const StepsSizes dd, const BoxSteps NN)
{
bool	isInitiated = false;

if(IInitiationType == 1)
{
list<Point>				LListOfCandidates;
list<Point> :: iterator iit;
Point					CCandidate;
int						ccptPoints(0);
int						RRandomPoint;
int						kk=0;

for(int ii=1; ii<NN.x-1; ii++) for(int jj=1; jj<NN.y-1; jj++) for(int kk=1; kk<NN.z-1; kk++)
{
/* Choose the points which are likely to launch a flash	*/
if(Eijk(ii,jj,kk, pphi,dd,NN)[0] >= (1+0*TThresholdOvershoot)*EEc.initiation[kk])
{
CCandidate.i=ii;
CCandidate.j=jj;
CCandidate.k=kk;
LListOfCandidates.push_back(CCandidate);
ccptPoints++;
}
};

if(ccptPoints>0)
{
/* Pick a point randomly in this list */
RRandomPoint = rand()%ccptPoints+1;
for(iit = LListOfCandidates.begin() ; iit!=LListOfCandidates.end() ; iit++)
{
if(kk==RRandomPoint)
{
IInitiationPoint = *iit;
UUn[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k] = 1;
pphi0 = pphi[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k];
break;
}
kk++;
}
isInitiated = true;
}
else
{
isInitiated = false;
cout<<"Cannot initiate flash (unsufficient field)\n";
}
}
else if(IInitiationType == 2)
{
/******************/
/* Init at Emax   */
/******************/
double Emax = 0;
for(int ii=1; ii<NN.x-1; ii++) for(int jj=1; jj<NN.y-1; jj++) for(int kk=1; kk<NN.z-1; kk++)
if(Eijk(ii,jj,kk, pphi,dd,NN)[0] >= EEc.initiation[kk])
if(Eijk(ii,jj,kk, pphi,dd,NN)[0] >= Emax)
{
Emax=Eijk(ii,jj,kk, pphi,dd,NN)[0];
IInitiationPoint.init(ii,jj,kk);
};
UUn[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k] = 1;
pphi0 = pphi[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k];
if(Eijk(IInitiationPoint.i,IInitiationPoint.j,IInitiationPoint.k, pphi,dd,NN)[0] < EEc.initiation[IInitiationPoint.k])
{
isInitiated = false;
cout<<"Cannot initiate flash (unsufficient field)\n";
}
else {isInitiated = true;}
}
else if(IInitiationType == 3)
{
/******************/
/* Fix Init Point */
/******************/
if(Eijk(IInitiationPoint.i, IInitiationPoint.j, IInitiationPoint.k, pphi,dd,NN)[0] >= EEc.initiation[IInitiationPoint.k])
{
UUn[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k] = 1;
pphi0 = pphi[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k];
isInitiated = true;
}
else
{
cout<<"Cannot initiate flash (unsufficient field)\n";
isInitiated = false;
}
}
else
{
cout<<"Wrong Choice for initiation procedure"<<endl;
exit(3);
}
/* Display coordinates of initiation point */
cout<<"Initx = "<<IInitiationPoint.i*dd.x*1e-3<<" km"<<endl;
cout<<"Inity = "<<IInitiationPoint.j*dd.y*1e-3<<" km"<<endl;
cout<<"Initz = "<<IInitiationPoint.k*dd.z*1e-3<<" km"<<endl;
//	SSOR.Solve(dd,NN,UUn,pphi);

/* Derive field and potential before discharge along the central vertical axis */
for(int kk=0 ; kk<NN.z ; kk++)
{
EEzNum[kk]	= Eijk((NN.x-1)/2,(NN.y-1)/2,kk,pphi,dd,NN)[3];
pphiNum[kk]	= pphi((NN.x-1)/2,(NN.y-1)/2,kk);
}

/* Store Initiation data */
pphiNum.fwrite("results/phiNumBF.dat");
EEzNum.fwrite("results/EnumBF.dat");
return isInitiated;
}
/**************************************************************************************/

/**************************************************************************************/
/* Grow the Tree																	  */
/**************************************************************************************/
void GrowTree(bool AAddNew, const bool iisBndXingPossible, const bool iisFlashAccoutedInBC,
const bool iisRsDeveloped,
const bool iisLinkXingPossible,  const bool iisChannelEquipotential,
const bool iisBCerrorCalculated, const bool iisEsEnergyCalculated,
const int BBCtype, const int rrho3dCalculationStep,
double& pphi0, double VVmin, double VVmax,
const double rrhoAmbMin, const double rrhoAmbMax,
double& QQchannelPlus, double& QQchannelMinus, double& QQtot_af,
CriticalFields& EEc, VoltageDrops& VVd,
const double TThresholdOvershoot, const double IInitR,
const Point& IInitiationPoint, ListLink& EEstablishedLinks,
ListVector& BBndUpdateErrors, ListVector& DDischargeDipoleMoment,
ListDouble& CCarriedCharge, ListDouble& CChannelPotential, ListDouble& EEsEnergy,
ListCMatrix1D& TTotalPotential, ListCMatrix1D& TTotalEfield,
SorSolution SSOR, const double eepsilon, const int MMaxStep,
CMatrix3D& pphi_cha, CMatrix3D& pphi_amb,
CMatrix3D& pphi, CMatrix3D& UUn,
const BoxLengths& LL, const StepsSizes dd, const BoxSteps NN, ListCharge& CChargeCfg)
{
bool					isConnectedToTheGround;			// Check connection to the ground
char					nname3D[50];					// File name for storage of charge density
int						ccptLinks(0);					// Current link iteration
double					BBndErrorTmp = 0;				// Error at the boundaries at the current step
double					EEsEnergyTmp = 0;				// Electrostatic energy at the current step
string					ffname3D;						// File name for storage of charge densitycurrent step
CMatrix1D				rrho3DStep(1);					// Variable used for storage of rrho3dCalculationStep
CMatrix1D				TTotalPotentialTemp(NN.z);		// Total potential on the central vertical axis at the current step
CMatrix1D				TTotalEfieldTemp(NN.z);			// Total eField on the central vertical axis at the current step
CMatrix2D				pphi2D_cha(NN.y,NN.z);			// Potential induced by the channel in the y-z plane
CMatrix3D				pphi_chaTmp;					// Potential due to the channel in the domain before update of the boundary conditions
CMatrix3D				rrho_cha(NN.x,NN.y,NN.z);		// 3-D Matrix for channel induced charge density
CMatrix1D				rrho3D(NN.x*NN.y*NN.z);			// Total charge density everywhere at the current step
int						nn = 0;							// Blind variable
ListLink::iterator		LLit;							// index on the list of established links
ListVector::iterator	it;								// index on a list of vectors
Charge					CC;								// Charge transfer at the current step
Potential				PP;								// Channel potential at the current step
Vector					pp;								// x-, y-, z-components and norm of the dipole moment at the current step
Vector					BBndError;						// Potential at the boundary before update, after update and relative difference between the two values

nn = 0;
for(int kk=0 ; kk<NN.z ; kk++)
{
TTotalEfieldTemp[kk]	= Eijk((NN.x-1)/2,(NN.y-1)/2,kk,pphi,dd,NN)[3];
TTotalPotentialTemp[kk]	= pphi((NN.x-1)/2,(NN.y-1)/2,kk);
if(iisEsEnergyCalculated == true && rrho3dCalculationStep == 0)
for(int jj=0 ; jj<NN.y ; jj++)
for(int ii=0 ; ii<NN.x ; ii++)
EEsEnergyTmp += PMC.eps0*pow(Eijk(ii,jj,kk,pphi,dd,NN)[0],2)/2*dd.x*dd.y*dd.z;
if(iisEsEnergyCalculated == false && rrho3dCalculationStep != 0)
for(int jj=0 ; jj<NN.y ; jj++) for(int ii=0 ; ii<NN.x ; ii++)
{
rrho3D[nn]		= rhoijk(ii,jj,kk,pphi,dd,NN)*1e+9; //_nC
nn++;
};
if(iisEsEnergyCalculated == true && rrho3dCalculationStep != 0)
for(int jj=0 ; jj<NN.y ; jj++) for(int ii=0 ; ii<NN.x ; ii++)
{
EEsEnergyTmp += PMC.eps0*pow(Eijk(ii,jj,kk,pphi,dd,NN)[0],2)/2*dd.x*dd.y*dd.z;
rrho3D[nn] = rhoijk(ii,jj,kk,pphi,dd,NN)*1e+9; //_nC
nn++;
};
}
if(rrho3dCalculationStep != 0)
{
rrho3DStep[0] = rrho3dCalculationStep;
rrho3D.fwrite("results/rho3d0.dat");
rrho3DStep.fwrite("results/rho3dStep.dat");

StoreIntermData(IInitiationPoint, IInitR, CChargeCfg,
EEstablishedLinks, EEc, pphi, pphi_amb,
dd, NN, ccptLinks);

}

pp = DipoleMoment(QQchannelPlus,pphi_cha,UUn,LL,dd,NN);
CChannelPotential.push_back(pphi0);
DDischargeDipoleMoment.push_back(pp);
CCarriedCharge.push_back(QQchannelPlus);
EEsEnergy.push_back(EEsEnergyTmp);
TTotalEfield.push_back(TTotalEfieldTemp);
TTotalPotential.push_back(TTotalPotentialTemp);

while(AAddNew==true && ccptLinks>=0)
{
ccptLinks++;
AAddNew	= AddNewLink(dd,NN,UUn,pphi, EEc,VVd, IInitiationPoint,EEstablishedLinks,
iisBndXingPossible,  iisRsDeveloped,
iisLinkXingPossible, iisChannelEquipotential);

/* Check Connection to the ground */
isConnectedToTheGround = false;
for(LLit = EEstablishedLinks.begin(); LLit != EEstablishedLinks.end() ; LLit++)
if(LLit->end.k == 0)
{
isConnectedToTheGround = true;
break;
}
if(iisChannelEquipotential==false)
{
UpdateBC(iisFlashAccoutedInBC,BBCtype,pphi,rrhoAmbMin,rrhoAmbMax,dd,NN);
SSOR.Solve(dd,NN,UUn,pphi);
};
if(iisChannelEquipotential==true)
{
if(isConnectedToTheGround == false)
pphi0 = fMinSearch(pphi0, QQchannelPlus, VVmin,VVmax , eepsilon,MMaxStep, pphi_cha,pphi_amb ,UUn, IInitiationPoint,EEstablishedLinks, dd,NN);
if (isConnectedToTheGround == true)
pphi0 = EqualizeAtGroundPotential(eepsilon,MMaxStep, pphi_cha,pphi_amb ,UUn, IInitiationPoint,EEstablishedLinks, dd,NN);
CChannelPotential.push_back(pphi0);

if(iisBCerrorCalculated == false)
{
UpdateBC(iisFlashAccoutedInBC,BBCtype,pphi_cha,rrhoAmbMin,rrhoAmbMax,dd,NN);
pphi			= pphi_amb+pphi_cha;
}

if(iisBCerrorCalculated ==true)
{
pphi_chaTmp		= pphi_cha;
UpdateBC(iisFlashAccoutedInBC,BBCtype,pphi_cha,rrhoAmbMin,rrhoAmbMax,dd,NN);
pphi			= pphi_amb+pphi_cha;

BBndError.x = 0;
BBndError.y = 0;
BBndError.z = 0;
for(int kk=0; kk<NN.z; kk++) for(int jj=0; jj<NN.y; jj++) for(int ii=0; ii<NN.x; ii++)
if( ii == 0 || ii == NN.x-1 ||jj == 0 || jj == NN.y-1 || kk == 0 || kk == NN.z-1)
{
BBndErrorTmp = fabs(pphi_cha(ii,jj,kk) - pphi_chaTmp(ii,jj,kk));
if(BBndErrorTmp >= fabs(BBndError.x-BBndError.y))
{
BBndError.x = pphi_chaTmp(ii,jj,kk);
BBndError.y = pphi_cha(ii,jj,kk);
BBndError.z = pphi_amb(ii,jj,kk);
}
}
BBndUpdateErrors.push_back(BBndError);
}

/**************************************************************************/
/* The update of the BC introduce sharp gradient and artificial charge    */
/* close to the boundary, which slightly violate the charge conservation. */
/* Because potential of the channel at stage N is influenced by the BC due*/
/* to the channel at step N-1. It is not possible to completely compensate*/
/* this effect. It is not possible to derive the potential of the channel */
/* and the potential at the boundary at the same time.					  */
/* We can try an approximation, using the potential of the channel derived*/
/* at step N, using BC at N-1. We then derive new BC, and re-run an SOR   */
/* algorithm to minimize the error.										  */
/**************************************************************************/
// Solution //
/*
PP.init(pphi_cha,UUn);
SSOR.init(pphi_cha, eepsilon,MMaxStep, dd, NN, PP, UUn);
SSOR.Solve(dd,NN,UUn,pphi_cha);
pphi = pphi_amb+pphi_cha;
*/

EEsEnergyTmp = 0;
if(rrho3dCalculationStep != 0 && ccptLinks%rrho3dCalculationStep==0)
{
sprintf(nname3D,"results/rho3d%d.dat", ccptLinks);
ffname3D	= nname3D;
nn          = 0;

StoreIntermData(IInitiationPoint, IInitR, CChargeCfg,
EEstablishedLinks, EEc, pphi, pphi_amb,
dd, NN, ccptLinks);

}
for(int kk=0 ; kk<NN.z ; kk++)
{
TTotalEfieldTemp[kk]	= Eijk((NN.x-1)/2,(NN.y-1)/2,kk,pphi,dd,NN)[3];
TTotalPotentialTemp[kk]	= pphi((NN.x-1)/2,(NN.y-1)/2,kk);
if(iisEsEnergyCalculated == true && rrho3dCalculationStep == 0)
for(int jj=0 ; jj<NN.y ; jj++) for(int ii=0 ; ii<NN.x ; ii++)
EEsEnergyTmp += PMC.eps0*pow(Eijk(ii,jj,kk,pphi,dd,NN)[0],2)/2*dd.x*dd.y*dd.z;
if(iisEsEnergyCalculated == false && rrho3dCalculationStep != 0 && ccptLinks%rrho3dCalculationStep==0)
for(int jj=0 ; jj<NN.y ; jj++) for(int ii=0 ; ii<NN.x ; ii++)
{
rrho3D[nn] = rhoijk(ii,jj,kk,pphi,dd,NN)*1e+9; //_nC
nn++;
};
if(iisEsEnergyCalculated == true && rrho3dCalculationStep != 0)
for(int jj=0 ; jj<NN.y ; jj++) for(int ii=0 ; ii<NN.x ; ii++)
{
EEsEnergyTmp += PMC.eps0*pow(Eijk(ii,jj,kk,pphi,dd,NN)[0],2)/2*dd.x*dd.y*dd.z;
if(ccptLinks%rrho3dCalculationStep==0)
{
rrho3D[nn] = rhoijk(ii,jj,kk,pphi,dd,NN)*1e+9; //_nC
nn++;
}
};
}
if(rrho3dCalculationStep != 0 && ccptLinks%rrho3dCalculationStep==0) rrho3D.fwrite(ffname3D);

pp = DipoleMoment(QQchannelPlus,pphi_cha,UUn,LL,dd,NN);
DDischargeDipoleMoment.push_back(pp);
CCarriedCharge.push_back(QQchannelPlus);
EEsEnergy.push_back(EEsEnergyTmp);
TTotalEfield.push_back(TTotalEfieldTemp);
TTotalPotential.push_back(TTotalPotentialTemp);
};
cout<<"------> Nb of established links   : "<<ccptLinks<<endl;
cout<<endl;
}
cout<<"Total Number of links: "<<ccptLinks<<endl;
rrho_cha		= Globalrho(pphi_cha,dd,NN);
QQchannelPlus	= ChannelChargePositive(rrho_cha,UUn,dd,NN);
QQchannelMinus	= ChannelChargeNegative(rrho_cha,UUn,dd,NN);
QQtot_af		= TotalCharge(Globalrho(pphi,dd,NN),dd,NN);

bool FFlagInit=false;
for(int ii=1; ii<NN.x-1; ii++) for(int jj=1; jj<NN.y-1; jj++) for(int kk=1; kk<NN.z-1; kk++)
if(Eijk(ii,jj,kk, pphi,dd,NN)[0] >= (1+TThresholdOvershoot)*EEc.initiation[kk])
{
//			cout<<"Initiation Threshold exceeded at: ["<<ii*dd.x<<" "<<jj*dd.y<<" "<<kk*dd.z<<"]\n";
FFlagInit=true;
break;
};
if(FFlagInit==false) cout<<"*** No further initiation possible. ***\n";
else{cout<<"*** Further initiation possible. ***\n";}
if (iisBCerrorCalculated == true)
cout<<"Error on the boundaries\n"<<setw(12)<<"phiCha(N)"<<setw(12)<<"phiCha(N+1)"<<setw(16)<<"phiAmb"<<endl;
};
/**************************************************************************************/

/**************************************************************************************/
/* Store main data																	  */
/**************************************************************************************/
void StoreData(const bool iisBCerrorCalculated, const bool iisEsEnergyCalculated,
const Point& IInitiationPoint, const double IInitR, ListCharge& CChargeCfg,
ListLink EEstablishedLinks,
ListVector BBndUpdateErrors, ListVector DDischargeDipoleMoment,
ListDouble CCarriedCharge, ListDouble CChannelPotential, ListDouble EEsEnergy,
ListCMatrix1D& TTotalPotential, ListCMatrix1D& TTotalEfield,
CriticalFields EEc, CMatrix1D& EEzNum, CMatrix1D& pphiNum,
CMatrix3D& EE, CMatrix3D& pphi, CMatrix3D& pphi_amb,
const StepsSizes dd, const BoxSteps NN)
{
CMatrix3D				rrho(NN.x, NN.y, NN.z);
CMatrix3D				rrho_amb(NN.x, NN.y, NN.z);
CMatrix2D				rrho_amb_yz( NN.y, NN.z);
CMatrix2D				rrho_amb_xz( NN.x, NN.z);
CMatrix2D				EE2D(NN.y, NN.z);
CMatrix2D				EEx2D(NN.y, NN.z);
CMatrix2D				EEy2D(NN.y, NN.z);
CMatrix2D				EEz2D(NN.y, NN.z);
CMatrix2D				pphi2D(NN.y, NN.z);
CMatrix1D				zz_gnd(1);
ListCMatrix1D			CChargeLayers;
ListCharge::iterator	iit;

// Derive required values before storage //
EE		= GlobalE(pphi,dd,NN);
rrho	= Globalrho(pphi,dd,NN);
rrho_amb= Globalrho(pphi_amb,dd,NN);

for(int kk=0 ; kk<NN.z ; kk++)
{
EEzNum[kk]	= Eijk((NN.x-1)/2,(NN.y-1)/2,kk,pphi,dd,NN)[3];
pphiNum[kk]	= pphi((NN.x-1)/2,(NN.y-1)/2,kk);
for(int ii=0 ; ii<NN.x ; ii++)
rrho_amb_xz[ii][kk]= rhoijk(ii,(NN.y-1)/2,kk,pphi_amb,dd,NN)*1e9;

for(int jj=0 ; jj<NN.y ; jj++)
{
rrho_amb_yz[jj][kk]= rhoijk((NN.x-1)/2,jj,kk,pphi_amb,dd,NN)*1e9;
EE2D[jj][kk]	   = Eijk((NN.x-1)/2,jj,kk,pphi,dd,NN)[0];
EEx2D[jj][kk]	   = Eijk((NN.x-1)/2,jj,kk,pphi,dd,NN)[1];
EEy2D[jj][kk]	   = Eijk((NN.x-1)/2,jj,kk,pphi,dd,NN)[2];
EEz2D[jj][kk]	   = Eijk((NN.x-1)/2,jj,kk,pphi,dd,NN)[3];
pphi2D[jj][kk]	   = pphi((NN.x-1)/2,jj,kk);
for(int ii=0 ; ii<NN.x ; ii++)
{
rrho(ii,jj,kk)		*= 1e9;
rrho_amb(ii,jj,kk)	*= 1e9;
}
}
}
for (iit=CChargeCfg.begin() ; iit!=CChargeCfg.end() ; iit++)
CChargeLayers.push_back(iit->getParams());

zz_gnd[0] = EEc.getParams()[3];

// Store altitude of the ground plane //
zz_gnd.fwrite("results/z_gnd.dat");

// Store Values in main plane //
rrho.fwrite("results/rho.dat");
rrho_amb.fwrite("results/rhoAmb.dat");
rrho_amb_yz.fwrite("results/rhoAmbYZ.dat");
rrho_amb_xz.fwrite("results/rhoAmbXZ.dat");
EE2D.fwrite("results/E2D.dat");
EEx2D.fwrite("results/Ex2D.dat");
EEy2D.fwrite("results/Ey2D.dat");
EEz2D.fwrite("results/Ez2D.dat");
pphi2D.fwrite("results/phi2D.dat");

/**************************************************************************
double	QchannelPlus;
double	QchannelMinus;
double	tmp;
ListDouble LL1(read("results/rho.dat"));
ListDouble LL2(read("results/rhoAmb.dat"));
ListDouble::iterator itLD1;
ListDouble::iterator itLD2;

QchannelPlus  = 0;
QchannelMinus = 0;
tmp			  = 0;

itLD2= LL2.begin();
for (itLD1=LL1.begin() ; itLD1!=LL1.end() ; itLD1++)
{
tmp = (*itLD1 - *itLD2)*1e-9;
if(tmp>=0)
QchannelPlus  += tmp*(float)dd.x*(float)dd.y*(float)dd.z;
if(tmp<=0)
QchannelMinus += tmp*(float)dd.x*(float)dd.y*(float)dd.z;
itLD2++;
}
cout.precision(10);
cout<<"Positive Charge in the channel (rounded)      : "<<QchannelPlus<<" C\n";
cout<<"Negative Charge in the channel (rounded)      : "<<QchannelMinus<<" C\n"<<endl;

QchannelPlus  = 0;
QchannelMinus = 0;
tmp			  = 0;
for(int ii=1 ; ii<NN.x-1 ; ii++) for(int jj=1 ; jj<NN.y-1 ; jj++) for(int kk=1 ; kk<NN.z-1 ; kk++)
{
tmp = rrho(ii,jj,kk)-rrho_amb(ii,jj,kk);
if (UUn(ii,jj,kk)!=0)
{
if(tmp>=0)
QchannelPlus  += tmp*dd.x*dd.y*dd.z*1e-9;
if(tmp<=0)
QchannelMinus += tmp*dd.x*dd.y*dd.z*1e-9;
}
};
cout<<"Positive Charge in the channel                : "<<QchannelPlus<<" C\n";
cout<<"Negative Charge in the channel                : "<<QchannelMinus<<" C\n"<<endl;
**************************************************************************/

// Store channel potential, charge and dipole moment //
write(CChannelPotential,"results/ChannelPotentials.dat");
write(CCarriedCharge,"results/CarriedCharge.dat");
write(DDischargeDipoleMoment,"results/DischargeDipoleMoment.dat");

// Store values of the electrostatic energy in the domain (if calculated) //
if(iisEsEnergyCalculated == true)
write(EEsEnergy,"results/EsEnergy.dat");

// Store values of the boundary errors in the domain (if calculated) //
if(iisBCerrorCalculated == true)
write(BBndUpdateErrors,"results/BndError.dat");

// Store d.z //
CMatrix1D ddxyz(3);
ddxyz[0] = dd.x;
ddxyz[1] = dd.y;
ddxyz[2] = dd.z;
ddxyz.fwrite("results/dxyz.dat");

// Store StepSizes //
CMatrix1D NNxyz(3);
NNxyz[0] = NN.x;
NNxyz[1] = NN.y;
NNxyz[2] = NN.z;
NNxyz.fwrite("results/Nxyz.dat");

// Store Init point //
CMatrix1D IInitPoint(4);
IInitPoint[0] = IInitiationPoint.i*dd.x;
IInitPoint[1] = IInitiationPoint.j*dd.y;
IInitPoint[2] = IInitiationPoint.k*dd.z;
IInitPoint[3] = IInitR;
IInitPoint.fwrite("results/InitPoint.dat");

// Store charge center parameters //
write(CChargeLayers,"results/ChargeLayers.dat");

// Store Values of the fields and potential on the principal axis //
pphiNum.fwrite("results/phiNumAF.dat");
EEzNum.fwrite("results/EnumAF.dat");
EEc.initiation.fwrite("results/Einitiation.dat");
EEc.negative.fwrite("results/EthNegative.dat");
EEc.positive.fwrite("results/EthPositive.dat");
write(TTotalPotential,"results/TotalPotential.dat");
write(TTotalEfield,"results/TotalEfield.dat");

// Store Values of the field and potential everywhere //
pphi.fwrite("results/phi.dat");
EE.fwrite("results/E.dat");

// Store the list of Established Links //
write(EEstablishedLinks,"results/EstablishedLinks.dat");
}
/**************************************************************************************/

void StoreIntermData(const Point& IInitiationPoint, const double IInitR, ListCharge& CChargeCfg,
ListLink EEstablishedLinks,
CriticalFields EEc, CMatrix3D& pphi, CMatrix3D& pphi_amb,
const StepsSizes dd, const BoxSteps NN, int cccptLinks)
{
CMatrix3D				rrho(NN.x, NN.y, NN.z);
CMatrix3D				rrho_amb(NN.x, NN.y, NN.z);
CMatrix2D				rrho_amb_yz( NN.y, NN.z);
CMatrix2D				rrho_amb_xz( NN.x, NN.z);
CMatrix2D				EE2D(NN.y, NN.z);
CMatrix2D				EEx2D(NN.y, NN.z);
CMatrix2D				EEy2D(NN.y, NN.z);
CMatrix2D				EEz2D(NN.y, NN.z);
CMatrix2D				pphi2D(NN.y, NN.z);
CMatrix1D				zz_gnd(1);
ListCMatrix1D			CChargeLayers;
ListCharge::iterator	iit;

for(int kk=0 ; kk<NN.z ; kk++)
{
for(int ii=0 ; ii<NN.x ; ii++)
rrho_amb_xz[ii][kk]= rhoijk(ii,(NN.y-1)/2,kk,pphi_amb,dd,NN)*1e9;

for(int jj=0 ; jj<NN.y ; jj++)
{
rrho_amb_yz[jj][kk]= rhoijk((NN.x-1)/2,jj,kk,pphi_amb,dd,NN)*1e9;
EE2D[jj][kk]	   = Eijk((NN.x-1)/2,jj,kk,pphi,dd,NN)[0];
EEx2D[jj][kk]	   = Eijk((NN.x-1)/2,jj,kk,pphi,dd,NN)[1];
EEy2D[jj][kk]	   = Eijk((NN.x-1)/2,jj,kk,pphi,dd,NN)[2];
EEz2D[jj][kk]	   = Eijk((NN.x-1)/2,jj,kk,pphi,dd,NN)[3];
pphi2D[jj][kk]	   = pphi((NN.x-1)/2,jj,kk);
}
}

if (cccptLinks==0) {
for (iit=CChargeCfg.begin() ; iit!=CChargeCfg.end() ; iit++)
CChargeLayers.push_back(iit->getParams());

zz_gnd[0] = EEc.getParams()[3];

// Store altitude of the ground plane //
zz_gnd.fwrite("results/z_gnd.dat");

// Store d.z //
CMatrix1D ddxyz(3);
ddxyz[0] = dd.x;
ddxyz[1] = dd.y;
ddxyz[2] = dd.z;
ddxyz.fwrite("results/dxyz.dat");

// Store StepSizes //
CMatrix1D NNxyz(3);
NNxyz[0] = NN.x;
NNxyz[1] = NN.y;
NNxyz[2] = NN.z;
NNxyz.fwrite("results/Nxyz.dat");

// Store Init point //
CMatrix1D IInitPoint(4);
IInitPoint[0] = IInitiationPoint.i*dd.x;
IInitPoint[1] = IInitiationPoint.j*dd.y;
IInitPoint[2] = IInitiationPoint.k*dd.z;
IInitPoint[3] = IInitR;
IInitPoint.fwrite("results/InitPoint.dat");

// Store charge center parameters //
write(CChargeLayers,"results/ChargeLayers.dat");
}


char    nnametmp[50];
string  ffnametmp;

sprintf(nnametmp,"results/rhoAmbYZ%d.dat", cccptLinks);
ffnametmp	= nnametmp;
rrho_amb_yz.fwrite(ffnametmp);

sprintf(nnametmp,"results/rhoAmbXZ%d.dat", cccptLinks);
ffnametmp	= nnametmp;
rrho_amb_xz.fwrite(ffnametmp);

sprintf(nnametmp,"results/E2D%d.dat", cccptLinks);
ffnametmp	= nnametmp;
EE2D.fwrite(ffnametmp);

sprintf(nnametmp,"results/Ex2D%d.dat", cccptLinks);
ffnametmp	= nnametmp;
EEx2D.fwrite(ffnametmp);

sprintf(nnametmp,"results/Ey2D%d.dat", cccptLinks);
ffnametmp	= nnametmp;
EEy2D.fwrite(ffnametmp);

sprintf(nnametmp,"results/Ez2D%d.dat", cccptLinks);
ffnametmp	= nnametmp;
EEz2D.fwrite(ffnametmp);

sprintf(nnametmp,"results/phi2D%d.dat", cccptLinks);
ffnametmp	= nnametmp;
pphi2D.fwrite(ffnametmp);


// Store the list of Established Links //
write(EEstablishedLinks,"results/EstablishedLinks.dat");
}
