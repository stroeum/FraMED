/*
 *  Input.cpp
 *  Created by Jeremy Riousset on 10/26/07.
 */

#include "Input.h"

/*Variables Declaration************************************************************/
SizeGrid				Var::N;														// Number of discretization points
SizeDomain				Var::L;														// Dimensions of the simulation domain
ResGrid					Var::d;														// Lengths of the discretization-grid


double					Var::z_gnd;													// Altitude of ground plane
double					Var::z_shift;												// Vertical displacement of the cloud
double					Var::y_shift;												// Horizontal displacement of the upper part of the cloud

const double			Var::epsilon	= 1e-10;									// SOR precision
const int				Var::MaxStep	= 7500;										// Allowed maximum number of point per SOR iteration

double					Var::InitX;													// X-coordinate of the initiation point
double					Var::InitY;													// Y-coordinate of the initiation point
double					Var::InitZ;													// Z-coordinate of the initiation point
Point					Var::InitiationPoint;										// i,j,k coordinates of the initiation point
double					Var::InitR;													// Radius of the initiation region (if applicable)

/* SAM variable. */
int						Var::NumLinks;												// Number of established links in the current simulation run


// The following values are only use to avoid warning at compilation //
double					Var::phi0			= 0;									// Channel potential after minization
double					Var::Vmin			= 0;									// Minimum in cloud potential
double					Var::Vmax			= 0;									// Maximum in cloud potential
double					Var::MaximumEfield  = 0;									// Maximum electric field magnitude throughout simulation
double					Var::rhoAmbMin		= 0;									// Minimum in cloud charge density
double					Var::rhoAmbMax		= 0;									// Maximum in cloud charge density
double					Var::QchannelPlus	= 0;									// Positive charge content of the channel (== charge transfer)
double					Var::QchannelMinus	= 0;									// Negative charge content of the channel
double					Var::Eps_bf			= 0;									// Electrostatic energy before the discharge
double					Var::Eps_af			= 0;									// Electrostatic energy after the discharge
double					Var::Qtot_bf		= 0;									// Total charge in the simulation domain before the discharge
double					Var::Qtot_af		= 0;									// Total charge in the simulation domain after the discharge
// The previous  values are only use to avoid warning at compilation //

double					Var::I1, Var::I2, Var::Iscreen;								// Loading currents I1 and I2, and screening current Iscreen
double                  Var::Q, Var::Xq,Var::Yq,Var::Zq, Var::Rq1,Var::Rq2,Var::Rq3;// Charge center parameters:
                                                                                    // * Charge (Q),
                                                                                    // * X-,Y-,Z-coordinates of a charge center (Xq,Yq,Zq),
                                                                                    // * 1st, 2nd and 3rd geometrical parameters (Rq1, Rq2, Rq3)
double                  Var::Eo,Var::Vo, Var::Xp,Var::Yp,Var::Zp, Var::Rp1,Var::Rp2,Var::Rp3;
                                                                                    // Potential center parameters:
                                                                                    // * Reference electric potential (Vo) and field (Eo)
                                                                                    // * X-,Y-,Z-coordinates of a potential center (Xp,Yp,Zp),
                                                                                    // * 1st, 2nd and 3rd geometrical parameters (Rp1, Rp2, Rp3)
double					Var::ThresholdOvershoot;									// % by which threshold must be exceeded to initiate discharge
bcType                  Var::BCtype;												// Boundary conditions
initType				Var::InitiationType;										// Initiation type
loadType                Var::LoadingType;                                           // Cloud charge method
int						Var::step3d;                                                // 3-d values calculated and store every step3d.
																					//(0: 3-D charge density is never calculated)

bool					Var::isInitiationPrevented;									// Only simulate cloud electrical structure
bool					Var::AddNew;												// Channel is allowed to propagate: Y/N
bool					Var::isBCerrorCalculated;									// Error at the boundary is calculated at each step: Y/N
bool					Var::isBndXingPossible;										// Channel is allowed to cross boundaries: Y/N
bool					Var::isQMinimized;								            // Channel potential is adapted to ensure charge neutrality: Y/N
bool					Var::isEsEnergyCalculated;									// Electrostatic energy is calculated: Y/N
bool					Var::isFlashAccoutedInBC;									// Channel charge is accounted for in derivation of BC: Y/N
bool					Var::isInitiationPossible;									// Initiation is possible in the simulation domain aipfter charge load: Y/N
bool					Var::isLinkXingPossible;									// Channels crosses are allowed: Y/N
bool					Var::isRsDeveloped;											// Return stroke development: Y/N
bool                    Var::isNewRun;                                              // Check if the run is a new run or a resumed simulation
bool                    Var::isVoltageDropped;                                      // Is there a voltage drop to account for?

Charge					Var::C;														// Charge center
Potential               Var::V;                                                     // Potential center
SorSolution				Var::SOR;													// Parameters of SOR algorithms

ListDouble				Var::ChannelPotential;										// Channel potential at each stage of development
ListDouble				Var::CarriedCharge;											// Channel transfer at each stage of development
ListDouble				Var::TransportedRhoEnd;								        // Charge density difference at end-node for each stage of development
ListDouble				Var::TransportedRhoNeg;								        // Global sum of negative net charge densities at each stage of development
ListDouble				Var::TransportedRhoPos;								        // Global sum of positive net charge densities at each stage of development
ListDouble				Var::EsEnergy;												// Electrostatic energy in the simulation domain
ListLink				Var::EstablishedLinks;										// List of links constituting the discharge tree
ListVector				Var::BndUpdateErrors;										// List of errors at the boundaries
ListVector				Var::DischargeDipoleMoment;									// Channel dipole moment at each stage of development
ListCMatrix1D			Var::TotalPotential;										// Potential on a vertical axis at the center of the simulation domain
ListCMatrix1D			Var::TotalEfield;											// E-field on a vertical axis at the center of the simulation domain
ListCharge				Var::ChargeCfg;												// Table with all parameters of the charge configuation

/* SAM variables. */
ListInt					Var::NumberOfCandidates;									// Number of candidates in current link addition step.
ListDouble				Var::MaximumCandidateOverreach;								// The greatest amount that a candidate exceeded the electric field
																					//  propagation threshold.

//ListCharge::iterator	Var::it;													// Table with all parameters of the charge configuation

CriticalFields			Var::Ec;													// Initiation, propagation of positive and negative channels
VoltageDrops			Var::Vd;													// Voltage drops in positive and negative channels

CMatrix3D				Var::E;														// _V/m	Total electric field
CMatrix3D				Var::phi;													// _V	Total electric potential
CMatrix3D				Var::phi_cha;												// _V	Channel electric potential
CMatrix3D				Var::phi_amb;												// _V	Cloud electric potential
CMatrix3D				Var::rho;													// _C/m3	total charge density
CMatrix3D				Var::rho_amb;												// _C/m3	ambient charge density
CMatrix3D				Var::Un;													// Map of occupied grid points
CMatrix1D				Var::phiNum;												// _V	Total electric potential on a vertical axis in the center of simulation domain
CMatrix1D				Var::EzNum;													// _V/m	Total electric field on a vertical axis in the center of simulation domain
CMatrix1D				Var::eFlux(7);												// _C Flux of electrostatic field through the boundaries

FILE *					Var::pSum;													// Pointer to the summary file.

clock_t					Var::instanceStartTime;										// Time at which the current instance of the program started.
clock_t					Var::runStartTime;											// Time at which the current simulation run started.

disType					Var::curType;

double					Var::maxAlt;

void Var::ClearLists()
{
	ChannelPotential.clear();										// Channel potential at each stage of development
	CarriedCharge.clear();											// Channel transfer at each stage of development
	EsEnergy.clear();												// Electrostatic energy in the simulation domain
	EstablishedLinks.clear();										// List of links constituting the discharge tree
	BndUpdateErrors.clear();										// List of errors at the boundaries
	DischargeDipoleMoment.clear();									// Channel dipole moment at each stage of development
	TotalPotential.clear();											// Potential on a vertical axis at the center of the simulation domain
	TotalEfield.clear();											// E-field on a vertical axis at the center of the simulation domain
	ChargeCfg.clear();												// Table with all parameters of the charge configuation
	/* SAM variables. */
	NumberOfCandidates.clear();										// Number of candidates in current link addition step.
	MaximumCandidateOverreach.clear();
}
/**********************************************************************************/
