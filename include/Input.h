/*
 *  Input.h
 *  Created by Jeremy Riousset on 10/25/07.
 */

#ifndef INPUT_H
#define INPUT_H
#include "VoltageDrops.h"
#include "CriticalFields.h"
#include "Links.h"
#include "SorSolution.h"
#include <cstring>

typedef enum {PROPAGATING, INTRA_CLOUD, CLOUD_TO_GROUND, JET, HORIZONTAL} disType;
typedef enum {TIN_CAN, OPEN_BC, G_G, FREE_SPACE} bcType;
typedef enum {RANDOM, AT_EMAX, AT_PREDEF_POS, AT_REL_EMAX} initType;
typedef enum {SET_POTENTIAL, SET_CHARGES, CURRENTS, FROM_FILE} loadType;

struct	Vector {double x,y,z;};
typedef list<Vector> ListVector;
typedef list<double> ListDouble;
typedef list<int> ListInt;
typedef	list<CMatrix1D> ListCMatrix1D;

/*DEFINE CONSTANTS*****************************************************************/
#define	eta 1.0
/**********************************************************************************/


/*GLOBAL VARIABLES DECLARATION*****************************************************/
class Var
{
public:

	static	SizeGrid				N;												// Number of discretization points
	static	SizeDomain				L;												// Dimensions of the simulation domain
	static	ResGrid					d;												// Lengths of the discretization-grid

	static	double					z_gnd;											// Altitude of ground plane
	static	double					z_shift;										// Vertical displacement of the cloud
	static	double					y_shift;										// Horizontal displacement of the upper part of the cloud

	static	const double			epsilon;										// SOR precision
	static	const int				MaxStep;										// Allowed maximum number of point per SOR iteration

	static	double					InitX;											// X-coordinate of the initiation point
	static	double					InitY;											// Y-coordinate of the initiation point
	static	double					InitZ;											// Z-coordinate of the initiation point
	static	Point					InitiationPoint;								// i,j,k coordinates of the initiation point
	static	double					InitR;											// Radius of the initiation region (if applicable)

	/* SAM variable. */
	static int						NumLinks;										// Number of established links in the current simulation run

	static	double					phi0;											// Channel potential after minization
	static	double					Vmin;											// Minimum in cloud potential
	static	double					Vmax;											// Maximum in cloud potential
	static	double					MaximumEfield;									// Maximum electric field magnitude throughout simulation
	static	double					rhoAmbMin;										// Minimum in cloud charge density
	static	double					rhoAmbMax;										// Maximum in cloud charge density
	static	double					QchannelPlus;									// Positive charge content of the channel (== charge transfer)
	static	double					QchannelMinus;									// Negative charge content of the channel
	static	double					Eps_bf;											// Electrostatic energy before the discharge
	static	double					Eps_af;											// Electrostatic energy after the discharge
	static	double					Qtot_bf;										// Total charge in the simulation domain before the discharge
	static	double					Qtot_af;										// Total charge in the simulation domain after the discharge
																					// The previous  values are only use to avoid warning at compilation //

	static	double					I1, I2, Iscreen;								// Loading currents I1 and I2, and screening current Iscreen
	static	double					Q, Xq,Yq,Zq, Rq1,Rq2,Rq3;						// Charge center parameters:
																					// * Charge (Q),
																					// * X-,Y-,Z-coordinates of a charge center (Xq,Yq,Zq),
																					// * 1st, 2nd and 3rd geometrical parameters (Rq1, Rq2, Rq3)
    static  double                  Eo,Vo, Xp,Yp,Zp, Rp1,Rp2,Rp3;                   // Potential center parameters:
                                                                                    // * Reference electric potential (Vo) and field (Eo)
                                                                                    // * X-,Y-,Z-coordinates of a potential center (Xp,Yp,Zp),
                                                                                    // * 1st, 2nd and 3rd geometrical parameters (Rp1, Rp2, Rp3)
	static	double					ThresholdOvershoot;								// %-age by which threshold must be exceeded to initiate discharge
    
	static	bcType					BCtype;											// Boundary conditions
	static	initType				InitiationType;									// Initiation type
    static	loadType                LoadingType;                                    // Choose charge content, loading currents, reading from file, etc.

    static	int						step3d;                                         // Charge density is calculated and store every rho3dCalculationStep steps.
																					//(0: 3-D charge density is never calculated)

	static	bool					isInitiationPrevented;							// Only simulate cloud electrical structure
	static	bool					AddNew;											// Channel is allowed to propagate: Y/N
	static	bool					isBCerrorCalculated;							// Error at the boundary is calculated at each step: Y/N
	static	bool					isBndXingPossible;								// Channel is allowed to cross boundaries: Y/N
	static	bool					isQMinimized;							        // Channel potential is adapted to ensure charge neutrality: Y/N
	static	bool					isEsEnergyCalculated;							// Electrostatic energy is calculated: Y/N
	static	bool					isFlashAccoutedInBC;							// Channel charge is accounted for in derivation of BC: Y/N
	static	bool					isInitiationPossible;							// Initiation is possible in the simulation domain after charge load: Y/N
	static	bool					isLinkXingPossible;								// Channels crosses are allowed: Y/N
	static	bool					isRsDeveloped;									// Return stroke development: Y/N
    static  bool                    isNewRun;                                       // Check if the run is a new run or a resumed simulation
    static	bool					isVoltageDropped;							    // Is there a voltage drop to account for: Y/N

	static	Charge					C;												// Charge center
    static  Potential               V;                                              // Potential
	static	SorSolution				SOR;											// Parameters of SOR algorithms

	static	ListDouble				ChannelPotential;								// Channel potential at each stage of development
	static	ListDouble				CarriedCharge;									// Channel transfer at each stage of development
	static	ListDouble				TransportedRhoEnd;								// Charge density difference at end-node for each stage of development
	static	ListDouble				TransportedRhoNeg;								// Global sum of negative net charge densities at each stage of development
	static	ListDouble				TransportedRhoPos;								// Global sum of positive net charge densities at each stage of development
	static	ListDouble				EsEnergy;										// Electrostatic energy in the simulation domain
	static	ListLink				EstablishedLinks;								// List of links constituting the discharge tree
	static	ListVector				BndUpdateErrors;								// List of errors at the boundaries
	static	ListVector				DischargeDipoleMoment;							// Channel dipole moment at each stage of development
	static	ListCMatrix1D			TotalPotential;									// Potential on a vertical axis at the center of the simulation domain
	static	ListCMatrix1D			TotalEfield;									// E-field on a vertical axis at the center of the simulation domain
	static	ListCharge				ChargeCfg;										// Table with all parameters of the charge configuation

	/* SAM variables. */
	static ListInt					NumberOfCandidates;								// Number of candidates in current link addition step.
	static ListDouble				MaximumCandidateOverreach;						// The greatest amount that a candidate exceeded the electric field
																					//  propagation threshold.

	static	CriticalFields			Ec;												// Initiation, propagation of positive and negative channels
	static	VoltageDrops			Vd;												// Voltage drops in positive and negative channels

	static	CMatrix3D				E;												// _V/m	Total electric field
	static	CMatrix3D				phi;											// _V	Total electric potential
	static	CMatrix3D				phi_cha;										// _V	Channel electric potential
	static	CMatrix3D				phi_amb;										// _V	Cloud electric potential
	static	CMatrix3D				rho;											// _C/m3	total charge density
	static	CMatrix3D				rho_amb;										// _C/m3	Cloud charge density
	static	CMatrix3D				Un;												// Map of occupied grid points
	static	CMatrix1D				phiNum;											// _V	Total electric potential on a vertical axis in the center of simulation domain
	static	CMatrix1D				EzNum;											// _V/m	Total electric field on a vertical axis in the center of simulation domain
	static	CMatrix1D				eFlux;											// _C Flux of electrostatic field through the boundaries

	static	FILE *					pSum;											// Pointer to the summary file.

	static	clock_t					instanceStartTime;								// Time at which the current instance of the program started.
	static 	clock_t					runStartTime;									// Time at which the current simulation run started.

	static	disType					curType;
	static	double					maxAlt;

	static	void					ClearLists();
};
/**********************************************************************************/

#endif
