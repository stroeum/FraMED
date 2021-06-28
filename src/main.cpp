/* main.cpp Sept. 5, 2006, 3.45PM  */

/**************************************************************************************/
/* Include libraries and header files												  */
/**************************************************************************************/
#include "BoundaryConditions.h"
#include "Trees.h"
//#include "Input.h"

#include <stdlib.h>
/**************************************************************************************/

/**************************************************************************************/
/* Declare functions																  */
/**************************************************************************************/
double	EoverEk(const double tt, CriticalFields OOverShotEc);
void	LoadTripoleModel(const double II1, const double II2, const double IIscreen);
bool	InitiateTree(const int IInitiationType, Point& IInitiationPoint);
void	GrowTree(bool AAddNew);

/* SAM functions. */
void	StoreData(); /* Writes intermediate as well as final data to files. */
void	SetPathName(); /* Sets the subfolder to which interim results are written. */
void 	SortPathName(int numLinks); /* Sets the subfolder to which final results are written. */

/**************************************************************************************/

/**************************************************************************************/
/* Main																				  */
/**************************************************************************************/
int main()
{
	/******************************************************************************/
	/* Initialisation & Preparation of Summary file								  */
	/******************************************************************************/
//	srand((unsigned)time(NULL));
	ListCharge::iterator	it;													// Table with all parameters of the charge configuration
	clock_t					startTime;
	FILE *					file;
	char fileName[200];

	startTime = clock(); /* Time current simulation run. */
	/* Set the interim results subfolder path name.  The path name is based on the starting time
	 * of the current simulation run.
	 */
	IO::setPathName("results");
	IO::getPathName(fileName);
	strcat(fileName, "/summary.txt");
	file = fopen(fileName, "w");
	printf(			"ii: Starting new equipotential lightning discharge simulation.\n");
	fprintf(file,	"ii: Starting new equipotential lightning discharge simulation.\n");

	Var::isChargeContentAssumed	= true;											// Choose either charge content or loading currents

	if(Var::isChargeContentAssumed == false)
	{
		Var::N.init(41,41,61);														// Number of discretization points
		Var::L.init(12e+3,12e+3,12e+3);												// _m	Size of the simulation domain
		Var::d.init(Var::L,Var::N);													// _m	Sizes of grid-steps
		// Initiate matrices dimensions //
		Var::E.init(Var::N.x,Var::N.y,Var::N.z);									// _V/m	Total electric field
		Var::phi.init(Var::N.x,Var::N.y,Var::N.z);									// _V	Total electric potential
		Var::phi_cha.init(Var::N.x,Var::N.y,Var::N.z);								// _V	Channel electric potential
		Var::phi_amb.init(Var::N.x,Var::N.y,Var::N.z);								// _V	Cloud electric potential
		Var::rho.init(Var::N.x,Var::N.y,Var::N.z);									// _C/m3	total charge density
		Var::Un.init(Var::N.x,Var::N.y,Var::N.z);									// Map of occupied grid points
		Var::phiNum.init(Var::N.z);													// _V	Total electric potential on a vertical axis in the center of simulation domain
		Var::EzNum.init(Var::N.z);													// _V/m	Total electric field on a vertical axis in the center of simulation domain
		// Endof Initiate matrices dimensions //

		Var::z_gnd   = 3e3;															// _m Altitude of ground level
		Var::z_shift = 0e3;															// _m Vertical displacement of the cloud
		Var::y_shift = 0e3;															// _m Horizontal influence of a windshear

		//	We assume that the first link is somehow established, then only the propagation threshold needs to be exceeded to develop the flash.
		Var::Ec.init(216e+3,216e+3,-216e+3, Var::z_gnd,Var::d,Var::N);				// _V/m Initiation-, Positive channel propagation-, Negative channel propagation-threshold
		Var::Vd.init(0e+5,-0e+5, Var::z_gnd,Var::d,Var::N);							// _V/m Voltage drop in positive, negative channels
		// CG //
		Var::InitX	= 27e3;															// _m X-coordinate of initiation point
		Var::InitY	= 20e3;															// _m Y-coordinate of initiation point
		Var::InitZ	= 5.55e3;														// _m Z-coordinate of initiation point
		Var::InitR	= 0;															// _m Radius of the initiation region (if applicable)

		Var::InitiationPoint.init((int)round(Var::InitX/Var::d.x), (int)round(Var::InitY/Var::d.y),(int)round(Var::InitZ/Var::d.z));
																					// Initiation point with coordinates expressed as i,j,k and not x,y,z
		Var::I1			= 1.5;														// _A Loading current I1
		Var::I2			= -90e-3;													// _A Loading current I2
		Var::Iscreen	= 0*664e-3;													// _A Screening current Iscreen

		// Charge center parameters //

		// Q (C) Charge content; Xq,Yq,Zq (m) Charge center coordinate; R,H (m) Size of the charge center //
		Var::Q =    0;	Var::Xq = Var::L.x/2;	Var::Yq = Var::L.y/2-Var::y_shift/2;	Var::Zq = 2.00e+3+Var::z_shift; Var::Rq1 = 1.50e+3;	Var::Rq3 = 1.50e+3;
		Var::C.disk(Var::Q, Var::Xq,Var::Yq,Var::Zq, Var::Rq1,Var::Rq3, Var::d,Var::N);
		Var::ChargeCfg.push_back(Var::C);
		Var::Q =	0;	Var::Xq = Var::L.x/2;	Var::Yq = Var::L.y/2-Var::y_shift/2;	Var::Zq = 3.75e+3+Var::z_shift; Var::Rq1 = 3.00e+3;	Var::Rq3 = 1.50e+3;
		Var::C.disk(Var::Q, Var::Xq,Var::Yq,Var::Zq, Var::Rq1,Var::Rq3, Var::d,Var::N);
		Var::ChargeCfg.push_back(Var::C);
		Var::Q =    0;	Var::Xq = Var::L.x/2;	Var::Yq = Var::L.y/2+Var::y_shift/2;	Var::Zq = 6.75e+3+Var::z_shift; Var::Rq1 = 4.00e+3;	Var::Rq3 = 1.50e+3;
		Var::C.disk(Var::Q, Var::Xq,Var::Yq,Var::Zq, Var::Rq1,Var::Rq3, Var::d,Var::N);
		Var::ChargeCfg.push_back(Var::C);
		Var::Q =    0;	Var::Xq = Var::L.x/2;	Var::Yq = Var::L.y/2+Var::y_shift/2;	Var::Zq = 8.00e+3+Var::z_shift; Var::Rq1 = 4.00e+3;	Var::Rq3 = 0.50e+3;
		Var::C.disk(Var::Q, Var::Xq,Var::Yq,Var::Zq, Var::Rq1,Var::Rq3, Var::d,Var::N);
		Var::ChargeCfg.push_back(Var::C);
		// Endof Charge center parameters //

		Var::ThresholdOvershoot		= 10;											// % by which the E-field initiation threshold must be exceeded in the simulation domain in order to start a flash(for initiation threshold)
		Var::BCtype					= 1;											// Boundary conditions:
																					// = 0 ``Tin Can"
																					// = 1 Open BC
																					// = 2 Moving Capacitor Plates
																					// = 3 Free Space
		Var::InitiationType			= 1;											// Initiation type:
																					// = 1 randomly in available region
																					// = 2 at point of maximum E-field magnitude
																					// = 3 at implemented position
		Var::rho3dCalculationStep	= 10;											// Charge density is calculated and store every rho3dCalculationStep steps.
																					// = 0 3-D charge density is never calculated
		Var::AddNew					= true;											// Channel is allowed to propagate: Y/N
		Var::isBCerrorCalculated	= true;											// Error at the boundary is calculated at each step: Y/N
		Var::isBndXingPossible		= false;										// Channel is allowed to cross boundaries: Y/N
		Var::isChannelEquipotential	= true;											// Channel potential is adapted to ensure charge neutrality: Y/N
		Var::isEsEnergyCalculated   = true;											// Electrostatic energy is calculated at each step: Y/N
		Var::isFlashAccoutedInBC 	= true;											// Channel charge is accounted for in derivation of BC: Y/N
		Var::isInitiationPossible	= false;										// Initiation is possible in the simulation domain after charge load: Y/N (automatically turned to true if initition is possible)
		Var::isLinkXingPossible		= false;										// Channels crosses are allowed: Y/N
		Var::isRsDeveloped			= true;											// Return stroke development: Y/N
		Var::isInitiationPrevented	= false;										// Only simulate cloud electrical structure

		printf("..: Loading charge layers...\n");
		fprintf(file,	"..: Loading charge layers...\n");

		Var::ThresholdOvershoot		/= 100; 										// Convert % into decimal

		LoadTripoleModel(Var::I1,Var::I2,Var::Iscreen);

		Var::Eps_bf=0;
		for(int ii=0 ; ii<Var::N.x ; ii++) for(int jj=0 ; jj<Var::N.y ; jj++) for(int kk=0 ; kk<Var::N.z ; kk++)
			Var::Eps_bf += eps0*pow(foo::Eijk(ii,jj,kk,Var::phi,Var::d,Var::N)[0],2)/2*Var::d.x*Var::d.y*Var::d.z;

		printf("++: Finished loading charge layers!\n");
		fprintf(file,	"\ndone\n\n");
	}
	else if(Var::isChargeContentAssumed == true)
	{
		double Lr = 9e3;
		Var::L.z  = 21e3;
		Var::C.init("results/rho2d280000.dat","results/rhos2d280000.dat","results/d.dat","results/N.dat","results/z_gnd.dat", Var::d, Var::N, Var::z_gnd, Lr, Var::L.z);
		Var::ChargeCfg.push_back(Var::C);
		
//		for(int ii=0; ii<Var::N.x ; ii++) for(int jj=0; jj<Var::N.y ; jj++) for(int kk=0; kk<Var::N.z ; kk++)
//			Var::C.rho[ii][jj][kk] *= 1e9;

		Var::L.init(Var::N.x*Var::d.x,	Var::N.y*Var::d.y,	Var::L.z);				// _m	Size of the simulation domain
		
		// Initiate matrices dimensions //
		Var::E.init(Var::N.x,Var::N.y,Var::N.z);									// _V/m	Total electric field
		Var::phi.init(Var::N.x,Var::N.y,Var::N.z);									// _V	Total electric potential
		Var::phi_cha.init(Var::N.x,Var::N.y,Var::N.z);								// _V	Channel electric potential
		Var::phi_amb.init(Var::N.x,Var::N.y,Var::N.z);								// _V	Cloud electric potential
		Var::rho.init(Var::N.x,Var::N.y,Var::N.z);									// _C/m3	total charge density
		Var::Un.init(Var::N.x,Var::N.y,Var::N.z);									// Map of occupied grid points
		Var::phiNum.init(Var::N.z);													// _V	Total electric potential on a vertical axis in the center of simulation domain
		Var::EzNum.init(Var::N.z);													// _V/m	Total electric field on a vertical axis in the center of simulation domain
		// Endof Initiate matrices dimensions //
		
		//	We assume that the first link is somehow established, then only the propagation threshold needs to be exceeded to develop the flash.
		Var::Ec.init(2.10e+5,2.16e+5,-2.16e+5, Var::z_gnd,Var::d,Var::N);			// _V/m Initiation-, Positive channel propagation-, Negative channel propagation-threshold
		Var::Vd.init(0e+5,-0e+5, Var::z_gnd,Var::d,Var::N);							// _V/m Voltage drop in positive, negative channels
		Var::InitX	= Var::L.x/2;													// _m X-coordinate of initiation point
		Var::InitY	= Var::L.y/2;													// _m Y-coordinate of initiation point
		Var::InitZ	= 9.25e3;														// _m Z-coordinate of initiation point
		
		Var::InitiationPoint.init((int)round(Var::InitX/Var::d.x), (int)round(Var::InitY/Var::d.y),(int)round(Var::InitZ/Var::d.z));
																					// Initiation point with coordinates expressed as i,j,k and not x,y,z
/*		
		IO::write(Var::C.rho, "rho.dat");
		IO::write(Var::z_gnd, "z_gnd.dat");
		IO::write(Var::d.x, Var::d.y, Var::d.z, "dxyz.dat");
		IO::write(Var::N.x, Var::N.y, Var::N.z, "Nxyz.dat");
	
		CMatrix1D rrho3D(Var::N.x*Var::N.y*Var::N.z);
		int nn = 0;
		for(int kk=0 ; kk<Var::N.z ; kk++) for(int jj=0 ; jj<Var::N.y ; jj++) for(int ii=0 ; ii<Var::N.x ; ii++)
		{
			rrho3D[nn]		= Var::C.rho[ii][jj][kk]*1e9; //_nC
			nn++;
		};
		IO::write(rrho3D,"rho3d0.dat");
 */		
		// Endof Charge center parameters //
		
		Var::ThresholdOvershoot		= -100;											// % by which the E-field initiation threshold must be exceeded 
																					// arbitrarily fixed to 1 in this case.
		Var::BCtype					= 1;											// Boundary conditions:
																					// = 0 ``Tin Can"
																					// = 1 Open BC
																					// = 2 Moving Capacitor Plates
																					// = 3 Free Space
		Var::InitiationType			= 1;											// Initiation type:
																					// = 1 randomly in available region
																					// = 2 at point of maximum E-field magnitude
																					// = 3 at implemented position
		Var::rho3dCalculationStep	= 0;											// Charge density is calculated and store every rho3dCalculationStep steps.
																					// = 0 3-D charge density is never calculated
		Var::AddNew					= true;											// Channel is allowed to propagate: Y/N
		Var::isBCerrorCalculated	= false;										// Error at the boundary is calculated at each step: Y/N
		Var::isBndXingPossible		= false;										// Channel is allowed to cross boundaries: Y/N
		Var::isChannelEquipotential	= true;											// Channel potential is adapted to ensure charge neutrality: Y/N
		Var::isEsEnergyCalculated	= false;										// Electrostatic energy is calculated at each step: Y/N
		Var::isFlashAccoutedInBC 	= true;											// Channel charge is accounted for in derivation of BC: Y/N
		Var::isInitiationPossible	= true;											// Initiation is possible in the simulation domain after charge load: Y/N
		Var::isLinkXingPossible		= false;										// Channels crosses are allowed: Y/N
		Var::isRsDeveloped			= false;										// Return stroke development: Y/N
		Var::isInitiationPrevented	= false;										// Only simulate cloud electrical structure

		printf("..: Loading charge layers.\n");
		fprintf(file, "..: Loading charge layers.\n");


		Var::ThresholdOvershoot /= 100;												// Convert % into decimal

		printf("ii:\t Grid dimensions: [%d, %d, %d].\n", Var::N.x, Var::N.y, Var::N.z);
		printf("ii:\t Discretized lengths: [%lf m, %lf m, %lf m].\n", Var::d.x, Var::d.y, Var::d.z);
		printf("ii:\t Total simulation size: [%lf km, %lf km, %lf km].\n", Var::L.x/1e3, Var::L.y/1e3, Var::L.z/1e3);
		Var::C.reset(Var::d,Var::N);

		printf("ii:\t Charge layers are as follows:\n");
		for(it=Var::ChargeCfg.begin() ; it!=Var::ChargeCfg.end() ; it++)
		{
			cout<<*it;
			Var::C += *it;
		}

		printf("..: Initializing domain...\n");
		Var::SOR.init(Var::phi,Var::epsilon, Var::MaxStep, Var::d, Var::N, Var::C, Var::Un);
		BC::Apply(Var::BCtype,Var::phi,Var::C.rho,Var::d,Var::N);
		Var::SOR.Solve(Var::d,Var::N,Var::Un,Var::phi);
		printf("++: Domain initialized!\n");


		// Derive Electrostatic energy before the discharge //
		Var::Eps_bf=0;
		for(int ii=0 ; ii<Var::N.x ; ii++) for(int jj=0 ; jj<Var::N.y ; jj++) for(int kk=0 ; kk<Var::N.z ; kk++)
			Var::Eps_bf += eps0*pow(foo::Eijk(ii,jj,kk,Var::phi,Var::d,Var::N)[0],2)/2*Var::d.x*Var::d.y*Var::d.z;

		printf("++: Finished loading charge layers!\n");
		fprintf(file,	"\ndone\n\n");
	}

	printf("ii: Finding potential extrema.\n");
	fprintf(file, "ii: Finding potential extrema.\n");

	Var::rho		= foo::Globalrho(Var::phi,Var::d,Var::N);
	Var::Qtot_bf	= foo::TotalCharge(Var::rho,Var::d,Var::N);
	Var::phi.MinMax(Var::Vmin,Var::Vmax);
	Var::rho.MinMax(Var::rhoAmbMin,Var::rhoAmbMax);
	Var::phi_amb	=Var::phi;

	printf("ii:\t Extrema: [%lf, %lf].\n", Var::Vmin, Var::Vmax);
	printf("..: Initiating the tree.\n");
	fprintf(file,	"ii:\t Extrema: [%lf, %lf].\n", Var::Vmin, Var::Vmax);
	fprintf(file,	"..: Initiating the tree.\n");

	Var::isInitiationPossible = InitiateTree(Var::InitiationType, Var::InitiationPoint);
	if(Var::isInitiationPrevented == true)
	{
		Var::isInitiationPossible = false;
		printf("ii:\t Initiation not allowed.\n");
		fprintf(file,	"Initiation not allowed.\n");
	}
	printf("++: Finished initiating the tree!\n");
	fprintf(file,	"\ndone\n\n");

	printf("..: Growing the tree.\n");
	fprintf(file, "..: Growing the tree.\n");

	// Var::NumLinks = 0;
	if(Var::isInitiationPossible)
	{
		GrowTree(Var::AddNew);
	}

	// Calculate electrostatic energy after the discharge //
	Var::Eps_af=0;
	for(int ii=0 ; ii<Var::N.x ; ii++) for(int jj=0 ; jj<Var::N.y ; jj++) for(int kk=0 ; kk<Var::N.z ; kk++)
		Var::Eps_af += eps0*pow(foo::Eijk(ii,jj,kk,Var::phi,Var::d,Var::N)[0],2)/2*Var::d.x*Var::d.y*Var::d.z;
	Var::eFlux = foo::eFieldFlux(Var::phi,Var::d,Var::N);
	printf("ii: Finished growing the tree.\n");
	fprintf(file,	"\ndone\n\n");

	printf("..: Estimating the field and storing results.\n");
	fprintf(file, "..: Estimating the field and storing results.\n");

	StoreData(); /* Store the final data (still in the interim results folder). */
	//SortPathName(Var::NumLinks); /* Create the final results folder. */
	/* Move the final results to the final results folder.  This function call also deletes the interim
	 * results folder.
	 */
	//IO::moveFromOldDirectory();
	printf("++: Finished estimating the field and storing results!\n");
	fprintf(file,	"\ndone\n\n");

	printf("..: Summarizing results.\n");
	fprintf(file, "..: Summarizing results.\n");

	if(Var::BCtype == 0)
	{
		printf("ii:\t Closed Boundary Conditions\n");
		fprintf(file,	"ii:\t Closed Boundary Conditions\n");
	}
	else if(Var::BCtype == 1)
	{
		printf("ii:\t Open Boundary Conditions\n");
		fprintf(file,	"ii:\t Open Boundary Conditions\n");
	}
	else if(Var::BCtype == 2)
	{
		printf("ii:\t Moving Capacitor Plate Boundary Conditions\n");
		fprintf(file,	"ii:\t Moving Capacitor Plate Boundary Conditions\n");
	}
	else if(Var::BCtype == 3)
	{
		printf("ii:\t Free Space Boundary Conditions\n");
		fprintf(file,	"ii:\t Free Space Boundary Conditions\n");
	}
	if(Var::isFlashAccoutedInBC==true)
	{
		printf("ii:\t Updated Boundary Conditions\n");
		fprintf(file,	"ii:\t Updated Boundary Conditions\n");
	}
	else
	{
		printf("ii:\t Fast Boundary Conditions (no update)\n");
		fprintf(file,	"ii:\t Fast Boundary Conditions (no update)\n");
	}

	/************************************************
	*	Changed "%f" to "%lf"						*
	*************************************************/

	printf(			"ii:\t Threshold Overshoot = %3.2lf %%\n"							 ,Var::ThresholdOvershoot*100);
	fprintf(file,	"ii:\t Threshold Overshoot = %3.2lf %%\n"							 ,Var::ThresholdOvershoot*100);
	printf(			"ii:\t N = [%d, %d, %d]\n"										 ,Var::N.x,Var::N.y,Var::N.z);
	fprintf(file,	"ii:\t N = [%d, %d, %d]\n"										 ,Var::N.x,Var::N.y,Var::N.z);
	printf(			"ii:\t d = [%6.2lf, %6.2lf, %6.2lf]\n"								 ,Var::d.x,Var::d.y,Var::d.z);
	fprintf(file,	"ii:\t d = [%6.2lf, %6.2lf, %6.2lf]\n"								 ,Var::d.x,Var::d.y,Var::d.z);
	printf(			"ii:\t El.Stat. Energy before channel propagation    : %12.1lf J\n" ,Var::Eps_bf);
	fprintf(file,	"ii:\t El.Stat. Energy before channel propagation    : %12.1lf J\n" ,Var::Eps_bf);
	printf(			"ii:\t El.Stat. Energy after  channel propagation    : %12.1lf J\n" ,Var::Eps_af);
	fprintf(file,	"ii:\t El.Stat. Energy after  channel propagation    : %12.1lf J\n" ,Var::Eps_af);
	printf(			"ii:\t Total    Charge before channel propagation    : %12.6lf C\n" ,Var::Qtot_bf);
	fprintf(file,	"ii:\t Total    Charge before channel propagation    : %12.6lf C\n" ,Var::Qtot_bf);
	printf(			"ii:\t Total    Charge after  channel propagation    : %12.6lf C\n" ,Var::Qtot_af);
	fprintf(file,	"ii:\t Total    Charge after  channel propagation    : %12.6lf C\n" ,Var::Qtot_af);
	printf(			"ii:\t Positive Charge in the channel                : %12.6lf C\n" ,Var::QchannelPlus);
	fprintf(file,	"ii:\t Positive Charge in the channel                : %12.6lf C\n" ,Var::QchannelPlus);
	printf(			"ii:\t Negative Charge in the channel                : %12.6lf C\n" ,Var::QchannelMinus);
	fprintf(file,	"ii:\t Negative Charge in the channel                : %12.6lf C\n" ,Var::QchannelMinus);
	printf(			"ii:\t Flux of electric field through the boundaries : %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf C\n" ,Var::eFlux(0), Var::eFlux(1), Var::eFlux(2), Var::eFlux(3), Var::eFlux(4), Var::eFlux(5), Var::eFlux(6));
	fprintf(file,	"ii:\t Flux of electric field through the boundaries : %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf %4.2lf C\n" ,Var::eFlux(0), Var::eFlux(1), Var::eFlux(2), Var::eFlux(3), Var::eFlux(4), Var::eFlux(5), Var::eFlux(6));
	printf(			"ii:\t Equivalent charge                             : %4.2lf C\n" ,Var::Qtot_af);
	fprintf(file,	"ii:\t Equivalent charge                             : %4.2lf C\n" ,Var::Qtot_af);
	printf(			"ii:\t Charge layers:\n");
	fprintf(file,	"ii:\t Charge layers:\n");
	fclose(file);

	for(it=Var::ChargeCfg.begin() ; it!=Var::ChargeCfg.end() ; it++)
	{
		cout<<*it;
		// write(*it, fname);
	}	
	
	clock_t endTime = clock();
	clock_t runTime = endTime - startTime;
	printf(			"ii: Total run time: %lf s.\n\n", (double)runTime/CLOCKS_PER_SEC);
	fprintf(file,	"ii: Total run time: %lf s.\n", (double)runTime/CLOCKS_PER_SEC);
	fclose(file);

	/* Clear the charge structure so that it can be set during the next simulation run.
	 * Hence, it is possible to modify the input file (from the 'buck' program') so that
	 * subsequent simulation runs can work with different source data.
	 */
	Var::ChargeCfg.clear();
	return 1;
}
/**********************************************************************************/

/**********************************************************************************/
/* Load a Tripolar Structure for the ThunderCloud								  */
/**********************************************************************************/
void LoadTripoleModel(const double II1, const double II2, const double IIscreen)
{
	clock_t	sstartTime	= clock();
	clock_t	eendTime	= clock();
	clock_t	rrunTime	= clock();
	ListCharge::iterator  it;														// Table with all parameters of the charge configuation

	CriticalFields OOverShotEc((1+Var::ThresholdOvershoot)*Var::Ec.getParams()[0],
								Var::Ec.getParams()[1],Var::Ec.getParams()[2],Var::Ec.getParams()[3],
								Var::d,Var::N);
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

		rr_tmp	= EoverEk(TTl, OOverShotEc);
		rr_l	= (rr_tmp>=1)*1+(rr_tmp<1)*-1;
		rr_tmp	= EoverEk(TTr, OOverShotEc);
		rr_r	= (rr_tmp>=1)*1+(rr_tmp<1)*-1;

		if(TTl > TTr) {Swap::DBL(TTl,TTr); Swap::INT(rr_l,rr_l);}
		while(fabs(2*(TTr-TTl)/(TTr+TTl))>Var::epsilon)
	//	while(rr_r>1e-5)
		{
			kk++;
			TTav	= (TTl+TTr)/2;
			rr_tmp	= EoverEk(TTav, OOverShotEc);
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

		while(dt>Var::epsilon)
		{
			kk++;
			rr_tmp	=  EoverEk(TT+dt, OOverShotEc);
			rr		= (rr_tmp>=1)*1+(rr_tmp<1)*-1;
			while(rr>=0 && dt>Var::epsilon)
			{
				dt /=2;
				rr_tmp	=  EoverEk(TT+dt, OOverShotEc);
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
		rr = EoverEk(TT, OOverShotEc);
		TT *= 1/rr;
	}



	/**********************************************************************************/
	/* Estimate how much the initiation threshold is exceeded						  */
	/**********************************************************************************/
	rr = EoverEk(TT, OOverShotEc);

	/**********************************************************************************/
	/* display results																  */
	/**********************************************************************************/
	cout<<"T = "<<setw(12)<<TT<<" ; max(E/((1+alpha)*Ek)) = "<<setw(13)<<rr<<endl;
	cout<<"   Nb of iterations    = "<<kk+1<<endl;;
	cout<<"ii:\t Initiation Threshold is exceeded by "<<(rr*(1+Var::ThresholdOvershoot)-1)*100<<"%\n";
	cout<<"ii:\t Final values of the charge domains:\n";
	cout<<"ii:\t\t [Q1 Q2 Q3 Q4]       = [ ";
	for (it=Var::ChargeCfg.begin() ; it!=Var::ChargeCfg.end() ; it++)
		cout<<it->getParams()[0]<<" ";
	cout<<"] C\n";
	eendTime = clock();
	rrunTime = eendTime - sstartTime;
	putchar('\n');
	printf("ii:\t Run time for layers load: %lf s.\n", (double)rrunTime/CLOCKS_PER_SEC);
};

double EoverEk(const double tt, CriticalFields OOverShotEc)
{
	double					MMaxRatio(0);
	double					rratio(0);
	int						nn;														// Counter of iteration
	double					QQ, XXq,YYq,ZZq, RRq1,RRq2,RRq3;
	Charge					CC(Var::d,Var::N);										// Temporary charge structure for SOR derivation
	ListCharge::iterator	it;														// Table with all parameters of the charge configuation

	nn=0; QQ=0;
	for (it=Var::ChargeCfg.begin() ; it!=Var::ChargeCfg.end() ; it++)
	{
		/* retrieve tripole geometrical parameters */
		XXq = it->getParams()[1];
		YYq = it->getParams()[2];
		ZZq = it->getParams()[3];
		RRq1= it->getParams()[4];
		RRq2= it->getParams()[5];
		RRq3= it->getParams()[6];

		/* evaluate layers charge */
		if(nn==0) QQ = -Var::I2*tt;
		if(nn==1) QQ = (-Var::I1+Var::I2)*tt;
		if(nn==2) QQ = (Var::I1*tt);
		if(nn==3) QQ = -Var::Iscreen*tt - ( -Var::I2*tt + (-Var::I1+Var::I2)*tt + (Var::I1*tt));
		if(nn>=4) QQ = 0;

		/* update layers */
		it->disk(QQ, XXq,YYq,ZZq, RRq1,RRq3, Var::d,Var::N);
		CC+=*it;
		nn++;
	}

	/* solve fields */
	Var::SOR.init(Var::phi,Var::epsilon, Var::MaxStep, Var::d, Var::N, CC, Var::Un);
	 BC::Apply(Var::BCtype,Var::phi,CC.rho,Var::d,Var::N);
	Var::SOR.Solve(Var::d,Var::N,Var::Un,Var::phi);

	/* Find ddeltaE */
	for(int kk=0; kk<Var::N.z; kk++) for(int jj=1; jj<Var::N.y-1; jj++) for(int ii=0; ii<Var::N.x; ii++)
	{
		rratio = foo::Eijk(ii,jj,kk,Var::phi,Var::d,Var::N)[0]/OOverShotEc.initiation[kk];
		if(rratio>=MMaxRatio)	MMaxRatio	= rratio;
	};
	return MMaxRatio;
}
/**************************************************************************************/

/**************************************************************************************/
/* Initiate Tree (Find initiation point)											  */
/**************************************************************************************/
bool InitiateTree(const int IInitiationType, Point& IInitiationPoint)
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
		
		for(int ii=1; ii<Var::N.x-1; ii++) for(int jj=1; jj<Var::N.y-1; jj++) for(int kk=1; kk<Var::N.z-1; kk++)
		{
			/* Choose the points which are likely to launch a flash	*/
			if(foo::Eijk(ii,jj,kk, Var::phi,Var::d,Var::N)[0] >= (1+0*Var::ThresholdOvershoot)*Var::Ec.initiation[kk])
			{
				CCandidate.i=ii;
				CCandidate.j=jj;
				CCandidate.k=kk;
				LListOfCandidates.push_back(CCandidate);
				ccptPoints++;
			}
		};

		printf("ii:\t Number of possible candidate endpoints: %d.\n", ccptPoints);
		//printf("ii:\t Maximum threshold overreach: %lf%%.\n", maxPerExc);

		if(ccptPoints>0)
		{
			/* Pick a point randomly in this list */
			RRandomPoint = rand()%ccptPoints+1;
			for(iit = LListOfCandidates.begin() ; iit!=LListOfCandidates.end() ; iit++)
			{
				if(kk==RRandomPoint)
				{
					IInitiationPoint = *iit;
					Var::Un[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k] = 1;
					Var::phi0 = Var::phi[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k];
					break;
				}
				kk++;
			}
			isInitiated = true;
		}
		else
		{
			isInitiated = false;
			cout<<"ii:\t Cannot initiate flash (insufficient field).\n";
		}
	}
	else if(IInitiationType == 2)
	{
		/******************/
		/* Init at Emax   */
		/******************/
		double Emax = 0;
		for(int ii=1; ii<Var::N.x-1; ii++) for(int jj=1; jj<Var::N.y-1; jj++) for(int kk=1; kk<Var::N.z-1; kk++)
			if(foo::Eijk(ii,jj,kk, Var::phi,Var::d,Var::N)[0] >= Var::Ec.initiation[kk])
				if(foo::Eijk(ii,jj,kk, Var::phi,Var::d,Var::N)[0] >= Emax)
				{
					Emax=foo::Eijk(ii,jj,kk, Var::phi,Var::d,Var::N)[0];
					IInitiationPoint.init(ii,jj,kk);
				};
		Var::Un[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k] = 1;
		Var::phi0 = Var::phi[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k];
		if(foo::Eijk(IInitiationPoint.i,IInitiationPoint.j,IInitiationPoint.k, Var::phi,Var::d,Var::N)[0] < Var::Ec.initiation[IInitiationPoint.k])
		{
			isInitiated = false;
			cout<<"ii:\t Cannot initiate flash (insufficient field).\n";
		}
		else {isInitiated = true;}
	}
	else if(IInitiationType == 3)
	{
		/******************/
		/* Fix Init Point */
		/******************/
		if(foo::Eijk(IInitiationPoint.i, IInitiationPoint.j, IInitiationPoint.k, Var::phi,Var::d,Var::N)[0] >= Var::Ec.initiation[IInitiationPoint.k])
		{
			Var::Un[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k] = 1;
			Var::phi0 = Var::phi[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k];
			isInitiated = true;
		}
		else
		{
			cout<<"ii:\t Cannot initiate flash (insufficient field).\n";
			isInitiated = false;
		}
	}
	else
	{
		cout<<"ee:\t Wrong Choice for initiation procedure."<<endl;
		exit(3);
	}
	/* Display coordinates of initiation point */
	printf("++:\t Discharge initated at [%lf km, %lf km, %lf km].\n",
			IInitiationPoint.i*Var::d.x/1e3, IInitiationPoint.j*Var::d.y/1e3, (IInitiationPoint.k*Var::d.z+Var::z_gnd)/1e3);
//	SOR.Solve(d,N,Un,phi);

	/* Derive field and potential before discharge along the central vertical axis */
	for(int kk=0 ; kk<Var::N.z ; kk++)
	{
		Var::EzNum[kk]	= foo::Eijk((Var::N.x-1)/2,(Var::N.y-1)/2,kk,Var::phi,Var::d,Var::N)[3];
		Var::phiNum[kk]	= Var::phi((Var::N.x-1)/2,(Var::N.y-1)/2,kk);
	}

	/* Store Initiation data */
	IO::write(Var::phiNum,	"phiNumBF.dat");
	IO::write(Var::EzNum,	"EnumBF.dat");
	return isInitiated;
}
/**********************************************************************************/

/**********************************************************************************/
/* Grow the Tree																  */
/**********************************************************************************/
void GrowTree(bool AAddNew)
{
	bool					isConnectedToTheGround;									// Check connection to the ground
	char					nname3D[50];											// File name for storage of charge density
	int						ccptLinks(0);											// Current link iteration
	double					BBndErrorTmp = 0;										// Error at the boundaries at the current step
	double					EEsEnergyTmp = 0;										// Electrostatic energy at the current step
	CMatrix1D				TTotalPotentialTemp(Var::N.z);							// Total potential on the central vertical axis at the current step
	CMatrix1D				TTotalEfieldTemp(Var::N.z);								// Total eField on the central vertical axis at the current step
	CMatrix2D				pphi2D_cha(Var::N.y,Var::N.z);							// Potential induced by the channel in the y-z plane
	CMatrix3D				pphi_chaTmp;											// Potential due to the channel in the domain before update of the boundary conditions
	CMatrix3D				rrho_cha(Var::N.x,Var::N.y,Var::N.z);					// 3-D Matrix for channel induced charge density
	CMatrix1D				rrho3D(Var::N.x*Var::N.y*Var::N.z);						// Total charge density everywhere at the current step
	int						nn = 0;													// Blind variable
	ListLink::iterator		LLit;													// index on the list of established links
	ListVector::iterator	it;														// index on a list of vectors
	Charge					CC;														// Charge transfer at the current step
	Potential				PP;														// Channel potential at the current step
	Vector					pp;														// x-, y-, z-components and norm of the dipole moment at the current step
	Vector					BBndError;												// Potential at the boundary before update, after update and relative difference between the two values

	nn = 0;
	for(int kk=0 ; kk<Var::N.z ; kk++)
	{
		TTotalEfieldTemp[kk]	= foo::Eijk((Var::N.x-1)/2,(Var::N.y-1)/2,kk,Var::phi,Var::d,Var::N)[3];
		TTotalPotentialTemp[kk]	= Var::phi((Var::N.x-1)/2,(Var::N.y-1)/2,kk);
		if(Var::isEsEnergyCalculated == true && Var::rho3dCalculationStep == 0)
			for(int jj=0 ; jj<Var::N.y ; jj++)
				for(int ii=0 ; ii<Var::N.x ; ii++)
					EEsEnergyTmp += eps0*pow(foo::Eijk(ii,jj,kk,Var::phi,Var::d,Var::N)[0],2)/2*Var::d.x*Var::d.y*Var::d.z;
		if(Var::isEsEnergyCalculated == false && Var::rho3dCalculationStep != 0)
			for(int jj=0 ; jj<Var::N.y ; jj++) for(int ii=0 ; ii<Var::N.x ; ii++)
			{
				rrho3D[nn]		= foo::rhoijk(ii,jj,kk,Var::phi,Var::d,Var::N)*1e+9; //_nC
				nn++;
			};
		if(Var::isEsEnergyCalculated == true && Var::rho3dCalculationStep != 0)
			for(int jj=0 ; jj<Var::N.y ; jj++) for(int ii=0 ; ii<Var::N.x ; ii++)
			{
				EEsEnergyTmp += eps0*pow(foo::Eijk(ii,jj,kk,Var::phi,Var::d,Var::N)[0],2)/2*Var::d.x*Var::d.y*Var::d.z;
				rrho3D[nn] = foo::rhoijk(ii,jj,kk,Var::phi,Var::d,Var::N)*1e+9; //_nC
				nn++;
			};
	}
	if(Var::rho3dCalculationStep != 0)
	{
		IO::write(Var::rho3dCalculationStep,"rho3dStep.dat");
		IO::write(rrho3D,					"rho3d0.dat");
	}

	pp = foo::DipoleMoment(Var::QchannelPlus,Var::phi_cha,Var::Un,Var::L,Var::d,Var::N);
	Var::ChannelPotential.push_back(Var::phi0);
	Var::DischargeDipoleMoment.push_back(pp);
	Var::CarriedCharge.push_back(Var::QchannelPlus);
	Var::EsEnergy.push_back(EEsEnergyTmp);
	Var::TotalEfield.push_back(TTotalEfieldTemp);
	Var::TotalPotential.push_back(TTotalPotentialTemp);

	/* Clear lists for multiple simulation runs */
	Var::EstablishedLinks.clear();
	Var::NumberOfCandidates.clear();
	Var::MaximumCandidateOverreach.clear();

	while(AAddNew==true && ccptLinks>=0)
	{
		ccptLinks++;
		AAddNew	= Tree::AddNewLink(Var::d,Var::N,Var::Un,Var::phi, Var::Ec,Var::Vd, Var::InitiationPoint,Var::EstablishedLinks,
							 Var::isBndXingPossible,  Var::isRsDeveloped,
							 Var::isLinkXingPossible, Var::isChannelEquipotential);

		/* Check Connection to the ground */
		isConnectedToTheGround = false;
		for(LLit = Var::EstablishedLinks.begin(); LLit != Var::EstablishedLinks.end() ; LLit++)
			if(LLit->end.k == 0)
			{
				isConnectedToTheGround = true;
				break;
			}
		if(Var::isChannelEquipotential==false)
		{
			 BC::Update(Var::isFlashAccoutedInBC,Var::BCtype,Var::phi,Var::rhoAmbMin,Var::rhoAmbMax,Var::d,Var::N);
			Var::SOR.Solve(Var::d,Var::N,Var::Un,Var::phi);
		};
		if(Var::isChannelEquipotential==true)
		{
			if(isConnectedToTheGround == false)
			{
				Var::phi0 = Tree::fMinSearch(Var::phi0, Var::QchannelPlus, Var::Vmin,Var::Vmax , Var::epsilon,Var::MaxStep, Var::phi_cha,Var::phi_amb ,Var::Un, Var::InitiationPoint,Var::EstablishedLinks, Var::d,Var::N);
			}
			if (isConnectedToTheGround == true)
			{
				Var::phi0 = Tree::EqualizeAtGroundPotential(Var::epsilon,Var::MaxStep, Var::phi_cha,Var::phi_amb ,Var::Un, Var::InitiationPoint,Var::EstablishedLinks, Var::d,Var::N);
			}
			Var::ChannelPotential.push_back(Var::phi0);

			if(Var::isBCerrorCalculated == false)
			{
				 BC::Update(Var::isFlashAccoutedInBC,Var::BCtype,Var::phi_cha,Var::rhoAmbMin,Var::rhoAmbMax,Var::d,Var::N);
				Var::phi	= Var::phi_amb+Var::phi_cha;
			}

			if(Var::isBCerrorCalculated ==true)
			{
				pphi_chaTmp	= Var::phi_cha;
				 BC::Update(Var::isFlashAccoutedInBC,Var::BCtype,Var::phi_cha,Var::rhoAmbMin,Var::rhoAmbMax,Var::d,Var::N);
				Var::phi	= Var::phi_amb+Var::phi_cha;

				BBndError.x = 0;
				BBndError.y = 0;
				BBndError.z = 0;
				for(int kk=0; kk<Var::N.z; kk++) for(int jj=0; jj<Var::N.y; jj++) for(int ii=0; ii<Var::N.x; ii++)
					if( ii == 0 || ii == Var::N.x-1 ||jj == 0 || jj == Var::N.y-1 || kk == 0 || kk == Var::N.z-1)
					{
						BBndErrorTmp = fabs(Var::phi_cha(ii,jj,kk) - pphi_chaTmp(ii,jj,kk));
						if(BBndErrorTmp >= fabs(BBndError.x-BBndError.y))
						{
							BBndError.x = pphi_chaTmp(ii,jj,kk);
							BBndError.y = Var::phi_cha(ii,jj,kk);
							BBndError.z = Var::phi_amb(ii,jj,kk);;
						}
					}
				Var::BndUpdateErrors.push_back(BBndError);
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
			PP.init(Var::phi_cha,Var::Un);
			Var::SOR.init(Var::phi_cha, Var::epsilon,Var::MaxStep, Var::d, Var::N, PP, Var::Un);
			Var::SOR.Solve(Var::d,Var::N,Var::Un,Var::phi_cha);
			Var::phi = Var::phi_amb+Var::phi_cha;
			*/

			EEsEnergyTmp = 0;
			if(Var::rho3dCalculationStep != 0 && ccptLinks%Var::rho3dCalculationStep==0)
			{
				sprintf(nname3D,"rho3d%d.dat", ccptLinks);
				nn          = 0;
			}
			for(int kk=0 ; kk<Var::N.z ; kk++)
			{
				TTotalEfieldTemp[kk]	= foo::Eijk((Var::N.x-1)/2,(Var::N.y-1)/2,kk,Var::phi,Var::d,Var::N)[3];
				TTotalPotentialTemp[kk]	= Var::phi((Var::N.x-1)/2,(Var::N.y-1)/2,kk);
				if(Var::isEsEnergyCalculated == true && Var::rho3dCalculationStep == 0)
					for(int jj=0 ; jj<Var::N.y ; jj++) for(int ii=0 ; ii<Var::N.x ; ii++)
							EEsEnergyTmp += eps0*pow(foo::Eijk(ii,jj,kk,Var::phi,Var::d,Var::N)[0],2)/2*Var::d.x*Var::d.y*Var::d.z;
				if(Var::isEsEnergyCalculated == false && Var::rho3dCalculationStep != 0 && ccptLinks%Var::rho3dCalculationStep==0)
					for(int jj=0 ; jj<Var::N.y ; jj++) for(int ii=0 ; ii<Var::N.x ; ii++)
					{
						rrho3D[nn] = foo::rhoijk(ii,jj,kk,Var::phi,Var::d,Var::N)*1e+9; //_nC
						nn++;
					};
				if(Var::isEsEnergyCalculated == true && Var::rho3dCalculationStep != 0)
					for(int jj=0 ; jj<Var::N.y ; jj++) for(int ii=0 ; ii<Var::N.x ; ii++)
					{
						EEsEnergyTmp += eps0*pow(foo::Eijk(ii,jj,kk,Var::phi,Var::d,Var::N)[0],2)/2*Var::d.x*Var::d.y*Var::d.z;
						if(ccptLinks%Var::rho3dCalculationStep==0)
						{
							rrho3D[nn] = foo::rhoijk(ii,jj,kk,Var::phi,Var::d,Var::N)*1e+9; //_nC
							nn++;
						}
					};
			}
			if(Var::rho3dCalculationStep != 0 && ccptLinks%Var::rho3dCalculationStep==0)
				IO::write(rrho3D, nname3D);

			pp = foo::DipoleMoment(Var::QchannelPlus,Var::phi_cha,Var::Un,Var::L,Var::d,Var::N);
			Var::DischargeDipoleMoment.push_back(pp);
			Var::CarriedCharge.push_back(Var::QchannelPlus);
			Var::EsEnergy.push_back(EEsEnergyTmp);
			Var::TotalEfield.push_back(TTotalEfieldTemp);
			Var::TotalPotential.push_back(TTotalPotentialTemp);
		};
		cout<<"ii:\t Nb of established links: "<<ccptLinks<<endl;

		/* Write interim results after every tenth link is added. */
		//if((ccptLinks % 10) == 0)
		//	StoreData();
	}
	cout<<"Total Number of links: "<<ccptLinks<<endl;
	rrho_cha			= foo::Globalrho(Var::phi_cha,Var::d,Var::N);
	Var::QchannelPlus	= foo::ChannelChargePositive(rrho_cha,Var::Un,Var::d,Var::N);
	Var::QchannelMinus	= foo::ChannelChargeNegative(rrho_cha,Var::Un,Var::d,Var::N);
	Var::Qtot_af		= foo::TotalCharge(foo::Globalrho(Var::phi,Var::d,Var::N),Var::d,Var::N);

	bool FFlagInit=false;
	for(int ii=1; ii<Var::N.x-1; ii++) for(int jj=1; jj<Var::N.y-1; jj++) for(int kk=1; kk<Var::N.z-1; kk++)
		if(foo::Eijk(ii,jj,kk, Var::phi,Var::d,Var::N)[0] >= (1+Var::ThresholdOvershoot)*Var::Ec.initiation[kk])
		{
//			cout<<"Initiation Threshold exceeded at: ["<<ii*d.x<<" "<<jj*d.y<<" "<<kk*d.z<<"]\n";
			FFlagInit=true;
			break;
		};
	if(FFlagInit==false) cout<<"*** No further initiation possible. ***\n";
	else{cout<<"*** Further initiation possible. ***\n";}
	if (Var::isBCerrorCalculated == true)
		cout<<"Error on the boundaries\n"<<setw(12)<<"phiCha(N)"<<setw(12)<<"phiCha(N+1)"<<setw(16)<<"phiAmb"<<endl;
};
/**************************************************************************************/

/**************************************************************************************/
/* Store main data																	  */
/**************************************************************************************/

/* SAM:	Remember, the IO class takes care of putting the results in the correct folder.
 *		Therefore, specifying file names such as "results/ChannelPotential.dat" is unnecessary.
 *
 * btw, SAM = "Slight Additions and Modifications"
 */
void StoreData()
{
	//clock_t		startTime = clock();
	//clock_t		endTime;
	//clock_t		runTime;

	CMatrix3D				rrho_amb(Var::N.x, Var::N.y, Var::N.z);
	CMatrix2D				rrho_amb_yz(Var::N.y, Var::N.z);
	CMatrix2D				rrho_amb_xz(Var::N.x, Var::N.z);
	CMatrix2D				EE2D(	Var::N.y, Var::N.z);
	CMatrix2D				EEx2D(	Var::N.y, Var::N.z);
	CMatrix2D				EEy2D(	Var::N.y, Var::N.z);
	CMatrix2D				EEz2D(	Var::N.y, Var::N.z);
	CMatrix2D				pphi2D(	Var::N.y, Var::N.z);
	ListCMatrix1D			CChargeLayers;
	ListCharge::iterator	iit;

	printf("..: Writing results.\n");

	// Derive required values before storage //
	Var::E		= foo::GlobalE(	Var::phi,		Var::d,Var::N);
	Var::rho	= foo::Globalrho(Var::phi,		Var::d,Var::N);
	rrho_amb	= foo::Globalrho(Var::phi_amb,	Var::d,Var::N);

	for(int kk=0 ; kk<Var::N.z ; kk++)
	{
		Var::EzNum[kk]	= foo::Eijk((Var::N.x-1)/2,(Var::N.y-1)/2,kk,Var::phi,Var::d,Var::N)[3];
		Var::phiNum[kk]	= Var::phi((Var::N.x-1)/2,(Var::N.y-1)/2,kk);
		for(int ii=0 ; ii<Var::N.x ; ii++)
			rrho_amb_xz[ii][kk]= foo::rhoijk(ii,(Var::N.y-1)/2,kk,Var::phi_amb,Var::d,Var::N)*1e9;

		for(int jj=0 ; jj<Var::N.y ; jj++)
		{
			rrho_amb_yz[jj][kk]= foo::rhoijk( (Var::N.x-1)/2,jj,kk, Var::phi_amb,	Var::d,Var::N)*1e9;
			EE2D[jj][kk]	   =   foo::Eijk( (Var::N.x-1)/2,jj,kk, Var::phi,		Var::d,Var::N)[0];
			EEx2D[jj][kk]	   =   foo::Eijk( (Var::N.x-1)/2,jj,kk, Var::phi,		Var::d,Var::N)[1];
			EEy2D[jj][kk]	   =   foo::Eijk( (Var::N.x-1)/2,jj,kk, Var::phi,		Var::d,Var::N)[2];
			EEz2D[jj][kk]	   =   foo::Eijk( (Var::N.x-1)/2,jj,kk, Var::phi,		Var::d,Var::N)[3];
			pphi2D[jj][kk]	   =    Var::phi( (Var::N.x-1)/2,jj,kk);
			for(int ii=0 ; ii<Var::N.x ; ii++)
			{
				Var::rho(ii,jj,kk)	*= 1e9;
				rrho_amb(ii,jj,kk)	*= 1e9;
			}
		}
	}
	for (iit=Var::ChargeCfg.begin() ; iit!=Var::ChargeCfg.end() ; iit++)
		CChargeLayers.push_back(iit->getParams());

	// Store altitude of the ground plane //
	IO::write(Var::Ec.getParams()[3],	"z_gnd.dat");

	// Store Values in main plane //
	IO::write(Var::rho,		"rho.dat");
	IO::write(rrho_amb,		"rhoAmb.dat");
	IO::write(rrho_amb_yz,	"rhoAmbYZ.dat");
	IO::write(rrho_amb_xz,	"rhoAmbXZ.dat");
	IO::write(EE2D,			"E2D.dat");
	IO::write(EEx2D,		"Ex2D.dat");
	IO::write(EEy2D,		"Ey2D.dat");
	IO::write(EEz2D,		"Ez2D.dat");
	IO::write(pphi2D,		"phi2D.dat");

	/**************************************************************************
	double	QchannelPlus;
	double	QchannelMinus;
	double	tmp;
	ListDouble LL1(IO::read("rho.dat"));
	ListDouble LL2(IO::read("rhoAmb.dat"));
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
			QchannelPlus  += tmp*(float)Var::d.x*(float)Var::d.y*(float)Var::d.z;
		if(tmp<=0)
			QchannelMinus += tmp*(float)Var::d.x*(float)Var::d.y*(float)Var::d.z;
		itLD2++;
	}
	cout.precision(10);
	cout<<"Positive Charge in the channel (rounded)      : "<<QchannelPlus<<" C\n";
	cout<<"Negative Charge in the channel (rounded)      : "<<QchannelMinus<<" C\n"<<endl;

	QchannelPlus  = 0;
	QchannelMinus = 0;
	tmp			  = 0;
	for(int ii=1 ; ii<Var::N.x-1 ; ii++) for(int jj=1 ; jj<Var::N.y-1 ; jj++) for(int kk=1 ; kk<Var::N.z-1 ; kk++)
	{
		tmp = Var::rho(ii,jj,kk)-rrho_amb(ii,jj,kk);
		if (Var::Un(ii,jj,kk)!=0)
		{
			if(tmp>=0)
				QchannelPlus  += tmp*Var::d.x*Var::d.y*Var::d.z*1e-9;
			if(tmp<=0)
				QchannelMinus += tmp*Var::d.x*Var::d.y*Var::d.z*1e-9;
		}
	};
	cout<<"Positive Charge in the channel                : "<<QchannelPlus<<" C\n";
	cout<<"Negative Charge in the channel                : "<<QchannelMinus<<" C\n"<<endl;
	**************************************************************************/

	// Store channel potential, charge and dipole moment //
	IO::write(Var::ChannelPotential,		"ChannelPotentials.dat");
	IO::write(Var::CarriedCharge,			"CarriedCharge.dat");
	IO::write(Var::DischargeDipoleMoment,	"DischargeDipoleMoment.dat");

	// Store values of the electrostatic energy in the domain (if calculated) //
	if(Var::isEsEnergyCalculated == true)
		IO::write(Var::EsEnergy,			"EsEnergy.dat");

	// Store values of the boundary errors in the domain (if calculated) //
	if(Var::isBCerrorCalculated == true)
		IO::write(Var::BndUpdateErrors,		"BndError.dat");

	// Store d, N and Initiation point //
	IO::write(Var::d.x, Var::d.y, Var::d.z,	"dxyz.dat");
	IO::write(Var::N.x, Var::N.y, Var::N.z,	"Nxyz.dat");
	IO::write(Var::InitiationPoint.i*Var::d.x, Var::InitiationPoint.j*Var::d.y, Var::InitiationPoint.k*Var::d.z, Var::InitR, "InitPoint.dat");

	// Store charge center parameters //
	IO::write(CChargeLayers,	"ChargeLayers.dat");

	// Store Values of the fields and potential on the principal axis //
	IO::write(Var::phiNum,			"phiNumAF.dat");
	IO::write(Var::EzNum,			"EnumAF.dat");
	IO::write(Var::Ec.initiation,	"Einitiation.dat");
	IO::write(Var::Ec.negative,		"EthNegative.dat");
	IO::write(Var::Ec.positive,		"EthPositive.dat");
	IO::write(Var::TotalPotential,	"TotalPotential.dat");
	IO::write(Var::TotalEfield,		"TotalEfield.dat");

	// Store Values of the field and potential everywhere //
	IO::write(Var::phi,				"phi.dat");
	IO::write(Var::E,				"E.dat");

	// Store the list of Established Links //
	IO::write(Var::EstablishedLinks,	"EstablishedLinks.dat");
}
/**************************************************************************************/
