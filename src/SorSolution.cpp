/* SorFunctions.cpp */

#include "SorSolution.h"
/**************************************************************************************/

SorSolution::SorSolution(CMatrix3D& pphi, double eepsilon, int MMaxStep, StepsSizes dd, BoxSteps NN, Charge& CC, const CMatrix3D& UUn)
{SorSolution::init(pphi, eepsilon,MMaxStep, dd,NN, CC, UUn);}
																				// Derivation of SORcoeff and Norm of the error.

void SorSolution::init(CMatrix3D& pphi, double eepsilon, int MMaxStep, StepsSizes dd, BoxSteps NN, Charge& CC, const CMatrix3D& UUn)
{
	double gghat = -2*(1/pow(dd.x,2)+1/pow(dd.y,2)+1/pow(dd.z,2));
	type	= ChargeDistribution;
	epsilon	= eepsilon;
	MaxStep	= MMaxStep;
	
	a		= 1 / pow(dd.x,2) / gghat;
    b		= 1 / pow(dd.x,2) / gghat;
    c		= 1 / pow(dd.y,2) / gghat;
    d		= 1 / pow(dd.y,2) / gghat;
    e		= 1 / pow(dd.z,2) / gghat;
    f		= 1 / pow(dd.z,2) / gghat;
    g		= 1 ;
	eErrDen = 0 ;

	h.init(NN.x,NN.y,NN.z);
	for(int  kk=0 ; kk<NN.z ; kk++) for(int  jj=0 ; jj<NN.y ; jj++) for(int  ii=0 ; ii<NN.x ; ii++)
	{
		if ( ii!=0 && ii!=NN.x-1 &&  jj!=0 &&  jj!=NN.y-1 &&  kk!=0 && kk!=NN.z-1)
		{h[ii][jj][kk] = (-1/PMC.eps0)*CC.rho[ii][jj][kk]/gghat;}
		else
		{h[ii][jj][kk]	= 0;}
		eErrDen += (UUn[ii][jj][kk]==0)*pow(h[ii][jj][kk],2);
	}
}// Init Charge

SorSolution::SorSolution(CMatrix3D& pphi, double eepsilon, int MMaxStep, StepsSizes dd, BoxSteps NN, Potential& PP, CMatrix3D& UUn)
{SorSolution::init(pphi, eepsilon,MMaxStep, dd,NN, PP,UUn);}

void SorSolution::init(CMatrix3D& pphi, double eepsilon, int MMaxStep, StepsSizes dd, BoxSteps NN, Potential& PP, CMatrix3D& UUn)
{
	double gghat	= -2*(1/pow(dd.x,2)+1/pow(dd.y,2)+1/pow(dd.z,2));
	type			= PotentialDistribution;
	epsilon			= eepsilon;
	MaxStep			= MMaxStep;
	
	//	UUn				= PP.Un;
	a				= 1 / pow(dd.x,2) / gghat;
    b				= 1 / pow(dd.x,2) / gghat;
    c				= 1 / pow(dd.y,2) / gghat;
    d				= 1 / pow(dd.y,2) / gghat;
    e				= 1 / pow(dd.z,2) / gghat;
    f				= 1 / pow(dd.z,2) / gghat;
    g				= 1 ;
	eErrDen			= 0 ;
	UUn				= PP.Un;
	
	if(PP.getEquiPotential()==true)
	{
		double VVo		= PP.getVo();	
//		h.init(1,1,1);
//		h(0,0,0) = VVo;									// Here h is only used to store the value
														// of the potential
		for(int ii = 0 ; ii<NN.x ; ii++) for(int jj = 0 ; jj<NN.y ; jj++) for(int kk = 0 ; kk<NN.z ; kk++)
		{
			if(PP.Un[ii][jj][kk] !=0) 
			{
				if ( ii!=0 && ii!=NN.x-1 &&  jj!=0 &&  jj!=NN.y-1 &&  kk!=0 && kk!=NN.z-1)
				{pphi[ii][jj][kk]	= VVo;}
				eErrDen	+= pow(VVo,2);
			}
		};
	}
	else if(PP.getEquiPotential()==false)
	{
		for(int ii = 0 ; ii<NN.x ; ii++) for(int jj = 0 ; jj<NN.y ; jj++) for(int kk = 0 ; kk<NN.z ; kk++)
		{
			if(PP.Un[ii][jj][kk]!=0)
			{
				if ( ii!=0 && ii!=NN.x-1 &&  jj!=0 &&  jj!=NN.y-1 &&  kk!=0 && kk!=NN.z-1)
				{pphi[ii][jj][kk]	= PP.rho[ii][jj][kk];}
				eErrDen	+= pow(PP.rho[ii][jj][kk],2);
//				if (fabs(PP.rho[ii][jj][kk])>=eErrDen) eErrDen = fabs(PP.rho[ii][jj][kk]);
			}
		};
	}	
} // init Potential 

void SorSolution::Solve(StepsSizes dd, BoxSteps NN, const CMatrix3D& UUn, CMatrix3D& pphi)
{
	clock_t sstartTime	= clock();
	clock_t eendTime	= sstartTime;
	clock_t rrunTime	= eendTime-sstartTime;
	
//	double*** UUn	= UUUn.pMatrix3d;
//	double*** pphi	= ppphi.pMatrix3d;
	
	/**********************************************************************************/
	/* Solution 2: Fast asymetric solution.											  */
	/**********************************************************************************/
	double		rres	= 0;
	double		rrb		= cos( PMC.pi / NN.max());
	double		wwb		= 2 / ( 1 + sqrt( 1 - pow(rrb,2) ) );
	double		eErrNum = 0;						// Numerator of the Remaining Error
	double		eErr	= 0;						// Remaining Error
	int			sstep	= 0;						// Step Counter
		
	/**********************************************************************************/
	/* NN.x represents the total number of points NN.x-1 is the maximum reachable     */ 
	/* index for i since L.x = (NN.x - 1)*dd.x                                        */
	/* "cout<<pphi[NN.x][NN.y][NN.z]<<endl<<endl<<endl;" returns a SigBus error.	  */
	/* Same comments for y and z direction											  */
	/**********************************************************************************/
	eErrNum	= eErrDen;
	eErr	= eErrNum/eErrDen;
	sstartTime = clock();
	while(eErr>epsilon && sstep<MaxStep)
	{
		eErrNum	= 0;
		for(int kk=1 ; kk<NN.z-1 ; kk++) for(int  jj=1 ; jj<NN.y-1 ; jj++) for(int  ii=1 ; ii<NN.x-1 ; ii++)
		{
			if (type == ChargeDistribution && UUn[ii][jj][kk] == 0)
			{
/*
				rres =
				a * pphi(ii-1,jj,kk)	+ b * pphi(ii+1,jj,kk) + 
				c * pphi(ii,jj-1,kk)	+ d * pphi(ii,jj+1,kk) + 
				e * pphi(ii,jj,kk-1)	+ f * pphi(ii,jj,kk+1) + 
				g * pphi(ii,jj,kk)	- h(ii,jj,kk);
				pphi(ii,jj,kk) -= wwb*rres;
				eErrNum += (rres*rres);
*/
				rres =
					a * pphi[ii-1][jj][kk]	+ b * pphi[ii+1][jj][kk] + 
					c * pphi[ii][jj-1][kk]	+ d * pphi[ii][jj+1][kk] + 
					e * pphi[ii][jj][kk-1]	+ f * pphi[ii][jj][kk+1] + 
					g * pphi[ii][jj][kk]	- h[ii][jj][kk];
				pphi[ii][jj][kk] -= wwb*rres;
				eErrNum			 += (rres*rres);
/*				
				rres =
					a * pphi.pMatrix3d[ii-1][jj][kk]	+ b * pphi.pMatrix3d[ii+1][jj][kk] + 
					c * pphi.pMatrix3d[ii][jj-1][kk]	+ d * pphi.pMatrix3d[ii][jj+1][kk] + 
					e * pphi.pMatrix3d[ii][jj][kk-1]	+ f * pphi.pMatrix3d[ii][jj][kk+1] + 
					g * pphi.pMatrix3d[ii][jj][kk]		- h.pMatrix3d[ii][jj][kk];
						pphi.pMatrix3d[ii][jj][kk]		-= wwb*rres;
				eErrNum += (rres*rres);
*/			
			}
			else if (type == PotentialDistribution && UUn[ii][jj][kk] == 0)
			{
				rres = 
					a * pphi[ii-1][jj][kk]	+ b * pphi[ii+1][jj][kk] + 
					c * pphi[ii][jj-1][kk]	+ d * pphi[ii][jj+1][kk] + 
					e * pphi[ii][jj][kk-1]	+ f * pphi[ii][jj][kk+1] + 
					g * pphi[ii][jj][kk]	/*- 0*/;
				pphi[ii][jj][kk] -= wwb*rres;
				eErrNum += pow(rres,2);
			}
		};
		
		eErr = eErrNum/eErrDen;
		sstep++;
	};
	
	eendTime	= clock();
	rrunTime	+= (eendTime - sstartTime);

	if(sstep==MaxStep)
	{
		cout<<"*** Allowed computation time exceeded.\n";
		cout<<"*** Maximum allowed iteration step reached.\n";
		cout<<"*** Precision of the result : "<< eErr<<".\n";
	}
	
	// cout<<"Max Steps = "<<MaxStep<<endl;
	// printf("epsilon: %e\n",eErr);
	// printf("steps in SOR solver    : %d\n",sstep);
	// printf("Run time for SOR solver: %fs\n",(double)rrunTime/100);

	// UUn.fwrite("results/UUn.dat");	
}
/**************************************************************************************/