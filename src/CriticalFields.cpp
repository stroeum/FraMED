/* CriticalFields.cpp */

#include "CriticalFields.h"

/**************************************************************************************/
ScaledFields::ScaledFields(double ppositive0, double nnegative0,
						   double zz_gnd, StepsSizes dd, BoxSteps NN)
{ScaledFields(ppositive0,nnegative0, zz_gnd, dd,NN);}

CMatrix1D ScaledFields::getParams()
{
	CMatrix1D PParams(3);
	PParams[0] = positive0;
	PParams[1] = negative0;
	PParams[2] = z_gnd;
	return PParams;
}

bool ScaledFields::init(double ppositive0, double nnegative0,
						double zz_gnd, StepsSizes dd, BoxSteps NN)
{
	z_gnd		= zz_gnd;
	positive0 = ppositive0;
	negative0 = nnegative0;
	positive.init(NN.z);
	negative.init(NN.z);
	
	int kk_gnd = (int)round(z_gnd/dd.z);
	for (int kk=0 ; kk<NN.z ; kk++)
	{
		positive[kk] = positive0*StdAtmosphere((kk+kk_gnd)*dd.z);
		negative[kk] = negative0*StdAtmosphere((kk+kk_gnd)*dd.z);
	}
	return true;
}
/**************************************************************************************/

/**************************************************************************************/
CriticalFields::CriticalFields(double iinitiation0, 
							   double ppositive0, double nnegative0,
							   double zz_gnd, StepsSizes dd, BoxSteps NN)
{CriticalFields::init(iinitiation0,ppositive0,nnegative0, zz_gnd, dd,NN);}

CMatrix1D CriticalFields::getParams()
{
	CMatrix1D PParams(4);
	PParams[0] = initiation0;
	PParams[1] = positive0;
	PParams[2] = negative0;
	PParams[3] = z_gnd;
	return PParams;
}

bool CriticalFields::init(double iinitiation0, double ppositive0, double nnegative0,
		  double zz_gnd, StepsSizes dd, BoxSteps NN)
{
	z_gnd		= zz_gnd;
	initiation0 = iinitiation0;
	positive0	= ppositive0;
	negative0	= nnegative0;
	initiation.init(NN.z);
	positive.init(NN.z);
	negative.init(NN.z);
	
	int kk_gnd = (int)round(z_gnd/dd.z);
	for (int kk=0 ; kk<NN.z ; kk++)
	{
		initiation[kk] = initiation0*StdAtmosphere((kk+kk_gnd)*dd.z);
		positive[kk] = positive0*StdAtmosphere((kk+kk_gnd)*dd.z);
		negative[kk] = negative0*StdAtmosphere((kk+kk_gnd)*dd.z);
	}
	return true;
}
/**************************************************************************************/

/**************************************************************************************/
double StdAtmosphere(double hh)
{
	/**********************************************************************************/
	/*							!!! Altitude in km !!!								  */
	/**********************************************************************************/
	hh*=1e-3;

	double ddh,rresult;
	int uupperBound, llowerBound, mmedian;
	
	/**********************************************************************************/
	/* US Standard Atmosphere:														  */
	/*	We replaced original 2.5e19cm-3 at 0 km altitude by 2.688e19 cm-3 which is our*/
	/*	reference number density at temperature 273 K.								  */
	/* This will be our altitude scaling. We prefer this to the exp(-z/8.4km) found in*/
	/* litertature, because the last one is temperature dependent.					  */
	/**********************************************************************************/
	
	double AAltitude[]={
		0.e+0,		5.e+0,		1.e+1,		1.5e+1,		2.e+1,
		2.5e+1,		3.e+1,		3.5e+1,		4.e+1,		4.5e+1,
		5e+1,		5.5e+1,		6e+1,		6.4e+1,		6.8e+1,
		7.2e+1,		7.6e+1,		8e+1,		8.4e+1,		8.8e+1,
		9.2e+1,		9.6e+1,		1e+2,		1.08e+2,	1.14e+2,
		1.2e+2,		1.26e+2,	1.32e+2,	1.4e+2,		1.5e+2};	
	
	double NNeutralDensity[]={
		2.688e+19,	1.53e+19,	8.59e+18,	4.05e+18,	1.85e+18,
		8.33e+17,	3.83e+17,	1.76e+17,	8.31e+16,	4.088e+16,
		2.13e+16,	1.181e+16,	6.33e+15,	3.93e+15,	2.39e+15,
		1.39e+15,	7.72e+14,	4.03e+14,	1.99e+14,	9.48e+13,
		4.37e+13,	2.07e+13,	1.04e+13,	3.18e+12,	1.43e+12,
		6.61e+11,	3.4e+11,	1.91e+11,	9.7e+10,	4.92e+10};
	
	if (hh<0 || hh>150e3)
	{
		cout<<"h is outside of 0-150 km interval"<<endl;
		exit(2);
	}
	
	/**********************************************************************************/
	/* We have 30 samples partioned in 30 intervals. First, we derive the upper/lower */
	/* bounds (uupperBound/llowerBound) of the interval								  */
	/**********************************************************************************/
		
	llowerBound=0;
	uupperBound=29;
	while(uupperBound!=llowerBound+1)
	{
		mmedian=(uupperBound+llowerBound)/2;
		ddh=hh-AAltitude[mmedian];
		if(ddh<0)  uupperBound=mmedian;
		if(ddh>=0) llowerBound=mmedian;
	}
	
	/**********************************************************************************/
	/* Now, we do an interpolation to estimate the neutral density at any height	  */
	/* between 0 and 150 km															  */
	/**********************************************************************************/

	rresult = (NNeutralDensity[llowerBound]/NNeutralDensity[0])*
		pow(NNeutralDensity[uupperBound]/NNeutralDensity[llowerBound],
			(hh-AAltitude[llowerBound])/(AAltitude[uupperBound]-AAltitude[llowerBound]));
//	double result=pow(10.,log10(NNeutralDensity[llowerBound]/NNeutralDensity[0]) + 
//			  (log10(NNeutralDensity[uupperBound])-log10(NNeutralDensity[llowerBound]))
//			  *(hh-AAltitude[llowerBound])/(AAltitude[uupperBound]-AAltitude[llowerBound]));
//	double c = 100*abs(rresult-result)/rresult;
	return rresult;
	
	/**********************************************************************************/
	/* We note that result == rresult analytically. However it does not seem to be the*/
	/* case here. The tests we tried with Matlab, show the same difference. We guess  */
	/* either an error in our reasoning or an approximation due to the use of log10 in*/
	/* the first case. We choose the first solution.								  */
	/**********************************************************************************/
}
/**************************************************************************************/
