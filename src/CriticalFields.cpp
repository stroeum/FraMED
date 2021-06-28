/*
 *  CriticalFields.cpp
 *  Created by Jeremy Riousset on 10/25/07.
 */

#include "CriticalFields.h"

/**************************************************************************************/
CriticalFields::CriticalFields(double iinitiation0, double ppositive0, double nnegative0, double zz_gnd, ResGrid _d, SizeGrid _N, int ScalingExponent)
{CriticalFields::init(iinitiation0,ppositive0,nnegative0, zz_gnd, _d,_N, ScalingExponent);}

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
						  double zz_gnd, ResGrid _d, SizeGrid _N, int ScalingExponent)
{
	z_gnd		= zz_gnd;
	initiation0 = iinitiation0;
	positive0	= ppositive0;
	negative0	= nnegative0;
	initiation.init(_N.z);
	positive.init(_N.z);
	negative.init(_N.z);
    
    int alpha = ScalingExponent;

	int kk_gnd = (int)round(z_gnd/_d.z);
	for (int kk=0 ; kk<_N.z ; kk++)
	{
		initiation[kk] = initiation0*pow(Scaling::StdAtmosphere((kk+kk_gnd)*_d.z),alpha);
		positive[kk] = positive0*pow(Scaling::StdAtmosphere((kk+kk_gnd)*_d.z),alpha);
		negative[kk] = negative0*pow(Scaling::StdAtmosphere((kk+kk_gnd)*_d.z),alpha);
	}
	return true;
}
/**************************************************************************************/
