/*
 *  CriticalFields.cpp
 *  Created by Jeremy Riousset on 10/25/07.
 */

#include "CriticalFields.h"

/**************************************************************************************/
CriticalFields::CriticalFields(double iinitiation0, double ppositive0, double nnegative0, double zz_gnd, ResGrid dd, SizeGrid NN)
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
						  double zz_gnd, ResGrid dd, SizeGrid NN)
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
		initiation[kk] = initiation0*Scaling::StdAtmosphere((kk+kk_gnd)*dd.z);
		positive[kk] = positive0*Scaling::StdAtmosphere((kk+kk_gnd)*dd.z);
		negative[kk] = negative0*Scaling::StdAtmosphere((kk+kk_gnd)*dd.z);
	}
	return true;
}
/**************************************************************************************/
