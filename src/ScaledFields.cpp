/*
 *  ScaledFields.cpp
 *  Created by Jeremy Riousset on 10/25/07.
 */

#include "ScaledFields.h"

/**************************************************************************************/
ScaledFields::ScaledFields(double ppositive0, double nnegative0, double zz_gnd, ResGrid _d, SizeGrid _N, int ScalingExponent)
{ScaledFields(ppositive0,nnegative0, zz_gnd, _d,_N, ScalingExponent);}
ScaledFields::ScaledFields(double ppositive0, double nnegative0, double zz_gnd, ResGrid _d, SizeGrid _N, int ScalingExponent, list<double> alt, list<double> ng)
{ScaledFields(ppositive0,nnegative0, zz_gnd, _d,_N, ScalingExponent,alt,ng);}
ScaledFields::ScaledFields(double ppositive0, double nnegative0, double zz_gnd, double H, ResGrid _d, SizeGrid _N)
{ScaledFields(ppositive0, nnegative0, zz_gnd, H, _d, _N);}

CMatrix1D ScaledFields::getParams()
{
	CMatrix1D PParams(3);
	PParams[0] = positive0;
	PParams[1] = negative0;
	PParams[2] = z_gnd;
	return PParams;
}

bool ScaledFields::init(double ppositive0, double nnegative0,
						double zz_gnd, ResGrid _d, SizeGrid _N, int ScalingExponent)
{
	z_gnd	  = zz_gnd;
	positive0 = ppositive0;
	negative0 = nnegative0;
	positive.init(_N.z);
	negative.init(_N.z);
    int alpha = ScalingExponent;

	int kk_gnd = (int)round(z_gnd/_d.z);
	for (int kk=0 ; kk<_N.z ; kk++)
	{
		positive[kk] = positive0*pow(Scaling::StdAtmosphere((kk+kk_gnd)*_d.z),alpha);
		negative[kk] = negative0*pow(Scaling::StdAtmosphere((kk+kk_gnd)*_d.z),alpha);
	}
	return true;
}

bool ScaledFields::init(double ppositive0, double nnegative0, double zz_gnd, ResGrid _d, SizeGrid _N, int ScalingExponent, list<double> alt, list<double> ng)
{
    z_gnd     = zz_gnd;
    positive0 = ppositive0;
    negative0 = nnegative0;
    positive.init(_N.z);
    negative.init(_N.z);
    int alpha = ScalingExponent;

    int kk_gnd = (int)round(z_gnd/_d.z);
    for (int kk=0 ; kk<_N.z ; kk++)
    {
        positive[kk] = positive0*pow(Scaling::StdAtmosphere((kk+kk_gnd)*_d.z,alt,ng),alpha);
        negative[kk] = negative0*pow(Scaling::StdAtmosphere((kk+kk_gnd)*_d.z,alt,ng),alpha);
    }
    return true;
}
/**************************************************************************************/

bool ScaledFields::init(double ppositive0, double nnegative0, double zz_gnd, double H, ResGrid _d, SizeGrid _N)
{ // Linear scaling !!! for testing purposes only !!!
	z_gnd		= zz_gnd;
	positive0	= ppositive0;
	negative0	= nnegative0;
	positive.init(_N.z);
	negative.init(_N.z);
	
	int kk_gnd = (int)round(z_gnd/_d.z);
	for (int kk=0 ; kk<_N.z ; kk++)
	{
		positive[kk] = positive0*(kk+kk_gnd)*_d.z/H;
		negative[kk] = negative0*(kk+kk_gnd)*_d.z/H;
	}
	return true;
}
