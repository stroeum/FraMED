/*
 *  VoltageDrops.h
 *  Created by Jeremy Riousset on 10/25/07.
 */

#ifndef VOLTAGEDROPS_H
#define VOLTAGEDROPS_H
#include "ScaledFields.h"

/**************************************************************************************/
class VoltageDrops : public ScaledFields
{
public:
	VoltageDrops(){};
	VoltageDrops(double plus0, double minus0, double gnd_alt, ResGrid d, SizeGrid N) :
		ScaledFields(plus0,minus0,gnd_alt,d,N){};
	// Constructor
	~VoltageDrops(){};								// Destructor
};
/**************************************************************************************/

#endif
