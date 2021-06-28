/*
 *  CriticalFields.h
 *  Created by Jeremy Riousset on 10/25/07.
 */

#ifndef CRITICALFIELDS_H
#define CRITICALFIELDS_H
#include "ScaledFields.h"

/**************************************************************************************/
class CriticalFields : public ScaledFields
{
protected:
	double initiation0;								// Value of the initiation field at z = z_gnd
public:
	CMatrix1D initiation;							// Initiation Threshold
	CriticalFields(){};
	CriticalFields(double initiation0, double positive0, double negative0, double z_gnd, ResGrid d, SizeGrid N, int ScalingExponent);
	// Constructor
	bool init(double initiation0, double positive0, double negative0, double z_gnd, ResGrid d, SizeGrid N, int ScalingExponent);
	// Initiate after declaration
	CMatrix1D getParams();							// Retrieve z_gnd, initiation0,positive0,Negative0
	~CriticalFields(){};							// Destructor
};
/**************************************************************************************/

#endif
