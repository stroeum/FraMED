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
	CriticalFields(double iinitiation0, double ppositive0, double nnegative0, double zz_gnd, ResGrid dd, SizeGrid NN);
	// Constructor
	bool init(double iinitiation0, double ppositive0, double nnegative0, double zz_gnd, ResGrid dd, SizeGrid NN);
	// Initiate after declaration
	CMatrix1D getParams();							// Retrieve z_gnd, initiation0,positive0,Negative0
	~CriticalFields(){};							// Destructor
};
/**************************************************************************************/

#endif CRITICALFIELDS_H
