/*
 *  ScaledFields.h
 *  Created by Jeremy Riousset on 10/25/07.
 */

#ifndef SCALED_FIELDS_H
#define SCALED_FIELDS_H

#include "ResGrid.h"
#include "Matrix.h"
#include "AtmScaling.h"

/**************************************************************************************/
class ScaledFields
{
protected:
	double z_gnd;   								// Altitude of the ground plane
	double positive0;								// Value of the propagation positive threshold at z = z_gnd
	double negative0;								// Value of the propagation negative threshold at z = z_gnd
public:
    CMatrix1D positive;								// Positive propagation threshold
	CMatrix1D negative;								// Negative propagation threshold

	ScaledFields(){};								// Default constructor
	ScaledFields(double plus0, double minus0, double gnd_alt, ResGrid d, SizeGrid N, int ScalingExponent);
	// Constructor surcharge
	CMatrix1D getParams();							// Retrieve z_gnd, initiation0,positive0,Negative0
	bool init(double plus0, double minus0, double gnd_alt, ResGrid d, SizeGrid N, int ScalingExponent);
	// Initiate after declaration
	~ScaledFields(){};								// Destructor
};
/**************************************************************************************/

#endif
