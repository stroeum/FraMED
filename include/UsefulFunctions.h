/*
 *  UsefulFunctions.h
 *  Created by Jeremy Riousset on 10/25/07.
 */

#ifndef USEFULFUNCTIONS_H
#define USEFULFUNCTIONS_H
#include "Input.h"
#include "IOFunctions.h"
#include <algorithm> 
class foo
{
public:
	/******************************************************************************/
	/* For Electric Field														  */
	/******************************************************************************/
	static CMatrix1D	Eijk(int, int, int,const CMatrix3D&, const ResGrid&, const SizeGrid&);	// Derive local Ex, Ey, Ez at point i,j,k as well as its magnitude.
	static CMatrix3D	GlobalE(const CMatrix3D&,const ResGrid&, const SizeGrid&, const int);	// Derive magnitude of E everywhere
	static CMatrix1D	eFieldFlux(const CMatrix3D& phi, const ResGrid& d, const SizeGrid& N);	// Derive flux of E-field through the boundaries in C.
	/******************************************************************************/

	/******************************************************************************/
	/* For Charge																  */
	/******************************************************************************/
	static double		rhoijk(int, int, int,const CMatrix3D&, const ResGrid&, const SizeGrid&);	// Derive the charge distribution at point i,j,k
	static CMatrix3D	Globalrho(const CMatrix3D&,const ResGrid&, const SizeGrid&);				// Derive the charge distribution
	static double		ChannelCharge(const CMatrix3D& rho, const CMatrix3D& Un, const ResGrid& d, const SizeGrid& N);
	static double		ChannelChargePositive(const CMatrix3D& rho, const CMatrix3D& Un, const ResGrid& d, const SizeGrid& N);
	static double		ChannelChargeNegative(const CMatrix3D& rho, const CMatrix3D& Un, const ResGrid& d, const SizeGrid& N);
	static double		TotalCharge(const CMatrix3D& rho, const ResGrid& d, const SizeGrid& N);
	static CMatrix1D	ChannelLinearDensity(const CMatrix3D& rho, const CMatrix3D& Un, const ResGrid& d, const SizeGrid& N);
	/******************************************************************************/

	/******************************************************************************/
	/* For Dipole Moment														  */
	/******************************************************************************/
	static Vector		DipoleMoment(double& CCarriedCharge, const CMatrix3D& phi, const CMatrix3D& Un, const SizeDomain& L, const ResGrid& d, const SizeGrid& N);
	// Also derive charge carried by the channel
	/******************************************************************************/

	/******************************************************************************/
	/* For E, phi, rho															  */
	/******************************************************************************/
	static bool		isfinite(const CMatrix3D& MM, const SizeGrid& N);				// Check finiteness of the matrix
	/******************************************************************************/
};
#endif
