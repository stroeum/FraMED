/* BoundaryConditions.h */

#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include "IOFunctions.h"
#include "UsefulFunctions.h"
class BC
{
public:
	static void Apply(int BCtype, CMatrix3D& phi, CMatrix3D rho, const ResGrid& d, const SizeGrid& N);
	static void Update(bool isFlashAccoutedForInBC, int BCtype, CMatrix3D& phi_cha, const double rhoAmbMin, const double rhoAmbMax, const ResGrid& d, const SizeGrid& N);
	static void AddUniformE(bool UniformE, const double& Eo, CMatrix3D& phi, CMatrix3D& Un, const SizeDomain& L, const ResGrid& d, const SizeGrid& N);
};
#endif
