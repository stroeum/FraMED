/* BoundaryConditions.h */

#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include <time.h>

#include "Constants.h"
#include "OutputFunctions.h"
#include "Matrix.h"
#include "SimParameters.h"

void ApplyBC(int BBCtype, CMatrix3D& pphi, CMatrix3D rrho, const StepsSizes& dd, const BoxSteps& NN);
void UpdateBC(bool iisFlashAccoutedForInBC, int BBCtype, CMatrix3D& pphi_cha, const double rrhoAmbMin, const double rrhoAmbMax, const StepsSizes& dd, const BoxSteps& NN);
void AddUniformE(bool uuniformE, const double& EEo, CMatrix3D& pphi, CMatrix3D& UUn, const BoxLengths& LL, const StepsSizes& dd, const BoxSteps& NN);
#endif //BOUNDARYCONDITIONS_H