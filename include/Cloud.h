/* Cloud.h */

#ifndef CLOUD_H
#define CLOUD_H

#include "IOFunctions.h"
#include "UsefulFunctions.h"
#include "BoundaryConditions.h"
#include "Swap.h"

#define ANOMALOUS_HEIGHT (11e3)
#define ANOMALOUS_OVERREACH (200.0)

class Cloud
{
public:
    static void	LoadTripole(FILE *, const double I1, const double I2, const double Iscreen);
    static double EoverEk(const double t, CriticalFields OverShotEc);

};

#endif
