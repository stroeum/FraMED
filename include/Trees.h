/* Trees.h */

#ifndef TREES_H
#define TREES_H

#include <list>
#include <time.h>

#include "Constants.h"
#include "CriticalFields.h"
#include "Links.h"
#include "Matrix.h"
#include "OutputFunctions.h"
#include "SimParameters.h"
#include "Sources.h"
#include "SorSolution.h"

/**************************************************************************************/
bool AddNewLink(StepsSizes dd, BoxSteps NN, 
				CMatrix3D& UUn, CMatrix3D& pphi, 
				CriticalFields& EEc, VoltageDrops& VVd,
				const Point& IInitiationPoint, ListLink& EEstablishedLinks, 
				bool iisBndXingPossible, bool iisRsDeveloped, 
				bool iisLinkXingPossible, bool iisChannelEquipotential);
/**************************************************************************************/

/**************************************************************************************/
void write(ListLink&, string);
/**************************************************************************************/

/**************************************************************************************/
double fMinSearch(const double VV, const double QQchannelPlus,
				  const double VVmin, const double VVmax, 
				  const double eepsilon, const int MMaxStep,
				  CMatrix3D& pphi_cha, CMatrix3D& pphi_amb, CMatrix3D& UUn,
				  const Point& IInitiationPoint, ListLink& EEstablishedLinks,
				  StepsSizes dd, const BoxSteps& NN);

double EqualizeAtGroundPotential(const double eepsilon, const int MMaxStep,
								 CMatrix3D& pphi_cha, CMatrix3D& pphi_amb, CMatrix3D& UUn,
								 const Point& IInitiationPoint, ListLink& EEstablishedLinks,
								 StepsSizes dd, const BoxSteps& NN);
/**************************************************************************************/

/**************************************************************************************/
double Qchannel(const double VV,
				const double eepsilon, const int MMaxStep,
				CMatrix3D& pphi_cha, CMatrix3D& pphi_amb, CMatrix3D& UUn,
				const Point& IInitiationPoint, ListLink& EEstablishedLinks,
				StepsSizes dd, const BoxSteps& NN);
/**************************************************************************************/

#endif TREES_H