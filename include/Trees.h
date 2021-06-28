/* Trees.h */

#ifndef TREES_H
#define TREES_H

#include "IOFunctions.h"
#include "UsefulFunctions.h"
#include "Swap.h"

#define ANOMALOUS_HEIGHT (11e3)
#define ANOMALOUS_OVERREACH (200.0)
class Tree
{
public:
	/******************************************************************************/
	static bool AddNewLink(ResGrid dd, SizeGrid NN,
						   CMatrix3D& UUn, CMatrix3D& pphi,
						   CriticalFields& EEc, VoltageDrops& VVd,
						   const Point& IInitiationPoint, ListLink& EEstablishedLinks,
						   bool iisBndXingPossible, bool iisRsDeveloped,
						   bool iisLinkXingPossible, bool iisChannelEquipotential);

	/******************************************************************************/

	static double fMinSearch(const double VV, const double QQchannelPlus,
							 const double VVmin, const double VVmax,
							 const double eepsilon, const int MMaxStep,
							 CMatrix3D& pphi_cha, CMatrix3D& pphi_amb, CMatrix3D& UUn,
							 const Point& IInitiationPoint, ListLink& EEstablishedLinks,
							 ResGrid dd, const SizeGrid& NN);

	/******************************************************************************/
	static double EqualizeAtGroundPotential(const double eepsilon, const int MMaxStep,
											CMatrix3D& pphi_cha, CMatrix3D& pphi_amb, CMatrix3D& UUn,
											const Point& IInitiationPoint, ListLink& EEstablishedLinks,
											ResGrid dd, const SizeGrid& NN);
	/******************************************************************************/

	static double Qchannel(const double VV,
						   const double eepsilon, const int MMaxStep,
						   CMatrix3D& pphi_cha, CMatrix3D& pphi_amb, CMatrix3D& UUn,
						   const Point& IInitiationPoint, ListLink& EEstablishedLinks,
						   ResGrid dd, const SizeGrid& NN);
	/******************************************************************************/
};

#endif TREES_H
