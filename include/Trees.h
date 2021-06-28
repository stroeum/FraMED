/* Trees.h */

#ifndef TREES_H
#define TREES_H

#include "IOFunctions.h"
#include "UsefulFunctions.h"
#include "BoundaryConditions.h"
#include "Swap.h"

#define ANOMALOUS_HEIGHT (11e3)
#define ANOMALOUS_OVERREACH (200.0)


class Tree
{
public:
    /******************************************************************************/
    static bool init(char *  folder, char * file, ListLink& EstablishedLinks);     // initiate instance of class

    /******************************************************************************/
    static bool Initiate(FILE * file, const int InitiationType, Point& InitiationPoint);       // find initiation point

    /******************************************************************************/
    static void Grow(FILE * file, bool AAddNew);
    
    /******************************************************************************/
    static void StoreData(FILE * file); /* Writes intermediate as well as final data to files. */

    /******************************************************************************/
	static bool AddNewLink(FILE * file, ResGrid d, SizeGrid N,
						   CMatrix3D& Un, CMatrix3D& phi,
						   CriticalFields& Ec, VoltageDrops& Vd,
						   const Point& InitiationPoint, ListLink& EstablishedLinks,
						   bool isBndXingPossible, bool isRsDeveloped,
						   bool isLinkXingPossible, bool isChannelEquipotential);

	/******************************************************************************/

	static double fMinSearch(FILE * file, const double V, const double QQchannelPlus,
							 const double Vmin, const double Vmax,
							 const double epsilon, const int MaxStep,
							 CMatrix3D& phi_cha, CMatrix3D& phi_amb, CMatrix3D& Un,
							 const Point& InitiationPoint, ListLink& EstablishedLinks,
							 ResGrid d, const SizeGrid& N);

	/******************************************************************************/
	static double EqualizeAtGroundPotential(const double epsilon, const int MaxStep,
											CMatrix3D& phi_cha, CMatrix3D& phi_amb, CMatrix3D& Un,
											const Point& InitiationPoint, ListLink& EstablishedLinks,
											ResGrid d, const SizeGrid& N);
	/******************************************************************************/

	static double Qchannel(const double V,
						   const double epsilon, const int MaxStep,
						   CMatrix3D& phi_cha, CMatrix3D& phi_amb, CMatrix3D& Un,
						   const Point& InitiationPoint, ListLink& EstablishedLinks,
						   ResGrid d, const SizeGrid& N);
	/******************************************************************************/
};

#endif
