/* Trees.cpp */
#include "Trees.h"

/**************************************************************************************/
bool Tree::AddNewLink(ResGrid dd, SizeGrid NN,
				CMatrix3D& UUn, CMatrix3D& pphi,
				CriticalFields& EEc, VoltageDrops& VVd,
				const Point& IInitiationPoint, ListLink& EEstablishedLinks,
				bool iisBndXingPossible,  bool iisRsDeveloped,
				bool iisLinkXingPossible, bool iisChannelEquipotential)
{
	clock_t		startTime = clock();
	clock_t		endTime;
	clock_t		runTime;

	Point				sstart;
	Point				eend;
	Link				CCandidate;
	Link				CChosenLink;
	ListLink			LListOfCandidates;
	ListLink::iterator	it1;
	ListLink::iterator	it2;
	double				CCounterOfCandidates;
	double				SSumProba;
	double				rr;							// random number in [0,1]

	double				linkAlt;					// Altitude (in meters) of endpoint of current link.
	bool				flag = false;				// is given the value "true" if the
													// candidate is crossing an establ-
													// ished link. Assume no crossing
													// link at start.
	bool				isFlash3D = true;			// equal true for  a 3-D flash
													// equal false for a 1-D flash used for testing purposes
	int cpt=0;

	/* SAM variables. */
	double overreach;								// The amound by which a candidate exceeds the electric field
													//  propagation threshold.
	double max_overreach(0);						// max_overreach = max(overreach), over all candidates.
	ListLink::iterator it=EEstablishedLinks.begin();

	printf("..: Attempting to add new link...\n");

	/* Seed the random number generator. */
	srand(time(NULL) * 101);

	sstart = IInitiationPoint;

	if(isFlash3D == true) do{
	/**********************************************************************************/
	/* At this point we are sure the starting point is available for linking		  */
	/**********************************************************************************/
	for(eend.i = sstart.i-1 ; eend.i <= sstart.i+1 ; eend.i++)
	{
		for(eend.j = sstart.j-1 ; eend.j <= sstart.j+1 ; eend.j++)
		{
			for(eend.k = sstart.k-1 ; eend.k <= sstart.k+1 ; eend.k++)
			{
				if( eend.i >= 0 && eend.i <= NN.x-1 &&
					eend.j >= 0 && eend.j <= NN.y-1 &&
					eend.k >= 0 && eend.k <= NN.z-1 )
					if(UUn[eend.i][eend.j][eend.k] == 0)
					{
					/**********************************************************************************/
					/* At this point we are sure the ending point is available for linking			  */
					/**********************************************************************************/
						CCandidate.init(sstart,eend, 0, dd,pphi[sstart.i][sstart.j][sstart.k], pphi[eend.i][eend.j][eend.k]);
					/**********************************************************************************/
					/* Now we establish a candidate link and check if crossing any established link	  */
					/**********************************************************************************/
						if(CCandidate.efield >= EEc.positive[CCandidate.end.k] ||
						   CCandidate.efield <= EEc.negative[CCandidate.end.k])
						{
							if(iisLinkXingPossible == false)
							{
								flag = false;
								if(CCandidate.type == x || CCandidate.type == y || CCandidate.type == z)
									flag = false;

								if(CCandidate.type == xy)
									for(it1 = EEstablishedLinks.begin() ; it1 != EEstablishedLinks.end() ; it1++)
										if(CCandidate.start.k == (*it1).start.k && CCandidate.start.k == (*it1).end.k &&
										   (CCandidate.start.i + CCandidate.end.i) == ((*it1).start.i + (*it1).end.i) &&
										   (CCandidate.start.j + CCandidate.end.j) == ((*it1).start.j + (*it1).end.j) )
											flag = true;

								if(CCandidate.type == yz)
									for(it1 = EEstablishedLinks.begin() ; it1 != EEstablishedLinks.end() ; it1++)
										if(CCandidate.start.i == (*it1).start.i && CCandidate.start.i == (*it1).end.i &&
										   (CCandidate.start.j + CCandidate.end.j) == ((*it1).start.j + (*it1).end.j) &&
										   (CCandidate.start.k + CCandidate.end.k) == ((*it1).start.k + (*it1).end.k) )
											flag = true;

								if(CCandidate.type == xz)
									for(it1 = EEstablishedLinks.begin() ; it1 != EEstablishedLinks.end() ; it1++)
										if(CCandidate.start.j == (*it1).start.j && CCandidate.start.j == (*it1).end.j &&
										   (CCandidate.start.i + CCandidate.end.i) == ((*it1).start.i + (*it1).end.i) &&
										   (CCandidate.start.k +CCandidate. end.k) == ((*it1).start.k + (*it1).end.k) )
											flag = true;

								if(CCandidate.type == xyz)
									for(it1 = EEstablishedLinks.begin() ; it1 != EEstablishedLinks.end() ; it1++)
										if( (CCandidate.start.i + CCandidate.end.i) == ((*it1).start.i + (*it1).end.i) &&
											(CCandidate.start.j + CCandidate.end.j) == ((*it1).start.j + (*it1).end.j) &&
											(CCandidate.start.k + CCandidate.end.k) == ((*it1).start.k + (*it1).end.k) )
											flag = true;

								if(flag == false)
								{
								/**********************************************************************************/
								/* At this stage CCandidate contains a viable candidate for further propagation	  */
								/**********************************************************************************/

								/* SAM procedure.  Used to determine the greatest amount by which a
								 * candidate exceeds the propagation threshold.
								 */
									if(CCandidate.efield >= EEc.positive[CCandidate.end.k])
										overreach = 100 * (CCandidate.efield - EEc.positive[CCandidate.end.k])/EEc.positive[CCandidate.end.k];
									else
										overreach = 100 * (CCandidate.efield - EEc.negative[CCandidate.end.k])/EEc.negative[CCandidate.end.k];

									if(overreach > max_overreach)
										max_overreach = overreach;


									// This sets the voltage drop wrt. the initiation point //
									if(CCandidate.start == IInitiationPoint)
										CCandidate.deltaV =
											(CCandidate.efield>=0)*VVd.positive[CCandidate.end.k]*CCandidate.l+
											(CCandidate.efield<=0)*VVd.negative[CCandidate.end.k]*CCandidate.l;
									else
										CCandidate.deltaV = it->deltaV +
											(CCandidate.efield>=0)*VVd.positive[CCandidate.end.k]*CCandidate.l+
											(CCandidate.efield<=0)*VVd.negative[CCandidate.end.k]*CCandidate.l;
							/*		// This section has only been designed for testing purposes //
									if(CCandidate.start == IInitiationPoint)
										CCandidate.deltaV =
										(CCandidate.efield>=0)*1+
										(CCandidate.efield<=0)*-1;
									else
										CCandidate.deltaV = it->deltaV +
										(CCandidate.efield>=0)*1+
										(CCandidate.efield<=0)*-1;
							*/		// Add the Link to the list of candidates //
									LListOfCandidates.push_back(CCandidate);
									// Increment counter of candidates //
									cpt++;
								}
							}
							else
							{
								LListOfCandidates.push_back(CCandidate);
								cpt++;
							}
						}
					}
				}
			}
		}

		if(sstart==IInitiationPoint)
		{
			it = EEstablishedLinks.begin();
			sstart = it->end;
		}
		else
		{
			it++;
			sstart = it->end;
		}
	}while(it!=EEstablishedLinks.end());

	if(isFlash3D == false) do{
	/**********************************************************************************/
	/* At this point we are sure the starting point is available for linking		  */
	/**********************************************************************************/
		eend.i = sstart.i;
		eend.j = sstart.j;
		for(eend.k = sstart.k-1 ; eend.k <= sstart.k+1 ; eend.k++)
			if(	eend.k >= 0 && eend.k <= NN.z-1 )
				if(UUn[eend.i][eend.j][eend.k] == 0)
				{
					/**********************************************************************************/
					/* At this point we are sure the ending point is available for linking			  */
					/**********************************************************************************/
					CCandidate.init(sstart,eend, 0, dd,pphi[sstart.i][sstart.j][sstart.k], pphi[eend.i][eend.j][eend.k]);
					/**********************************************************************************/
					/* Now we establish a candidate link and check if crossing any established link	  */
					/**********************************************************************************/
						if(CCandidate.efield >= EEc.positive[CCandidate.end.k] ||
						   CCandidate.efield <= EEc.negative[CCandidate.end.k])
						{
							if(iisLinkXingPossible == false)
							{
								flag = false;
								if(flag == false)
								{
								/**********************************************************************************/
								/* At this stage CCandidate contains a viable candidate for further propagation	  */
								/**********************************************************************************/
									// This sets the voltage drop wrt. the initiation point //
									if(CCandidate.start == IInitiationPoint)
										CCandidate.deltaV =
											(CCandidate.efield>=0)*VVd.positive[CCandidate.end.k]*CCandidate.l+
											(CCandidate.efield<=0)*VVd.negative[CCandidate.end.k]*CCandidate.l;
									else
										CCandidate.deltaV = it->deltaV +
											(CCandidate.efield>=0)*VVd.positive[CCandidate.end.k]*CCandidate.l+
											(CCandidate.efield<=0)*VVd.negative[CCandidate.end.k]*CCandidate.l;
									// Add the Link to the list of candidates //
									LListOfCandidates.push_back(CCandidate);
									// Increment counter of candidates //
									cpt++;
								}
							}
							else
							{
								LListOfCandidates.push_back(CCandidate);
								cpt++;
							}
						}
					};

		if(sstart==IInitiationPoint)
		{
			it = EEstablishedLinks.begin();
			sstart = it->end;
		}
		else
		{
			it++;
			sstart = it->end;
		}
	}while(it!=EEstablishedLinks.end());

	if(max_overreach > ANOMALOUS_OVERREACH)
	{
		printf("AA: Anomalous maximum candidate overreach encountered.\n");
		Var::pSum = IO::openFile("summary.txt", "a");
		fprintf(Var::pSum, "AA: Anomalous maximum candidate overreach of %lf%% encountered.\n", max_overreach);
		fprintf(Var::pSum, "ii:\t Number of links added so far: %d.\n", Var::NumLinks);
		fclose(Var::pSum);

	}

	/**********************************************************************************/
	/* Derive probability of propagation for each candidate.						  */
	/**********************************************************************************/
	SSumProba = 0;
	CCounterOfCandidates = 0;
	for (it1 = LListOfCandidates.begin() ; it1 != LListOfCandidates.end() ; it1++)
	{
		if( (*it1).efield >= EEc.positive[(*it1).end.k] )
			(*it1).proba = pow(fabs((*it1).efield -EEc.positive[(*it1).end.k]),eta);
		else if ( (*it1).efield <= EEc.negative[(*it1).end.k] )
			(*it1).proba = pow(fabs((*it1).efield -EEc.negative[(*it1).end.k]),eta);
		else{
			cout<<"ee:\t This link should not exist. There is an error in the code!!!\n";
			Var::pSum = IO::openFile("summary.txt", "a");
			fprintf(Var::pSum, "ee:\t Bad link encountered.\n");
			fprintf(Var::pSum, "ii:\t\t Number of links added so far: %d.\n", Var::NumLinks);
			fprintf(Var::pSum, "xx:\t Program exiting after fatal error.\n");
			fclose(Var::pSum);
			exit(1);
		}
		SSumProba	+= (*it1).proba;
		CCounterOfCandidates++;
	}

	if(CCounterOfCandidates == 0)
	{
		endTime = clock();
		runTime = endTime-startTime;
		printf("ii:\t No more candidates.\n");
		cout<<"ii:\t Run time for Link addition: "<<(double)runTime/CLOCKS_PER_SEC<<" s."<<endl;
		Var::pSum = IO::openFile("summary.txt", "a");
		fprintf(Var::pSum, "ii:\t Simulation has run out of candidates.\n");
		fprintf(Var::pSum, "ii:\t\t Number of links added so far: %d.\n", Var::NumLinks);

		if(Var::maxAlt > ANOMALOUS_HEIGHT)
		{
			Var::curType = JET;
			printf("ii:\t Discharge exceeded anomalous height of %lf km.\n", ANOMALOUS_HEIGHT/1e3);
			printf("..:\t\t Classifying discharge as a jet.\n");
			fprintf(Var::pSum, "ii:\t Discharge exceeded anomalous height of %lf km.\n", ANOMALOUS_HEIGHT/1e3);
			fprintf(Var::pSum, "..:\t\t Classifying discharge as a jet.\n");
		}
		else
		{
			Var::curType = INTRA_CLOUD;
			printf("ii:\t Classifying discharge as intracloud.\n");
			fprintf(Var::pSum, "ii:\t Classifying discharge as intracloud.\n");
		}
		fclose(Var::pSum);

		return false;
	}

	/* Keep track of the number of candidates and maximum candidate overreach.
	 * Keeping track of the number of candidates was the motivation to define the
	 * 'ListInt' data type.
	 */
	cout<<"ii:\t Number of candidates: "<<CCounterOfCandidates<<"."<<endl;
	Var::NumberOfCandidates.push_back(CCounterOfCandidates);
	printf("ii:\t Maximum candidate overreach: %lf%%.\n", max_overreach);
	Var::MaximumCandidateOverreach.push_back(max_overreach);

	/**********************************************************************************/
	/* Convert probability to a	segment of length between 0 and 1, using the proba as */
	/* defined till now would increase the difficulty to choose a weightered by the   */
	/* probas. Hence we convert the proba to ease the choice of a link. Link i is	  */
	/* given the proba sum(p(j), j=1..i), consequently each proba is still between 0  */
	/* and 1.																		  */
	/* We randomly choose a point in this interval, we choose the link which proba is */
	/* the smallest greater than the random value previously defined				  */
	/**********************************************************************************/
	it2 = LListOfCandidates.begin();
	(*it2).proba /= SSumProba;
	for(it1=LListOfCandidates.begin() ; it2 != LListOfCandidates.end() ; it1++)
	{
		it2 = it1;
		it2++;
		(*it2).proba /= SSumProba;
		(*it2).proba += (*it1).proba;
	}

	rr = random()/(double)RAND_MAX;
	it1 = LListOfCandidates.begin();
	while(it1 != LListOfCandidates.end() && (*it1).proba <= rr) {it1++;};
	CChosenLink = *it1;
	LListOfCandidates.clear();

	/**********************************************************************************/
	/* The chosen link is always the one propagating in the highest electric field,   */
	/* (i.e., the one with the highest propability to propagate). This way, the		  */
	/* stochasticity is eleminated.													  */
	/**********************************************************************************/
/*
		rr = 0;
		it1 = LListOfCandidates.begin();
		for(it1 = LListOfCandidates.begin() ; it1 != LListOfCandidates.end() ; it1++)
		{
			if( it1->proba >= rr)
			{
				rr			= it1->proba;
				CChosenLink = *it1;
			};
		}
		LListOfCandidates.clear();
*/
	/**********************************************************************************/
	/* We finally fix the link in Un and add it to the list of established links and  */
	/* update phi.																	  */
	/**********************************************************************************/
	UUn[CChosenLink.end.i][CChosenLink.end.j][CChosenLink.end.k] =
		UUn[CChosenLink.start.i][CChosenLink.start.j][CChosenLink.start.k];
	if (iisChannelEquipotential == true)
	{
		// The two following lines are not important since potential is rederived
		// afterwards to account for potential drop and overall neutrality in command
		// pphi0 = fMinSearch(pphi0, QQchannelPlus, VVmin,VVmax , eepsilon,MMaxStep, pphi_cha,pphi_amb ,UUn, IInitiationPoint,EEstablishedLinks, dd,NN);
		// in main.cpp
//		pphi[CChosenLink.end.i][CChosenLink.end.j][CChosenLink.end.k] =
//		pphi[CChosenLink.start.i][CChosenLink.start.j][CChosenLink.start.k];
	}
	else if (iisChannelEquipotential == false)
	{
		pphi[CChosenLink.end.i][CChosenLink.end.j][CChosenLink.end.k] =
		pphi[CChosenLink.start.i][CChosenLink.start.j][CChosenLink.start.k] -
		( (CChosenLink.efield<0) * (EEc.negative[CChosenLink.end.k]) +
		  (CChosenLink.efield>0) * (EEc.positive[CChosenLink.end.k]) ) * CChosenLink.l;
	}
	EEstablishedLinks.push_back(CChosenLink);

	/**********************************************************************************/
	/* Send authorization of propagation if a link has been established and that it   */
	/* does not reach the upper/lower boundary.										  */
	/* If we have NN.z points, the matrix is indexed from 0 to NN.z-1, as we here     */
	/* consider that boundaries cannot be linked, the accessible range is 1..NN.z-2   */
	/**********************************************************************************/
	endTime = clock();
	runTime = endTime - startTime;

	/**********************************************************************************/
	/* Check boundary crossing														  */
	/**********************************************************************************/
	if(NN.IsOnBoundary(CChosenLink.end))
		if(iisBndXingPossible == false)
		{
			printf("ii:\t Channel reached a boundary at [%d, %d, %d].\n", CChosenLink.end.i, CChosenLink.end.j, CChosenLink.end.k);
			cout<<"ii:\t Run time for Link addition: "<<(double)runTime/CLOCKS_PER_SEC<<" s"<<endl;

			Var::pSum = IO::openFile("summary.txt", "a");
			fprintf(Var::pSum, "ii:\t Channel reached a boundary at [%d, %d, %d].\n", CChosenLink.end.i, CChosenLink.end.j, CChosenLink.end.k);
			fprintf(Var::pSum, "ii:\t\t Number of links added so far: %d.\n", Var::NumLinks);

			if(CChosenLink.end.k == (Var::N.z - 1))
			{
				Var::curType = JET;
				printf("ii:\t Discharge reached top of domain.\n");
				printf("..:\t\t Classifying discharge as a jet.\n");
				fprintf(Var::pSum, "ii:\t Discharge reached top of domain.\n");
				fprintf(Var::pSum, "..:\t\t Classifying discharge as a jet.\n");
			}
			else if(CChosenLink.end.k == 0)
			{
				Var::curType = CLOUD_TO_GROUND;
				printf("ii:\t Discharge reached bottom of domain.\n");
				printf("..:\t\t Classifying discharge as a cloud-to-ground strike.\n");
				fprintf(Var::pSum, "ii:\t Discharge reached bottom of domain.\n");
				fprintf(Var::pSum, "..:\t\t Classifying discharge as a cloud-to-ground strike.\n");
			}
			else
			{
				Var::curType = HORIZONTAL;
				printf("ii:\t Discharge reached a side of the domain.\n");
				printf(".,:\t\t Classifying discharge as \"horizontal\".\n");
				fprintf(Var::pSum, "ii:\t Discharge reached a side of the domain.\n");
				fprintf(Var::pSum, ".,:\t\t Classifying discharge as \"horizontal\".\n");
			}

			fclose(Var::pSum);

			return false;
		};

	/**********************************************************************************/
	/* Allow/Prevent Return Stroke DeVeloPmenT										  */
	/**********************************************************************************/
	linkAlt = CChosenLink.end.GetZ();
	printf("ii:\t Link termination position: [%3.1f %3.1f %3.1f] km\n",(CChosenLink.end.i*Var::d.x)/1e3, (CChosenLink.end.j*Var::d.y)/1e3,(CChosenLink.end.k*Var::d.z+Var::z_gnd)/1e3);
	if(linkAlt > Var::maxAlt)
		Var::maxAlt = linkAlt;

	if(linkAlt > ANOMALOUS_HEIGHT)
	{
		printf("AA:\t\t Altitude of link indicates possibility of jet development.\n");
		Var::pSum = IO::openFile("summary.txt", "a");
		fprintf(Var::pSum, "AA:\t\t Altitude of link indicates possibility of jet development.\n");
		fclose(Var::pSum);
	}

	if(CChosenLink.end.k == 0 || CChosenLink.end.k == NN.z-1)
	{
		if(iisRsDeveloped == true)
		{
			for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++)
			{
				// Turn channel potential to 0 //
				for(int kk=0 ; kk<NN.z ; kk++)
					if(UUn[ii][jj][kk] == UUn[CChosenLink.end.i][CChosenLink.end.j][CChosenLink.end.k])
						pphi[ii][jj][kk] = 0;
				// Prevent changes in ground/ionospheric potential //
				UUn[ii][jj][CChosenLink.end.k] = UUn[CChosenLink.end.i][CChosenLink.end.j][CChosenLink.end.k];
				pphi[ii][jj][CChosenLink.end.k] = 0;
			};
			cout<<"ii:\t Return Stroke Developed!\n";
			Var::pSum = IO::openFile("summary.txt", "a");
			fprintf(Var::pSum, "ii:\t Return stroke developed!\n");
			fprintf(Var::pSum, "ii:\t\t Number of links added so far: %d.\n", Var::NumLinks);
			fclose(Var::pSum);
		};
	}

	cout<<"ii:\t Run time for Link addition: "<<(double)runTime/CLOCKS_PER_SEC<<" s."<<endl;
	return true;
}

/**************************************************************************************/
/* Adjust potential in the channel													  */
/**************************************************************************************/
double Tree::Qchannel(const double VV,
					  const double eepsilon, const int MMaxStep,
					  CMatrix3D& pphi_cha, CMatrix3D& pphi_amb, CMatrix3D& UUn,
					  const Point& IInitiationPoint, ListLink& EEstablishedLinks,
					  ResGrid dd, const SizeGrid& NN)
{
	ListLink::iterator it;

	Potential	P1;//(pphi_cha,UUn);
	SorSolution	SSOR;//(pphi_cha, eepsilon,MMaxStep, dd, NN, P1, UUn);

	pphi_cha[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k]	= VV - pphi_amb[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k];
	for(it=EEstablishedLinks.begin() ; it!= EEstablishedLinks.end() ; it++)
		pphi_cha[it->end.i][it->end.j][it->end.k]	= VV - it->deltaV - pphi_amb[it->end.i][it->end.j][it->end.k];

	P1.init(pphi_cha,UUn);
	SSOR.init(pphi_cha, eepsilon,MMaxStep, dd, NN, P1, UUn);
	SSOR.Solve(dd,NN,UUn,pphi_cha);

	double QQ = foo::ChannelCharge(foo::Globalrho(pphi_cha,dd,NN),UUn,dd,NN);
	return QQ;
}
/**************************************************************************************/

/**************************************************************************************/
/* Derive Channel Potential to minimize total Charge - Dichotomy vs. Nelder-Mead	  */
/**************************************************************************************/
double Tree::fMinSearch(const double VV, const double QQchannelPlus,
						const double VVmin, const double VVmax,
						const double eepsilon, const int MMaxStep,
						CMatrix3D& pphi_cha, CMatrix3D& pphi_amb, CMatrix3D& UUn,
						const Point& IInitiationPoint, ListLink& EEstablishedLinks,
						ResGrid dd, const SizeGrid& NN)
{
	double choice = 0;
	if (choice == 0) // Bisection Method //
	{
		clock_t startTime = clock();
		double	CC;
		double	CCl	= VVmin;
		double	CCr = VVmax;
		double  CCav((CCr+CCl)/2);
		double	QQr = Qchannel( CCr , eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, IInitiationPoint,EEstablishedLinks, dd,NN);
		double	QQl = Qchannel( CCl , eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, IInitiationPoint,EEstablishedLinks, dd,NN);
		double	QQav= Qchannel( CCav, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, IInitiationPoint,EEstablishedLinks, dd,NN);
		int		kk	= 0;

		if(CCl > CCr)
		{
			Swap::DBL(CCl,CCr);
			Swap::DBL(QQl,QQr);
		}

//		while(fabs(2*(CCr-CCl)/(CCr+CCl))>eepsilon)
		while(fabs(QQav/QQchannelPlus)>eepsilon && fabs(2*(CCr-CCl)/(CCr+CCl))>eepsilon)
		{
			kk++;
			CCav = (CCr+CCl)/2;
			QQav = Qchannel( CCav, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, IInitiationPoint,EEstablishedLinks, dd,NN);
			if( QQl*QQav > 0 )
			{
				CCl = CCav;
				QQl = QQav;
			}
			else
			{
				CCr = CCav;
				QQr = QQav;
			}
//			cout<<"[Vl Vr] = ["<<setw(12)<<CCl<<" "<<setw(12)<<CCr<<"] ; [Ql Qr] = ["<<setw(12)<<QQl<<setw(13)<<QQr<<"] ; Nb of iterations = "<<kk<<endl;
//			cout<<"|QQav/QQchannelPlus| = "<<setw(12)<<fabs(QQav/QQchannelPlus)<<endl;
		}
//		if(fabs(QQav/QQchannelPlus)<=eepsilon)		cout<<"==============================>> break on CHARGE   "<<endl;
//		if(fabs(2*(CCr-CCl)/(CCr+CCl))<=eepsilon)	cout<<"==============================>> break on POTENTIAL"<<endl;

		CC		= (fabs(QQr)>fabs(QQl))*CCl+(fabs(QQl)>fabs(QQr))*CCr; //choose value which minimize QQ
//		QQav	= (fabs(QQr)>fabs(QQl))*QQl+(fabs(QQl)>fabs(QQr))*QQr;
		QQav = Qchannel( CCav, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, IInitiationPoint,EEstablishedLinks, dd,NN);
		clock_t endTime = clock();
		clock_t runTime = endTime - startTime;

//		double QQcha_cha = ChannelCharge(Globalrho(pphi_cha,dd,NN),UUn,dd,NN);
//		double QQcha_tot = TotalCharge(Globalrho(pphi_cha,dd,NN),dd,NN);
//		cout<<"QQcha_cha = "<<QQcha_cha<<endl;
//		cout<<"QQcha_tot = "<<QQcha_tot<<endl;
//		cout<<"Qdiff     = " << QQcha_tot-QQcha_cha<<endl;

		cout<<"ii:\t Run time for Qminimization: "<<(double)runTime/CLOCKS_PER_SEC<<"s\n";
//		cout<<"Vc = "<<setw(12)<<CC<<" ; Qc/Qc+ = "<<setw(12)<<QQav<<" ; Nb of iterations = "<<kk<<endl;
		return CC;
	}
	else // Nelder-Mead //
	{
		clock_t startTime = clock();
		double	CC;
		double	x1(VVmin)	, f1;
		double	x2(VVmax)	, f2;
//		double	x1((1-1e-5)*VV)	, f1;
//		double	x2((1+1e-5)*VV)	, f2;
		double	xr(0)		, fr(0);
		double	xe(0)		, fe(0);
		double	xc(0)		, fc(0);
		double	xcc(0)		, fcc(0);
		int		kk = 0;


		while(fabs(2*(x1-x2)/(x1+x2))>eepsilon /*&& kk<25*/)
//		while(fabs(Qchannel((x1+x2)/2, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, IInitiationPoint,EEstablishedLinks, dd,NN))>=1e-6)
		{
			kk++;
			if		(kk != 1 && x1 == xr )	{f1 = fr ;}
			else if	(kk != 1 && x1 == xe )	{f1 = fe ;}
			else if	(kk != 1 && x1 == xc )	{f1 = fc ;}
			else if	(kk != 1 && x1 == xcc)	{f1 = fcc;}
			else							{f1 = fabs(Qchannel(x1, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, IInitiationPoint,EEstablishedLinks, dd,NN));}

			if		(kk != 1 && x2 == xr )	{f2 = fr ;}
			else if	(kk != 1 && x2 == xe )	{f2 = fe ;}
			else if	(kk != 1 && x2 == xc )	{f2 = fc ;}
			else if	(kk != 1 && x2 == xcc)	{f2 = fcc;}
			else							{f2 = fabs(Qchannel(x2, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, IInitiationPoint,EEstablishedLinks, dd,NN));}
			// Step 1: Order //
			if(f1 > f2)
			{
				Swap::DBL(x1,x2);
				Swap::DBL(f1,f2);
			}
			// Step 2: Reflect //
			xr = 2.*x1-x2;
			if		(kk != 1 && xr == x1 )	{fr = f1 ;}
			else if	(kk != 1 && xr == x2 )	{fr = f2 ;}
			else if	(kk != 1 && xr == xe )	{fr = fe ;}
			else if	(kk != 1 && xr == xc )	{fr = fc ;}
			else if	(kk != 1 && xr == xcc)	{fr = fcc;}
			else							{fr = fabs(Qchannel(xr, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, IInitiationPoint,EEstablishedLinks, dd,NN));}
			// Step 3: Expand //
			if(fr < f1)
			{
				xe = 3.*x1-2.*x2;
				if		(kk != 1 && xe == x1 )	{fe = f1 ;}
				else if	(kk != 1 && xe == x2 )	{fe = f2 ;}
				else if	(kk != 1 && xe == xr )	{fe = fr ;}
				else if	(kk != 1 && xe == xc )	{fe = fc ;}
				else if	(kk != 1 && xe == xcc)	{fe = fcc;}
				else							{fe = fabs(Qchannel(xe, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, IInitiationPoint,EEstablishedLinks, dd,NN));}

				if(fe<fr) {x2 = xe;}
				else {x2 = xr;} // (fe>=fr)
			}
			// Step 4: Contract //
			else
			{
				// Step 4a) Outside //
				if(f1<=fr && fr<f2)
				{
					xc = 3/2.*x1-1/2.*x2;
					if		(kk != 1 && xc == x1 )	{fc = f1 ;}
					else if	(kk != 1 && xc == x2 )	{fc = f2 ;}
					else if	(kk != 1 && xc == xr )	{fc = fr ;}
					else if	(kk != 1 && xc == xe )	{fc = fe ;}
					else if	(kk != 1 && xc == xcc)	{fc = fcc;}
					else							{fc = fabs(Qchannel(xc, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, IInitiationPoint,EEstablishedLinks, dd,NN));}

					if(fc<=fr) {x2=xc;}
					// Step 5: Shrink Step //
					else {x2 = x1 + 1/2.*(x2-x1);}
				}
				// Step 4b) Inside //
				if(fr>=f2)
				{
					xcc = 1/2.*x1+1/2.*x2;
					if		(kk != 1 && xcc == x1)	{fcc = f1 ;}
					else if	(kk != 1 && xcc == x2)	{fcc = f2 ;}
					else if	(kk != 1 && xcc == xr)	{fcc = fr ;}
					else if	(kk != 1 && xcc == xe)	{fcc = fe ;}
					else if	(kk != 1 && xcc == xc)	{fcc = fc;}
					else							{fcc = fabs(Qchannel(xcc, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, IInitiationPoint,EEstablishedLinks, dd,NN));}

					if(fcc<f2) {x2=xcc;}
					// Step 5: Shrink Step //
					else {x2 = x1 + 1/2*(x2-x1);}
				}
			}
//			cout<<"[x1 x2] = ["<<setw(12)<<x1<<" "<<setw(12)<<x2<<"] ; Q = "<<setw(12)<<fabs(Qchannel((x1+x2)/2, eepsilon,MMaxStep, pphi_cha,pphi_amb,UUn, IInitiationPoint,EEstablishedLinks, dd,NN))<<" ; Nb of iterations = "<<kk<<endl;
		}
//		cout<<"*** FinaL ***"<<endl<<"[x1 x2] = ["<<setw(12)<<x1<<" "<<setw(12)<<x2<<"] ;\nQ = "<<setw(12)<<ChannelCharge(Globalrho(pphi_cha,dd,NN), UUn, dd,NN)<<endl;
//		cout<<"Qt = "<<setw(12)<<TotalCharge(Globalrho(pphi_cha+pphi_amb,dd,NN), dd,NN)<<endl;
		clock_t endTime = clock();
		clock_t runTime = endTime - startTime;

		cout<<"ii:\t Run time for Qminimization: "<<(double)runTime/CLOCKS_PER_SEC<<" s\n";
		cout<<"Vc = "<<setw(12)<<x1<<" ; Qc = "<<setw(12)<<f1<<" ; Nb of iterations = "<<kk<<endl;
		return CC=x1;
	}
}
/**************************************************************************************/

/**************************************************************************************/
double Tree::EqualizeAtGroundPotential(const double eepsilon, const int MMaxStep,
									   CMatrix3D& pphi_cha, CMatrix3D& pphi_amb, CMatrix3D& UUn,
									   const Point& IInitiationPoint, ListLink& EEstablishedLinks,
									   ResGrid dd, const SizeGrid& NN)
{
	ListLink::iterator	it;
	Potential			P1;		//(pphi_cha,UUn);
	SorSolution			SSOR;	//(pphi_cha, eepsilon,MMaxStep, dd, NN, P1, UUn);

	pphi_cha[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k]	= -pphi_amb[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k];
	for(it=EEstablishedLinks.begin() ; it!= EEstablishedLinks.end() ; it++)
		pphi_cha[it->end.i][it->end.j][it->end.k]	= -pphi_amb[it->end.i][it->end.j][it->end.k];
	P1.init(pphi_cha,UUn);
	SSOR.init(pphi_cha, eepsilon,MMaxStep, dd, NN, P1, UUn);
	SSOR.Solve(dd,NN,UUn,pphi_cha);
	return 0; // Total potential of the channel phi = phi_amb + phi_cha = phi_amb + (-phi_amb) = 0;
}
/**************************************************************************************/
