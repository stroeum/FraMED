/* Trees.cpp */
#include "Trees.h"

/**************************************************************************************/
/* List of links																	  */
/**************************************************************************************/
void write(ListLink& LL, string ss)
{
	char *fname = &ss[0];
	
	ListLink::iterator	it;
	FILE *file;
	char *relativePath;
	
	/**********************************************************************************/
	/* Open the file to be writen.													  */
	/**********************************************************************************/
#ifdef __MWERKS__
	char *flags  = "w";
	int hierarchy;
	char MacPath[256],*Mname,*Uname;
	/* Parse the file name for the Macintosh */
	hierarchy=0;
	MacPath[0]=':';
	Mname=MacPath+1;
	for(Uname=fname;*Uname!='\0';Uname++){
		if(*Uname=='/'){
			*(Mname++) =':';
			hierarchy++;
		}else if(*Uname==':')
			nrerror("Sorry, ':' is not allowed in filenames");
		else
			*(Mname++)=*Uname;
	}
	*Mname='\0';
	if(hierarchy==0)
		relativePath=MacPath+1;
	else
		relativePath=MacPath;
#else /* UNIX */
	relativePath=fname;
#endif
	/* printf("Opening file '%s'\n",relativePath); */
	if(!(file=fopen(relativePath,"w")))
		cout<<"cannot open file";
	//nrerror("cannot open file"); 
	
	/**********************************************************************************/
	/* Write datas.																	  */
	/**********************************************************************************/
	for (it=LL.begin() ; it!=LL.end() ; it++)
	{
		fprintf(file,"%d %d %d %d %d %d %f %f %f %f\n", it->start.i, it->start.j, it->start.k, it->end.i, it->end.j, it->end.k, it->l, it->efield, it->deltaV, it->proba ); 
	}
	fclose(file);
}
/**************************************************************************************/

/**************************************************************************************/
bool AddNewLink(StepsSizes dd, BoxSteps NN, 
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
	bool				flag = false;				// is given the value "true" if the
													// candidate is crossing an establ-
													// ished link. Assume no crossing
													// link at start.
	bool				isFlash3D = true;			// equal true for  a 3-D flash
													// equal false for a 1-D flash used for testing purposes 
	int cpt=0;
	ListLink::iterator it=EEstablishedLinks.begin();
	
	sstart = IInitiationPoint;
	
	if(isFlash3D == true) do{
	/**********************************************************************************/
	/* At this point we are sure the starting point is available for linking		  */
	/**********************************************************************************/
	for(eend.i = sstart.i-1 ; eend.i <= sstart.i+1 ; eend.i++)
		for(eend.j = sstart.j-1 ; eend.j <= sstart.j+1 ; eend.j++)
			for(eend.k = sstart.k-1 ; eend.k <= sstart.k+1 ; eend.k++)
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
	
	/**********************************************************************************/
	/* Derive probability of probagation for each candidate.						  */
	/**********************************************************************************/		
	SSumProba = 0;
	CCounterOfCandidates = 0;
	for (it1 = LListOfCandidates.begin() ; it1 != LListOfCandidates.end() ; it1++)
	{
		if( (*it1).efield >= EEc.positive[(*it1).end.k] )
			(*it1).proba = pow(fabs((*it1).efield -EEc.positive[(*it1).end.k]),PMC.eta);	
		else if ( (*it1).efield <= EEc.negative[(*it1).end.k] )
			(*it1).proba = pow(fabs((*it1).efield -EEc.negative[(*it1).end.k]),PMC.eta);
		else{
			cout<<"!!! This link should not exist. There is an error in the code !!!\n";
			exit(1);
		}
		SSumProba	+= (*it1).proba;
		CCounterOfCandidates++;
	}

	if(CCounterOfCandidates == 0)
	{
		endTime = clock();
		runTime = endTime-startTime;
		cout<<"------> No more candidate"<<endl;
		cout<<"------> Run time for Link addition: "<<(double)runTime/CLOCKS_PER_SEC<<"s"<<endl;
//        cout<<"------> Run time for Link addition: "<<(double)runTime/100<<"s"<<endl;
		return false;
	}
	cout<<"------> Number of candidates      : "<<CCounterOfCandidates<<endl;
		
	/**********************************************************************************/
	/* Convert probabilite to a	segment of length between 0 and 1, using the proba as */
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
	
	rr = rand()/(double)RAND_MAX;
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
	if( CChosenLink.end.i == 0 || CChosenLink.end.i == NN.x-1 ||
		CChosenLink.end.j == 0 || CChosenLink.end.j == NN.y-1 ||
		CChosenLink.end.k == 0 || CChosenLink.end.k == NN.z-2)
		if(iisBndXingPossible == false)
		{
			cout<<"Channel reached a boundary (near a boundary if Nz)"<<endl;
			cout<<"------> Run time for Link addition: "<<(double)runTime/CLOCKS_PER_SEC<<"s"<<endl;
//            cout<<"------> Run time for Link addition: "<<(double)runTime/100<<"s"<<endl;
			return false;
		};
	
	/**********************************************************************************/
	/* Allow/Prevent Return Stroke DeVeloPmenT										  */
	/**********************************************************************************/
	cout<<"------> Link termination position : ["<<CChosenLink.end.i*dd.x*1e-3<<" "<<CChosenLink.end.j*dd.y*1e-3<<" "<<CChosenLink.end.k*dd.z*1e-3<<"]"<<endl;
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
			cout<<"! Return Stroke Developped !\n";
		};
	}
	
	cout<<"------> Run time for Link addition: "<<(double)runTime/CLOCKS_PER_SEC<<"s"<<endl;
//    cout<<"------> Run time for Link addition: "<<(double)runTime/100<<"s"<<endl;
	return true;
}

/**************************************************************************************/
/* Adjust potential in the channel													  */
/**************************************************************************************/
double Qchannel(const double VV,
				const double eepsilon, const int MMaxStep,
				CMatrix3D& pphi_cha, CMatrix3D& pphi_amb, CMatrix3D& UUn,
				const Point& IInitiationPoint, ListLink& EEstablishedLinks,
				StepsSizes dd, const BoxSteps& NN)
{
	bool phiStorage = false;
	ListLink::iterator it;
	
	if(phiStorage==false)
	{
		Potential	P1;//(pphi_cha,UUn);
		SorSolution	SSOR;//(pphi_cha, eepsilon,MMaxStep, dd, NN, P1, UUn);

		pphi_cha[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k]	= VV - pphi_amb[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k];
		for(it=EEstablishedLinks.begin() ; it!= EEstablishedLinks.end() ; it++)
			pphi_cha[it->end.i][it->end.j][it->end.k]	= VV - it->deltaV - pphi_amb[it->end.i][it->end.j][it->end.k];
		
		P1.init(pphi_cha,UUn);
		SSOR.init(pphi_cha, eepsilon,MMaxStep, dd, NN, P1, UUn);
		SSOR.Solve(dd,NN,UUn,pphi_cha);
		
		double QQ = ChannelCharge(Globalrho(pphi_cha,dd,NN),UUn,dd,NN);
		// double QQ = TotalCharge(Globalrho(pphi_cha,dd,NN),dd,NN);
		// cout<<"QQ tot = "<<QQ<<endl;
		return QQ;
	}
	
	else if(phiStorage==true)
	{
		CMatrix1D	pphi_Ca(NN.z);
		CMatrix1D	pphi_Cb(NN.z);
		CMatrix2D	pphi_C2D(NN.y,NN.z);
		CMatrix2D	pphi_A2D(NN.y,NN.z);
		CMatrix1D	pphi_A(NN.z);
		Potential	P1;//(pphi_cha,UUn);
		SorSolution	SSOR;//(pphi_cha, eepsilon,MMaxStep, dd, NN, P1, UUn);
		
		pphi_cha[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k]	= VV - pphi_amb[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k];
		for(it=EEstablishedLinks.begin() ; it!= EEstablishedLinks.end() ; it++)
			pphi_cha[it->end.i][it->end.j][it->end.k]	= VV - it->deltaV - pphi_amb[it->end.i][it->end.j][it->end.k];
		
/*		
		for(int ii= 1 ; ii<NN.x-1 ; ii++) for(int jj=1 ; jj<NN.y-1 ; jj++) for(int kk=1 ; kk<NN.z-1 ; kk++)
		{	
			if(UUn[ii][jj][kk]!=0) 
			{pphi_cha[ii][jj][kk]		= VV - pphi_amb[ii][jj][kk];}
		};
*/		
		for(int kk=0; kk<NN.z ; kk++){
			pphi_Ca[kk]	= pphi_cha((NN.x-1)/2,(NN.y-1)/2,kk);
			pphi_A[kk]	= pphi_amb((NN.x-1)/2,(NN.y-1)/2,kk);
		}
		
		P1.init(pphi_cha,UUn);
		SSOR.init(pphi_cha, eepsilon,MMaxStep, dd, NN, P1, UUn);
		SSOR.Solve(dd,NN,UUn,pphi_cha);
		
		for(int kk=0; kk<NN.z ; kk++){
			pphi_Cb[kk]	= pphi_cha((NN.x-1)/2,(NN.y-1)/2,kk);
			for(int jj=0 ; jj<NN.y-1 ; jj++)
			{
				pphi_C2D[jj][kk] = pphi_cha((NN.x-1)/2,jj,kk);
				pphi_A2D[jj][kk] = pphi_amb((NN.x-1)/2,jj,kk);
			}
		}	
		pphi_Ca.fwrite("results/phiC1Da.dat");
		pphi_Cb.fwrite("results/phiC1Db.dat");
		pphi_A.fwrite("results/phiA1D.dat");
		pphi_C2D.fwrite("results/phiC2D.dat");
		pphi_A2D.fwrite("results/phiA2D.dat");
		
		double QQ = ChannelCharge(Globalrho(pphi_cha,dd,NN),UUn,dd,NN);
//		double QQ = TotalCharge(Globalrho(pphi_cha,dd,NN), dd,NN);
		return QQ;
	}
	return 0;
}
/**************************************************************************************/

/**************************************************************************************/
/* Derive Channel Potential to minimize total Charge - Dichotomy vs. Nelder-Mead	  */
/**************************************************************************************/
double fMinSearch(const double VV, const double QQchannelPlus,
				  const double VVmin, const double VVmax, 
				  const double eepsilon, const int MMaxStep,
				  CMatrix3D& pphi_cha, CMatrix3D& pphi_amb, CMatrix3D& UUn,
				  const Point& IInitiationPoint, ListLink& EEstablishedLinks,
				  StepsSizes dd, const BoxSteps& NN)
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
			SwitchValues(CCl,CCr);
			SwitchValues(QQl,QQr);
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
		
		cout<<"------> Run time for Qminimization: "<<(double)runTime/CLOCKS_PER_SEC<<"s\n";
//        cout<<"------> Run time for Qminimization: "<<(double)runTime/100<<"s\n";
//		cout<<"Vc = "<<setw(12)<<CC<<" ; Qc/Qc+ = "<<setw(12)<<QQav<<" ; Nb of iterations = "<<kk<<endl;
		return CC;
	}
	else // Nelder-Mead //
	{
		clock_t startTime = clock();
		double	CC;
		double	x1(VVmin)	, f1(0.);
		double	x2(VVmax)	, f2(0.);
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
				SwitchValues(x1,x2);
				SwitchValues(f1,f2);
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
		
		cout<<"------> Run time for Qminimization: "<<(double)runTime/CLOCKS_PER_SEC<<"s\n";
//        cout<<"------> Run time for Qminimization: "<<(double)runTime/100<<"s\n";
		cout<<"Vc = "<<setw(12)<<x1<<" ; Qc = "<<setw(12)<<f1<<" ; Nb of iterations = "<<kk<<endl;
		return CC=x1;
	}
}
/**************************************************************************************/

/**************************************************************************************/
double EqualizeAtGroundPotential(const double eepsilon, const int MMaxStep,
								 CMatrix3D& pphi_cha, CMatrix3D& pphi_amb, CMatrix3D& UUn,
								 const Point& IInitiationPoint, ListLink& EEstablishedLinks,
								 StepsSizes dd, const BoxSteps& NN)
{
	bool phiStorage = false;
	ListLink::iterator it;
	
	if(phiStorage==false)
	{
		Potential	P1;//(pphi_cha,UUn);
		SorSolution	SSOR;//(pphi_cha, eepsilon,MMaxStep, dd, NN, P1, UUn);
		
		pphi_cha[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k]	= -pphi_amb[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k];
		for(it=EEstablishedLinks.begin() ; it!= EEstablishedLinks.end() ; it++)
			pphi_cha[it->end.i][it->end.j][it->end.k]	= -pphi_amb[it->end.i][it->end.j][it->end.k];	
		P1.init(pphi_cha,UUn);
		SSOR.init(pphi_cha, eepsilon,MMaxStep, dd, NN, P1, UUn);
		SSOR.Solve(dd,NN,UUn,pphi_cha);
		return 0; // Total potential of the channel phi = phi_amb + phi_cha = phi_amb + (-phi_amb) = 0;
	}
	
	else if(phiStorage==true)
	{
		CMatrix1D	pphi_Ca(NN.z);
		CMatrix1D	pphi_Cb(NN.z);
		CMatrix2D	pphi_C2D(NN.y,NN.z);
		CMatrix2D	pphi_A2D(NN.y,NN.z);
		CMatrix1D	pphi_A(NN.z);
		Potential	P1;//(pphi_cha,UUn);
		SorSolution	SSOR;//(pphi_cha, eepsilon,MMaxStep, dd, NN, P1, UUn);
		
		pphi_cha[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k]	= -pphi_amb[IInitiationPoint.i][IInitiationPoint.j][IInitiationPoint.k];
		for(it=EEstablishedLinks.begin() ; it!= EEstablishedLinks.end() ; it++)
			pphi_cha[it->end.i][it->end.j][it->end.k]	= -pphi_amb[it->end.i][it->end.j][it->end.k];
		
		/*		
			for(int ii= 1 ; ii<NN.x-1 ; ii++) for(int jj=1 ; jj<NN.y-1 ; jj++) for(int kk=1 ; kk<NN.z-1 ; kk++)
		{	
				if(UUn[ii][jj][kk]!=0) 
				{pphi_cha[ii][jj][kk]		= VV - pphi_amb[ii][jj][kk];}
		};
		*/		
		for(int kk=0; kk<NN.z ; kk++){
			pphi_Ca[kk]	= pphi_cha((NN.x-1)/2,(NN.y-1)/2,kk);
			pphi_A[kk]	= pphi_amb((NN.x-1)/2,(NN.y-1)/2,kk);
		}
		
		P1.init(pphi_cha,UUn);
		SSOR.init(pphi_cha, eepsilon,MMaxStep, dd, NN, P1, UUn);
		SSOR.Solve(dd,NN,UUn,pphi_cha);
		
		for(int kk=0; kk<NN.z ; kk++){
			pphi_Cb[kk]	= pphi_cha((NN.x-1)/2,(NN.y-1)/2,kk);
			for(int jj=0 ; jj<NN.y-1 ; jj++)
			{
				pphi_C2D[jj][kk] = pphi_cha((NN.x-1)/2,jj,kk);
				pphi_A2D[jj][kk] = pphi_amb((NN.x-1)/2,jj,kk);
			}
		}	
		pphi_Ca.fwrite("results/phiC1Da.dat");
		pphi_Cb.fwrite("results/phiC1Db.dat");
		pphi_A.fwrite("results/phiA1D.dat");
		pphi_C2D.fwrite("results/phiC2D.dat");
		pphi_A2D.fwrite("results/phiA2D.dat");
		return 0;
	}
	return 0;	
}
/**************************************************************************************/
