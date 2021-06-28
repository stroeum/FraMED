/* BoundaryConditions.cpp */

#include "BoundaryConditions.h"
/**************************************************************************************/
void BC::Apply(int BBCtype, CMatrix3D& pphi, CMatrix3D rrho, const ResGrid& dd, const SizeGrid& NN)
{
	clock_t startTime = clock();
	/**********************************************************************************/
	/* This function modifies the field at boundaries else than z = 0, however, once  */
	/* fixed, this potential is not modified.										  */
	/* If BBCtype =																	  */
	/*		0. ``Tin Can Model": the potential on the boundaries is fixed at 0		  */
	/*		1. ``Open BC Model": PEC ground, open BC everywhere on the sides		  */
	/*      2. ``Capacitor Model": PEC ground and electrosphere, open side boundaries */
	/*      3. ``Free Space Model": open side, bottom and top boundaries			  */
	/**********************************************************************************/

	if(BBCtype!=0 && BBCtype!=1 && BBCtype!=2 && BBCtype!=3)
	{
		printf("ee: Incorrect type of boundary conditions.\n");
		printf("xx: Program terminating.\n");
		exit(4);
	}
	if(BBCtype == 0) // PEC boundaries //
	{
		// No modifications required //
	}
	else if(BBCtype == 1) // PEC ground, side and upper boundaries open. //
	{
		/******************************************************************************/
		/* To reduce the time of computation, only the source points exceeding xx% of */
		/* the maximal charge density -in absolute value- will be considered.		  */
		/******************************************************************************/
		double rhoMax = 0;
		for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
			if(fabs(rrho[ii][jj][kk]) > rhoMax) rhoMax = fabs(rrho[ii][jj][kk]);

		/******************************************************************************/
		/* We neglect the ambient Laplacian field on Earth (100 V/m = 1e-3kV/cm)	  */
		/* The potential due to the ambient Laplacian field VL = 0.					  */
		/* To take into account this potential, VL = 100 * z (in V)					  */
		/* NB: subscript "s" refers to the source, which volume we are integrating on.*/
		/*	   Hence "s" refers to variable of integration and no subscript refers to */
		/*     the position of the boundary.										  */
		/******************************************************************************/
		double	Eambient = 0;
		double	VL;

		for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
			/**************************************************************************/
			/* Check if point belongs to a boundary									  */
			/**************************************************************************/
			if(
				(kk == NN.z-1													) || // Up    //
				(ii == 0		&& kk > 0 && kk<NN.z-1							) || // Back  //
				(ii == NN.x-1	&& kk > 0 && kk<NN.z-1							) || // Front //
				(jj == 0		&& kk > 0 && kk<NN.z-1  && ii > 0 && ii<NN.x-1	) || // Left  //
				(jj == NN.y-1	&& kk > 0 && kk<NN.z-1  && ii > 0 && ii<NN.x-1	)    // Right //
			  )
			{
				VL = Eambient * kk *dd.z;
				pphi[ii][jj][kk] = VL;
				for(int is=0 ; is<NN.x ; is++) for(int js=0 ; js<NN.y ; js++) for(int ks=0 ; ks<NN.z ; ks++)
					if(!(is == ii && js == jj && ks == kk))
						if(fabs(rrho[is][js][ks]) > .01*rhoMax)
							pphi[ii][jj][kk] +=
								dd.x*dd.y*dd.z/(4*M_PI*eps0) *
								(rrho[is][js][ks]/sqrt(pow( (ii-is)*dd.x,2 ) + pow( (jj-js)*dd.y,2 ) + pow( (kk-ks)*dd.z,2 )) -
								 rrho[is][js][ks]/sqrt(pow( (ii-is)*dd.x,2 ) + pow( (jj-js)*dd.y,2 ) + pow( (kk+ks)*dd.z,2 )) );
			};
	}
	else if(BBCtype == 2) // PEC ground and electrosphere, side boundaries open.//
	{
		/******************************************************************************/
		/* To reduce the time of computation, only the source points exceeding xx% of */
		/* the maximal charge density -in absolute value- will be considered.		  */
		/******************************************************************************/
		double rhoMax = 0;
		for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
			if(fabs(rrho[ii][jj][kk]) > rhoMax) rhoMax = fabs(rrho[ii][jj][kk]);

		/******************************************************************************/
		/* We neglect the ambient Laplacian field on Earth (100 V/m = 1e-3kV/cm)	  */
		/* The potential due to the ambient Laplacian field VL = 0.					  */
		/* To take into account this potential, VL = 100 * z (in V)					  */
		/* NB: subscript "s" refers to the source, which volume we are integrating on.*/
		/*	   Hence "s" refers to variable of integration and no subscript refers to */
		/*     the position of the boundary.										  */
		/******************************************************************************/
		double	Eambient = 0;
		double	VL;

		/******************************************************************************/
		/* Define number of images and store their positions						  */
		/******************************************************************************/
		int	M(25);				// Account for M ground images and M ionospheric images
		int k_GndImg = 0;		// Altitude of ground images
		int k_IonImg = 0;		// Altitude of iono/electrosphere images
		int	k_Ion    = NN.z-1;	// Altitude coordinate of the iono/electrosphere

		/******************************************************************************/
		/* Derive potential at the boundaries										  */
		/******************************************************************************/
		for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
			/**************************************************************************/
			/* Check if point belongs to a side boundary							  */
			/**************************************************************************/
			if(
				(ii == 0		&& kk > 0 && kk<NN.z-1							) || // Back  //
				(ii == NN.x-1	&& kk > 0 && kk<NN.z-1							) || // Front //
				(jj == 0		&& kk > 0 && kk<NN.z-1  && ii > 0 && ii<NN.x-1	) || // Left  //
				(jj == NN.y-1	&& kk > 0 && kk<NN.z-1  && ii > 0 && ii<NN.x-1	)    // Right //
				)
			{
				VL = Eambient * kk *dd.z;
				pphi[ii][jj][kk] = VL;
				for(int is=0 ; is<NN.x ; is++) for(int js=0 ; js<NN.y ; js++) for(int ks=0 ; ks<NN.z ; ks++)
					if(!(is == ii && js == jj && ks == kk))
						if(fabs(rrho[is][js][ks]) > .01*rhoMax)
						{
							k_GndImg = ks; // Altitude of ground images			   //
							k_IonImg = ks; // Altitude of iono/electrosphere images //

							pphi[ii][jj][kk] += dd.x*dd.y*dd.z/(4*M_PI*eps0) *
								( rrho[is][js][ks]/sqrt(pow( (ii-is)*dd.x,2 ) + pow( (jj-js)*dd.y,2 ) + pow( (kk- ks)*dd.z,2 ))); // Charge
							for(int mm=1; mm<=M; mm++)
							{
								k_GndImg = k_GndImg - ( mm%2*2*ks + (mm-1)%2*2*(k_Ion-ks) );
								k_IonImg = k_IonImg + ( (mm-1)%2*2*ks + mm%2*2*(k_Ion-ks) );
								pphi[ii][jj][kk] += dd.x*dd.y*dd.z/(4*M_PI*eps0) * pow(-1.0,mm)*
									(rrho[is][js][ks]/sqrt(pow( (ii-is)*dd.x,2 ) + pow( (jj-js)*dd.y,2 ) + pow( (kk- k_GndImg)*dd.z,2 )) +	// Ground Images
									 rrho[is][js][ks]/sqrt(pow( (ii-is)*dd.x,2 ) + pow( (jj-js)*dd.y,2 ) + pow( (kk- k_IonImg)*dd.z,2 )) ); // Ionosphere Images
							}
						};
			};
	}
	else if(BBCtype == 3) // Free space, side, upper, and lower boundaries open. //
	{
		/******************************************************************************/
		/* To reduce the time of computation, only the source points exceeding xx% of */
		/* the maximal charge density -in absolute value- will be considered.		  */
		/******************************************************************************/
		double rhoMax = 0;
		for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
			if(fabs(rrho[ii][jj][kk]) > rhoMax) rhoMax = fabs(rrho[ii][jj][kk]);

		/******************************************************************************/
		/* We neglect the ambient Laplacian field on Earth (100 V/m = 1e-3kV/cm)	  */
		/* The potential due to the ambient Laplacian field VL = 0.					  */
		/* To take into account this potential, VL = 100 * z (in V)					  */
		/* NB: subscript "s" refers to the source, which volume we are integrating on.*/
		/*	   Hence "s" refers to variable of integration and no subscript refers to */
		/*     the position of the boundary.										  */
		/******************************************************************************/
		double	Eambient = 0;
		double	VL;

		for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
			/**************************************************************************/
			/* Check if point belongs to a boundary									  */
			/**************************************************************************/
			if(
				(kk == 0														) || // Down  //
				(kk == NN.z-1													) || // Up    //
				(ii == 0		&& kk > 0 && kk<NN.z-1							) || // Back  //
				(ii == NN.x-1	&& kk > 0 && kk<NN.z-1							) || // Front //
				(jj == 0		&& kk > 0 && kk<NN.z-1  && ii > 0 && ii<NN.x-1	) || // Left  //
				(jj == NN.y-1	&& kk > 0 && kk<NN.z-1  && ii > 0 && ii<NN.x-1	)    // Right //
			  )
			{
				VL = Eambient * kk *dd.z;
				pphi[ii][jj][kk] = VL;
				for(int is=0 ; is<NN.x ; is++) for(int js=0 ; js<NN.y ; js++) for(int ks=0 ; ks<NN.z ; ks++)
					if(!(is == ii && js == jj && ks == kk))
						if(fabs(rrho[is][js][ks]) > .01*rhoMax)
							pphi[ii][jj][kk] +=
								dd.x*dd.y*dd.z/(4*M_PI*eps0) *
								rrho[is][js][ks]/sqrt(pow( (ii-is)*dd.x,2 ) + pow( (jj-js)*dd.y,2 ) + pow( (kk-ks)*dd.z,2 ));
			};
	}
	clock_t endTime = clock();
	clock_t runTime;
	runTime = endTime - startTime;
//	printf("Run time for application of open BC        : %fs\n",(double)runTime/100);
}

/**************************************************************************************/
void BC::Update(bool iisFlashAccoutedForInBC, int BBCtype, CMatrix3D& pphi_cha, const double rrhoAmbMin, const double rrhoAmbMax, const ResGrid& dd, const SizeGrid& NN)
{
	clock_t startTime	= clock();
	clock_t endTime		= clock();
	clock_t runTime		= endTime - startTime;	;

	if(iisFlashAccoutedForInBC==true)
	{
		/******************************************************************************/
		/* This function modifies the field at boundaries else than z = 0, however,   */
		/* once fixed, this potential is not modified.								  */
		/* If BBCtype =																  */
		/*		0. ``Tin Can Model": the potential on the boundaries is fixed at 0	  */
		/*		1. ``Open BC Model": PEC ground, open BC everywhere on the sides	  */
		/*      2. ``Capacitor Model": PEC ground and electrosphere, open side		  */
		/*         boundaries														  */
		/*      3. ``Free Space Model": open side, bottom and top boundaries		  */
		/******************************************************************************/

		if(BBCtype == 0) // PEC boundaries //
		{
			cout<<"ii: Cannot proceed to inclusion of all images in the Tin Can Model yet..."<<endl;
			cout<<"xx: Program terminating.\n"<<endl;
			exit(5);
		}
		else if(BBCtype == 1) // PEC ground, side and upper boundaries open. //
		{
			CMatrix3D	rrho_cha(NN.x,NN.y,NN.z);
			double		AAbsRhoAmbMax = max(fabs(rrhoAmbMin),fabs(rrhoAmbMax));
			double		AAbsRhoChaMax = 0;
			for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
			{
				rrho_cha[ii][jj][kk] = foo::rhoijk(ii,jj,kk,pphi_cha,dd,NN);
				if(fabs(rrho_cha[ii][jj][kk])>= AAbsRhoChaMax)
					AAbsRhoChaMax = fabs(rrho_cha[ii][jj][kk]);
			};
			if(AAbsRhoChaMax>.01*AAbsRhoAmbMax)
			{
				for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
					/**********************************************************************/
					/* Check if point belongs to a boundary								  */
					/**********************************************************************/
					if(
						(kk == NN.z-1													) || // Up    //
						(ii == 0		&& kk > 0 && kk<NN.z-1							) || // Back  //
						(ii == NN.x-1	&& kk > 0 && kk<NN.z-1							) || // Front //
						(jj == 0		&& kk > 0 && kk<NN.z-1  && ii > 0 && ii<NN.x-1	) || // Left  //
						(jj == NN.y-1	&& kk > 0 && kk<NN.z-1  && ii > 0 && ii<NN.x-1	)    // Right //
						)
					{
						pphi_cha[ii][jj][kk] = 0;
						for(int is=0 ; is<NN.x ; is++) for(int js=0 ; js<NN.y ; js++) for(int ks=0 ; ks<NN.z ; ks++)
							if(!(is == ii && js == jj && ks == kk))
								if(fabs(rrho_cha[is][js][ks]) > .01*AAbsRhoAmbMax)
									pphi_cha[ii][jj][kk] +=
										dd.x*dd.y*dd.z/(4*M_PI*eps0) *
										(rrho_cha[is][js][ks]/sqrt(pow( (ii-is)*dd.x,2 ) + pow( (jj-js)*dd.y,2 ) + pow( (kk-ks)*dd.z,2 )) -
										 rrho_cha[is][js][ks]/sqrt(pow( (ii-is)*dd.x,2 ) + pow( (jj-js)*dd.y,2 ) + pow( (kk+ks)*dd.z,2 )) );
					};
				endTime = clock();
				runTime = endTime - startTime;
				cout<<"ii:\t Run time for BC update: "<<(double)runTime/CLOCKS_PER_SEC<<" s."<<endl;
			}
			else {cout<<"ii:\t No need for BC update.\n";};
		}
		else if(BBCtype == 2) // PEC ground and electrosphere, side boundaries open.//
		{
			/******************************************************************************/
			/* Define variables for the charge density									  */
			/******************************************************************************/
			CMatrix3D	rrho_cha(NN.x,NN.y,NN.z);
			double		AAbsRhoAmbMax = max(fabs(rrhoAmbMin),fabs(rrhoAmbMax));
			double		AAbsRhoChaMax = 0;

			/******************************************************************************/
			/* Define number of images and store their positions						  */
			/******************************************************************************/
			int	M(25);				// Account for M ground images and M ionospheric images
			int k_GndImg = 0;		// Altitude of ground images
			int k_IonImg = 0;		// Altitude of iono/electrosphere images
			int	k_Ion    = NN.z-1;	// Altitude coordinate of the iono/electrosphere

			for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
			{
				rrho_cha[ii][jj][kk] = foo::rhoijk(ii,jj,kk,pphi_cha,dd,NN);
				if(fabs(rrho_cha[ii][jj][kk])>= AAbsRhoChaMax)
					AAbsRhoChaMax = fabs(rrho_cha[ii][jj][kk]);
			};
			if(AAbsRhoChaMax>.01*AAbsRhoAmbMax)
			{
				for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
					/******************************************************************/
					/* Check if point belongs to a boundary							  */
					/******************************************************************/
					if(
						(ii == 0		&& kk > 0 && kk<NN.z-1							) || // Back  //
						(ii == NN.x-1	&& kk > 0 && kk<NN.z-1							) || // Front //
						(jj == 0		&& kk > 0 && kk<NN.z-1  && ii > 0 && ii<NN.x-1	) || // Left  //
						(jj == NN.y-1	&& kk > 0 && kk<NN.z-1  && ii > 0 && ii<NN.x-1	)    // Right //
						)
					{
						pphi_cha[ii][jj][kk] = 0;
						for(int is=0 ; is<NN.x ; is++) for(int js=0 ; js<NN.y ; js++) for(int ks=0 ; ks<NN.z ; ks++)
						{
							if(!(is == ii && js == jj && ks == kk))
								if(fabs(rrho_cha[is][js][ks]) > .01*AAbsRhoAmbMax)
								{
									k_GndImg = ks; // Altitude of ground images			   //
									k_IonImg = ks; // Altitude of iono/electrosphere images //

									pphi_cha[ii][jj][kk] += dd.x*dd.y*dd.z/(4*M_PI*eps0) *
										( rrho_cha[is][js][ks]/sqrt(pow( (ii-is)*dd.x,2 ) + pow( (jj-js)*dd.y,2 ) + pow( (kk- ks)*dd.z,2 ))); // Charge
									for(int mm=1; mm<=M; mm++)
									{
										k_GndImg = k_GndImg - ( mm%2*2*ks + (mm-1)%2*2*(k_Ion-ks) );
										k_IonImg = k_IonImg + ( (mm-1)%2*2*ks + mm%2*2*(k_Ion-ks) );
										pphi_cha[ii][jj][kk] += dd.x*dd.y*dd.z/(4*M_PI*eps0) * pow(-1.0,mm)*
											(rrho_cha[is][js][ks]/sqrt(pow( (ii-is)*dd.x,2 ) + pow( (jj-js)*dd.y,2 ) + pow( (kk- k_GndImg)*dd.z,2 )) +	// Ground Images
											 rrho_cha[is][js][ks]/sqrt(pow( (ii-is)*dd.x,2 ) + pow( (jj-js)*dd.y,2 ) + pow( (kk- k_IonImg)*dd.z,2 )) ); // Ionosphere Images
									};
								};
						};
					};
				endTime = clock();
				runTime = endTime - startTime;
				cout<<"ii:\t Run time for BC update: "<<(double)runTime/CLOCKS_PER_SEC<<" s."<<endl;
			}
			else cout<<"ii:\t No need for BC update.\n";
		}
		else if(BBCtype == 3) // Free Space, side, lower, and upper boundaries open. //
		{
			CMatrix3D	rrho_cha(NN.x,NN.y,NN.z);
			double		AAbsRhoAmbMax = max(fabs(rrhoAmbMin),fabs(rrhoAmbMax));
			double		AAbsRhoChaMax = 0;
			for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
			{
				rrho_cha[ii][jj][kk] = foo::rhoijk(ii,jj,kk,pphi_cha,dd,NN);
				if(fabs(rrho_cha[ii][jj][kk])>= AAbsRhoChaMax)
					AAbsRhoChaMax = fabs(rrho_cha[ii][jj][kk]);
			};
			if(AAbsRhoChaMax>.01*AAbsRhoAmbMax)
			{
				for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
					/**********************************************************************/
					/* Check if point belongs to a boundary								  */
					/**********************************************************************/
					if(
						(kk == 0														) || // Down  //
						(kk == NN.z-1													) || // Up    //
						(ii == 0		&& kk > 0 && kk<NN.z-1							) || // Back  //
						(ii == NN.x-1	&& kk > 0 && kk<NN.z-1							) || // Front //
						(jj == 0		&& kk > 0 && kk<NN.z-1  && ii > 0 && ii<NN.x-1	) || // Left  //
						(jj == NN.y-1	&& kk > 0 && kk<NN.z-1  && ii > 0 && ii<NN.x-1	)    // Right //
						)
					{
						pphi_cha[ii][jj][kk] = 0;
						for(int is=0 ; is<NN.x ; is++) for(int js=0 ; js<NN.y ; js++) for(int ks=0 ; ks<NN.z ; ks++)
							if(!(is == ii && js == jj && ks == kk))
								if(fabs(rrho_cha[is][js][ks]) > .01*AAbsRhoAmbMax)
									pphi_cha[ii][jj][kk] +=
										dd.x*dd.y*dd.z/(4*M_PI*eps0) *
										rrho_cha[is][js][ks]/sqrt(pow( (ii-is)*dd.x,2 ) + pow( (jj-js)*dd.y,2 ) + pow( (kk-ks)*dd.z,2 ));
					};
				endTime = clock();
				runTime = endTime - startTime;
				cout<<"ii:\t Run time for BC update: "<<(double)runTime/CLOCKS_PER_SEC<<" s."<<endl;
			}
			else {cout<<"ii:\t No need for BC update.\n";};
		};
	}
	else cout<<"ii:\t No update done."<<endl;
}



/**************************************************************************************/
void BC::AddUniformE(bool uuniformE, const double& EEo, CMatrix3D& pphi, CMatrix3D& UUn, const SizeDomain& LL, const ResGrid& dd, const SizeGrid& NN)
{
	clock_t startTime = clock();
	if(uuniformE == true)
	{
		cout<<"Uniform E-field Added\n";
		/******************************************************************************/
		/* Upper and Lower boundaries												  */
		/******************************************************************************/
		for(int ii=0 ; ii<NN.x ; ii++)
			for(int jj=0 ; jj<NN.y ; jj++)
			{
				UUn(ii,jj,0)		= 1;
				pphi(ii,jj,0)		+= 0;
				UUn(ii,jj,NN.z-1)	= 1;
				pphi(ii,jj,NN.z-1)	+= -EEo*(dd.z*(NN.z-1));
			}

		/******************************************************************************/
		/* Front and Back boundaries												  */
		/******************************************************************************/
		for(int jj=0 ; jj<NN.y ; jj++)
			for(int kk=1 ; kk<NN.z-1 ; kk++)
			{
				UUn(0,jj,kk)		= 1;
				pphi(0,jj,kk)		+= -EEo*(kk*dd.z);;
				UUn(NN.x-1,jj,kk)	= 1;
				pphi(NN.x-1,jj,kk)	+= -EEo*(kk*dd.z);
			}

		/******************************************************************************/
		/* Left and Right boundaries												  */
		/******************************************************************************/
		for(int ii=1 ; ii<NN.x-1 ; ii++)
			for(int kk=1 ; kk<NN.z-1 ; kk++)
			{
				UUn(ii,0,kk)		= 1;
				pphi(ii,0,kk)		+= -EEo*(kk*dd.z);;
				UUn(ii,NN.y-1,kk)	= 1;
				pphi(ii,NN.y-1,kk)	+= -EEo*(kk*dd.z);
			}
	}
	else cout<<"No Uniform E-field Applied.\n";

	clock_t endTime = clock();
	clock_t runTime = endTime - startTime;
	printf("Run time for application of BC: %lf s.\n",(double)runTime/CLOCKS_PER_SEC);
}
/**************************************************************************************/
