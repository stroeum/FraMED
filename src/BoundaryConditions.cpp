/* BoundaryConditions.cpp */

#include "BoundaryConditions.h"
/**************************************************************************************/
void BC::Apply(int _BCtype, CMatrix3D& _phi, CMatrix3D _rho, const ResGrid& _d, const SizeGrid& _N)
{
	clock_t _StartTime = clock();
	/**********************************************************************************/
	/* This function modifies the field at boundaries else than z = 0, however, once  */
	/* fixed, this potential is not modified.										  */
	/* If _BCtype =																	  */
	/*		0. ``Tin Can Model": the potential on the boundaries is fixed at 0		  */
	/*		1. ``Open BC Model": PEC ground, open BC everywhere on the sides		  */
	/*      2. ``Capacitor Model": PEC ground and electrosphere, open side boundaries */
	/*      3. ``Free Space Model": open side, bottom and top boundaries			  */
	/**********************************************************************************/

	if(_BCtype!=TIN_CAN && _BCtype!=OPEN_BC && _BCtype!=G_G && _BCtype!=FREE_SPACE)
	{
		printf("ee: Incorrect type of boundary conditions.\n");
		printf("xx: Program terminating.\n");
		exit(4);
	}
	if(_BCtype == TIN_CAN) // PEC boundaries //
	{
		// No modifications required //
	}
	else if(_BCtype == OPEN_BC) // PEC ground, side and upper boundaries open. //
	{
		/******************************************************************************/
		/* To reduce the time of computation, only the source points exceeding xx% of */
		/* the maximal charge density -in absolute value- will be considered.		  */
		/******************************************************************************/
		double _RhoMax = 0;
		for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
			if(fabs(_rho[ii][jj][kk]) > _RhoMax) _RhoMax = fabs(_rho[ii][jj][kk]);

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

		for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
			/**************************************************************************/
			/* Check if point belongs to a boundary									  */
			/**************************************************************************/
			if(
				(kk == _N.z-1													) || // Up    //
				(ii == 0		&& kk > 0 && kk<_N.z-1							) || // Back  //
				(ii == _N.x-1	&& kk > 0 && kk<_N.z-1							) || // Front //
				(jj == 0		&& kk > 0 && kk<_N.z-1  && ii > 0 && ii<_N.x-1	) || // Left  //
				(jj == _N.y-1	&& kk > 0 && kk<_N.z-1  && ii > 0 && ii<_N.x-1	)    // Right //
			  )
			{
				VL = Eambient * kk *_d.z;
				_phi[ii][jj][kk] = VL;
				for(int is=0 ; is<_N.x ; is++) for(int js=0 ; js<_N.y ; js++) for(int ks=0 ; ks<_N.z ; ks++)
					if(!(is == ii && js == jj && ks == kk))
						if(fabs(_rho[is][js][ks]) > .01*_RhoMax)
							_phi[ii][jj][kk] +=
								_d.x*_d.y*_d.z/(4*PMC.pi*PMC.eps0) *
								(_rho[is][js][ks]/sqrt(pow( (ii-is)*_d.x,2 ) + pow( (jj-js)*_d.y,2 ) + pow( (kk-ks)*_d.z,2 )) -
								 _rho[is][js][ks]/sqrt(pow( (ii-is)*_d.x,2 ) + pow( (jj-js)*_d.y,2 ) + pow( (kk+ks)*_d.z,2 )) );
			};
	}
	else if(_BCtype == G_G) // PEC ground and electrosphere, side boundaries open.//
	{
		/******************************************************************************/
		/* To reduce the time of computation, only the source points exceeding xx% of */
		/* the maximal charge density -in absolute value- will be considered.		  */
		/******************************************************************************/
		double _RhoMax = 0;
		for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
			if(fabs(_rho[ii][jj][kk]) > _RhoMax) _RhoMax = fabs(_rho[ii][jj][kk]);

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
		int	k_Ion    = _N.z-1;	// Altitude coordinate of the iono/electrosphere

		/******************************************************************************/
		/* Derive potential at the boundaries										  */
		/******************************************************************************/
		for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
			/**************************************************************************/
			/* Check if point belongs to a side boundary							  */
			/**************************************************************************/
			if(
				(ii == 0		&& kk > 0 && kk<_N.z-1							) || // Back  //
				(ii == _N.x-1	&& kk > 0 && kk<_N.z-1							) || // Front //
				(jj == 0		&& kk > 0 && kk<_N.z-1  && ii > 0 && ii<_N.x-1	) || // Left  //
				(jj == _N.y-1	&& kk > 0 && kk<_N.z-1  && ii > 0 && ii<_N.x-1	)    // Right //
				)
			{
				VL = Eambient * kk *_d.z;
				_phi[ii][jj][kk] = VL;
				for(int is=0 ; is<_N.x ; is++) for(int js=0 ; js<_N.y ; js++) for(int ks=0 ; ks<_N.z ; ks++)
					if(!(is == ii && js == jj && ks == kk))
						if(fabs(_rho[is][js][ks]) > .01*_RhoMax)
						{
							k_GndImg = ks; // Altitude of ground images			   //
							k_IonImg = ks; // Altitude of iono/electrosphere images //

							_phi[ii][jj][kk] += _d.x*_d.y*_d.z/(4*PMC.pi*PMC.eps0) *
								( _rho[is][js][ks]/sqrt(pow( (ii-is)*_d.x,2 ) + pow( (jj-js)*_d.y,2 ) + pow( (kk- ks)*_d.z,2 ))); // Charge
							for(int mm=1; mm<=M; mm++)
							{
								k_GndImg = k_GndImg - ( mm%2*2*ks + (mm-1)%2*2*(k_Ion-ks) );
								k_IonImg = k_IonImg + ( (mm-1)%2*2*ks + mm%2*2*(k_Ion-ks) );
								_phi[ii][jj][kk] += _d.x*_d.y*_d.z/(4*PMC.pi*PMC.eps0) * pow(-1.0,mm)*
									(_rho[is][js][ks]/sqrt(pow( (ii-is)*_d.x,2 ) + pow( (jj-js)*_d.y,2 ) + pow( (kk- k_GndImg)*_d.z,2 )) +	// Ground Images
									 _rho[is][js][ks]/sqrt(pow( (ii-is)*_d.x,2 ) + pow( (jj-js)*_d.y,2 ) + pow( (kk- k_IonImg)*_d.z,2 )) ); // Ionosphere Images
							}
						};
			};
	}
	else if(_BCtype == FREE_SPACE) // Free space, side, upper, and lower boundaries open. //
	{
		/******************************************************************************/
		/* To reduce the time of computation, only the source points exceeding xx% of */
		/* the maximal charge density -in absolute value- will be considered.		  */
		/******************************************************************************/
		double _RhoMax = 0;
		for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
			if(fabs(_rho[ii][jj][kk]) > _RhoMax) _RhoMax = fabs(_rho[ii][jj][kk]);

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

		for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
			/**************************************************************************/
			/* Check if point belongs to a boundary									  */
			/**************************************************************************/
			if(
				(kk == 0														) || // Down  //
				(kk == _N.z-1													) || // Up    //
				(ii == 0		&& kk > 0 && kk<_N.z-1							) || // Back  //
				(ii == _N.x-1	&& kk > 0 && kk<_N.z-1							) || // Front //
				(jj == 0		&& kk > 0 && kk<_N.z-1  && ii > 0 && ii<_N.x-1	) || // Left  //
				(jj == _N.y-1	&& kk > 0 && kk<_N.z-1  && ii > 0 && ii<_N.x-1	)    // Right //
			  )
			{
				VL = Eambient * kk *_d.z;
				_phi[ii][jj][kk] = VL;
				for(int is=0 ; is<_N.x ; is++) for(int js=0 ; js<_N.y ; js++) for(int ks=0 ; ks<_N.z ; ks++)
					if(!(is == ii && js == jj && ks == kk))
						if(fabs(_rho[is][js][ks]) > .01*_RhoMax)
							_phi[ii][jj][kk] +=
								_d.x*_d.y*_d.z/(4*PMC.pi*PMC.eps0) *
								_rho[is][js][ks]/sqrt(pow( (ii-is)*_d.x,2 ) + pow( (jj-js)*_d.y,2 ) + pow( (kk-ks)*_d.z,2 ));
			};
	}
	clock_t _EndTime = clock();
	clock_t _RunTime;
	_RunTime = _EndTime - _StartTime;
//	printf("Run time for application of open BC        : %fs\n",(double)_RunTime/100);
}

/**************************************************************************************/
void BC::Update(bool iisFlashAccoutedForInBC, int _BCtype, CMatrix3D& pphi_cha, const double rrhoAmbMin, const double rrhoAmbMax, const ResGrid& _d, const SizeGrid& _N)
{
	clock_t _StartTime	= clock();
	clock_t _EndTime	= clock();
	clock_t _RunTime	= _EndTime - _StartTime;	;
    string msg;
	if(iisFlashAccoutedForInBC==true)
	{
		/******************************************************************************/
		/* This function modifies the field at boundaries else than z = 0, however,   */
		/* once fixed, this potential is not modified.								  */
		/* If _BCtype =																  */
		/*		0. ``Tin Can Model": the potential on the boundaries is fixed at 0	  */
		/*		1. ``Open BC Model": PEC ground, open BC everywhere on the sides	  */
		/*      2. ``Capacitor Model": PEC ground and electrosphere, open side		  */
		/*         boundaries														  */
		/*      3. ``Free Space Model": open side, bottom and top boundaries		  */
		/******************************************************************************/

		if(_BCtype == TIN_CAN) // PEC boundaries //
		{
			printf("ii: Cannot proceed to inclusion of all images in the Tin Can Model yet...\n");
			printf("xx: Program terminating.\n");
			exit(5);
		}
		else if(_BCtype == OPEN_BC) // PEC ground, side and upper boundaries open. //
		{
			CMatrix3D	rrho_cha(_N.x,_N.y,_N.z);
			double		AAbsRhoAmbMax = max(fabs(rrhoAmbMin),fabs(rrhoAmbMax));
			double		AAbsRhoChaMax = 0;
			for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
			{
				rrho_cha[ii][jj][kk] = foo::rhoijk(ii,jj,kk,pphi_cha,_d,_N);
				if(fabs(rrho_cha[ii][jj][kk])>= AAbsRhoChaMax)
					AAbsRhoChaMax = fabs(rrho_cha[ii][jj][kk]);
			};
			if(AAbsRhoChaMax>.01*AAbsRhoAmbMax)
			{
				for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
					/**********************************************************************/
					/* Check if point belongs to a boundary								  */
					/**********************************************************************/
					if(
						(kk == _N.z-1													) || // Up    //
						(ii == 0		&& kk > 0 && kk<_N.z-1							) || // Back  //
						(ii == _N.x-1	&& kk > 0 && kk<_N.z-1							) || // Front //
						(jj == 0		&& kk > 0 && kk<_N.z-1  && ii > 0 && ii<_N.x-1	) || // Left  //
						(jj == _N.y-1	&& kk > 0 && kk<_N.z-1  && ii > 0 && ii<_N.x-1	)    // Right //
						)
					{
						pphi_cha[ii][jj][kk] = 0;
						for(int is=0 ; is<_N.x ; is++) for(int js=0 ; js<_N.y ; js++) for(int ks=0 ; ks<_N.z ; ks++)
							if(!(is == ii && js == jj && ks == kk))
								if(fabs(rrho_cha[is][js][ks]) > .01*AAbsRhoAmbMax)
									pphi_cha[ii][jj][kk] +=
										_d.x*_d.y*_d.z/(4*PMC.pi*PMC.eps0) *
										(rrho_cha[is][js][ks]/sqrt(pow( (ii-is)*_d.x,2 ) + pow( (jj-js)*_d.y,2 ) + pow( (kk-ks)*_d.z,2 )) -
										 rrho_cha[is][js][ks]/sqrt(pow( (ii-is)*_d.x,2 ) + pow( (jj-js)*_d.y,2 ) + pow( (kk+ks)*_d.z,2 )) );
					};
				_EndTime = clock();
				_RunTime = _EndTime - _StartTime;
                msg = "ii:\t Run time for BC update: " + to_string((double)_RunTime/CLOCKS_PER_SEC)  + " s.\n";
				printf("%s",msg.c_str());
			}
			else
                printf("ii:\t No need for BC update.\n");
		}
		else if(_BCtype == G_G) // PEC ground and electrosphere, side boundaries open.//
		{
			/******************************************************************************/
			/* Define variables for the charge density									  */
			/******************************************************************************/
			CMatrix3D	rrho_cha(_N.x,_N.y,_N.z);
			double		AAbsRhoAmbMax = max(fabs(rrhoAmbMin),fabs(rrhoAmbMax));
			double		AAbsRhoChaMax = 0;

			/******************************************************************************/
			/* Define number of images and store their positions						  */
			/******************************************************************************/
			int	M(25);				// Account for M ground images and M ionospheric images
			int k_GndImg = 0;		// Altitude of ground images
			int k_IonImg = 0;		// Altitude of iono/electrosphere images
			int	k_Ion    = _N.z-1;	// Altitude coordinate of the iono/electrosphere

			for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
			{
				rrho_cha[ii][jj][kk] = foo::rhoijk(ii,jj,kk,pphi_cha,_d,_N);
				if(fabs(rrho_cha[ii][jj][kk])>= AAbsRhoChaMax)
					AAbsRhoChaMax = fabs(rrho_cha[ii][jj][kk]);
			};
			if(AAbsRhoChaMax>.01*AAbsRhoAmbMax)
			{
				for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
					/******************************************************************/
					/* Check if point belongs to a boundary							  */
					/******************************************************************/
					if(
						(ii == 0		&& kk > 0 && kk<_N.z-1							) || // Back  //
						(ii == _N.x-1	&& kk > 0 && kk<_N.z-1							) || // Front //
						(jj == 0		&& kk > 0 && kk<_N.z-1  && ii > 0 && ii<_N.x-1	) || // Left  //
						(jj == _N.y-1	&& kk > 0 && kk<_N.z-1  && ii > 0 && ii<_N.x-1	)    // Right //
						)
					{
						pphi_cha[ii][jj][kk] = 0;
						for(int is=0 ; is<_N.x ; is++) for(int js=0 ; js<_N.y ; js++) for(int ks=0 ; ks<_N.z ; ks++)
						{
							if(!(is == ii && js == jj && ks == kk))
								if(fabs(rrho_cha[is][js][ks]) > .01*AAbsRhoAmbMax)
								{
									k_GndImg = ks; // Altitude of ground images			   //
									k_IonImg = ks; // Altitude of iono/electrosphere images //

									pphi_cha[ii][jj][kk] += _d.x*_d.y*_d.z/(4*PMC.pi*PMC.eps0) *
										( rrho_cha[is][js][ks]/sqrt(pow( (ii-is)*_d.x,2 ) + pow( (jj-js)*_d.y,2 ) + pow( (kk- ks)*_d.z,2 ))); // Charge
									for(int mm=1; mm<=M; mm++)
									{
										k_GndImg = k_GndImg - ( mm%2*2*ks + (mm-1)%2*2*(k_Ion-ks) );
										k_IonImg = k_IonImg + ( (mm-1)%2*2*ks + mm%2*2*(k_Ion-ks) );
										pphi_cha[ii][jj][kk] += _d.x*_d.y*_d.z/(4*PMC.pi*PMC.eps0) * pow(-1.0,mm)*
											(rrho_cha[is][js][ks]/sqrt(pow( (ii-is)*_d.x,2 ) + pow( (jj-js)*_d.y,2 ) + pow( (kk- k_GndImg)*_d.z,2 )) +	// Ground Images
											 rrho_cha[is][js][ks]/sqrt(pow( (ii-is)*_d.x,2 ) + pow( (jj-js)*_d.y,2 ) + pow( (kk- k_IonImg)*_d.z,2 )) ); // Ionosphere Images
									};
								};
						};
					};
				_EndTime = clock();
				_RunTime = _EndTime - _StartTime;
                msg = "ii:\t Run time for BC update: " + to_string((double)_RunTime/CLOCKS_PER_SEC)  + " s.\n";
                printf("%s",msg.c_str());
			}
			else printf("ii:\t No need for BC update.\n");
		}
		else if(_BCtype == FREE_SPACE) // Free Space, side, lower, and upper boundaries open. //
		{
			CMatrix3D	rrho_cha(_N.x,_N.y,_N.z);
			double		AAbsRhoAmbMax = max(fabs(rrhoAmbMin),fabs(rrhoAmbMax));
			double		AAbsRhoChaMax = 0;
			for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
			{
				rrho_cha[ii][jj][kk] = foo::rhoijk(ii,jj,kk,pphi_cha,_d,_N);
				if(fabs(rrho_cha[ii][jj][kk])>= AAbsRhoChaMax)
					AAbsRhoChaMax = fabs(rrho_cha[ii][jj][kk]);
			};
			if(AAbsRhoChaMax>.01*AAbsRhoAmbMax)
			{
				for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
					/**********************************************************************/
					/* Check if point belongs to a boundary								  */
					/**********************************************************************/
					if(
						(kk == 0														) || // Down  //
						(kk == _N.z-1													) || // Up    //
						(ii == 0		&& kk > 0 && kk<_N.z-1							) || // Back  //
						(ii == _N.x-1	&& kk > 0 && kk<_N.z-1							) || // Front //
						(jj == 0		&& kk > 0 && kk<_N.z-1  && ii > 0 && ii<_N.x-1	) || // Left  //
						(jj == _N.y-1	&& kk > 0 && kk<_N.z-1  && ii > 0 && ii<_N.x-1	)    // Right //
						)
					{
						pphi_cha[ii][jj][kk] = 0;
						for(int is=0 ; is<_N.x ; is++) for(int js=0 ; js<_N.y ; js++) for(int ks=0 ; ks<_N.z ; ks++)
							if(!(is == ii && js == jj && ks == kk))
								if(fabs(rrho_cha[is][js][ks]) > .01*AAbsRhoAmbMax)
									pphi_cha[ii][jj][kk] +=
										_d.x*_d.y*_d.z/(4*PMC.pi*PMC.eps0) *
										rrho_cha[is][js][ks]/sqrt(pow( (ii-is)*_d.x,2 ) + pow( (jj-js)*_d.y,2 ) + pow( (kk-ks)*_d.z,2 ));
					};
				_EndTime = clock();
				_RunTime = _EndTime - _StartTime;
                msg = "ii:\t Run time for BC update: " + to_string((double)_RunTime/CLOCKS_PER_SEC)  + " s.\n";
                printf("%s",msg.c_str());
			}
			else printf("ii:\t No need for BC update.\n");
		};
	}
	else printf("ii:\t No update done.\n");
}



/**************************************************************************************/
void BC::AddUniformE(bool uuniformE, const double& EEo, CMatrix3D& _phi, CMatrix3D& UUn, const SizeDomain& LL, const ResGrid& _d, const SizeGrid& _N)
{
	clock_t _StartTime = clock();
	if(uuniformE == true)
	{
		printf("Uniform E-field Added\n");
		/******************************************************************************/
		/* Upper and Lower boundaries												  */
		/******************************************************************************/
		for(int ii=0 ; ii<_N.x ; ii++)
			for(int jj=0 ; jj<_N.y ; jj++)
			{
				UUn(ii,jj,0)		= 1;
				_phi(ii,jj,0)		+= 0;
				UUn(ii,jj,_N.z-1)	= 1;
				_phi(ii,jj,_N.z-1)	+= -EEo*(_d.z*(_N.z-1));
			}

		/******************************************************************************/
		/* Front and Back boundaries												  */
		/******************************************************************************/
		for(int jj=0 ; jj<_N.y ; jj++)
			for(int kk=1 ; kk<_N.z-1 ; kk++)
			{
				UUn(0,jj,kk)		= 1;
				_phi(0,jj,kk)		+= -EEo*(kk*_d.z);;
				UUn(_N.x-1,jj,kk)	= 1;
				_phi(_N.x-1,jj,kk)	+= -EEo*(kk*_d.z);
			}

		/******************************************************************************/
		/* Left and Right boundaries												  */
		/******************************************************************************/
		for(int ii=1 ; ii<_N.x-1 ; ii++)
			for(int kk=1 ; kk<_N.z-1 ; kk++)
			{
				UUn(ii,0,kk)		= 1;
				_phi(ii,0,kk)		+= -EEo*(kk*_d.z);;
				UUn(ii,_N.y-1,kk)	= 1;
				_phi(ii,_N.y-1,kk)	+= -EEo*(kk*_d.z);
			}
	}
	else printf("No Uniform E-field Applied.\n");

	clock_t _EndTime = clock();
	clock_t _RunTime = _EndTime - _StartTime;
	printf("Run time for application of BC: %lf s.\n",(double)_RunTime/CLOCKS_PER_SEC);
}
/**************************************************************************************/
