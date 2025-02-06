/*
 *  UsefulFunctions.cpp
 *  Created by Jeremy Riousset on 10/25/07.
 */

#include "UsefulFunctions.h"

/**************************************************************************************/
/* For Electric Field																  */
/**************************************************************************************/
CMatrix1D foo::Eijk(int ii, int jj, int kk,const CMatrix3D& phi,const ResGrid& d, const SizeGrid& N)
{
	CMatrix1D EE(4);
	// x component //
	if (ii!=0 && ii!=N.x-1)	EE[1] = - (phi[ii+1][jj][kk]-phi[ii-1][jj][kk]) / (2*d.x);
	else if (ii==0)			EE[1] = - (phi[ii+1][jj][kk]-phi[ii][jj][kk])   / d.x;
	else if (ii==N.x-1)		EE[1] = - (phi[ii][jj][kk]-phi[ii-1][jj][kk])   / d.x;

	// y component //
	if (jj!=0 && jj!=N.y-1)	EE[2] = - (phi[ii][jj+1][kk]-phi[ii][jj-1][kk]) / (2*d.y);
	else if (jj==0)			EE[2] = - (phi[ii][jj+1][kk]-phi[ii][jj][kk])   / d.y;
	else if (jj==N.y-1)		EE[2] = - (phi[ii][jj][kk]-phi[ii][jj-1][kk])   / d.y;

	// z component //
	if (kk!=0 && kk!=N.z-1)	EE[3] = - (phi[ii][jj][kk+1]-phi[ii][jj][kk-1]) / (2*d.z);
	else if (kk==0)			EE[3] = - (phi[ii][jj][kk+1]-phi[ii][jj][kk])   / d.z;
	else if (kk==N.z-1)		EE[3] = - (phi[ii][jj][kk]-phi[ii][jj][kk-1])   / d.z;

	// magnitude //
    EE[0] = sqrt(pow(EE[1],2)+pow(EE[2],2)+pow(EE[3],2));
	return EE;
};

CMatrix3D foo::GlobalE(const CMatrix3D& phi,const ResGrid& d, const SizeGrid& N, const int nn)
{
	double		_Ex(0), _Ey(0), _Ez(0);
	CMatrix3D	EE(N.x,N.y,N.z);
	CMatrix3D	Ex(N.x,N.y,N.z);
	CMatrix3D	Ey(N.x,N.y,N.z);
	CMatrix3D	Ez(N.x,N.y,N.z);
    char		_strEx3D[50];
    char		_strEy3D[50];
    char		_strEz3D[50];
    sprintf(_strEx3D,"Ex3d%d.dat", nn);
    sprintf(_strEy3D,"Ey3d%d.dat", nn);
    sprintf(_strEz3D,"Ez3d%d.dat", nn);

	for(int ii=0 ; ii<N.x ; ii++) for(int jj=0 ; jj<N.y ; jj++) for(int kk=0 ; kk<N.z ; kk++)
	{
		// x component //
		if (ii!=0 && ii!=N.x-1)	_Ex = - (phi[ii+1][jj][kk]-phi[ii-1][jj][kk])   / (2*d.x);
		else if (ii==0)			_Ex = - (phi[ii+1][jj][kk]-phi[ii][jj][kk])     / d.x;
		else if (ii==N.x-1)		_Ex = - (phi[ii][jj][kk]-phi[ii-1][jj][kk])     / d.x;

		// y component //
		if (jj!=0 && jj!=N.y-1)	_Ey = - (phi[ii][jj+1][kk]-phi[ii][jj-1][kk])   / (2*d.y);
		else if (jj==0)			_Ey = - (phi[ii][jj+1][kk]-phi[ii][jj][kk])     / d.y;
		else if (jj==N.y-1)		_Ey = - (phi[ii][jj][kk]-phi[ii][jj-1][kk])     / d.y;

		// z component //
		if (kk!=0 && kk!=N.z-1)	_Ez = - (phi[ii][jj][kk+1]-phi[ii][jj][kk-1])   / (2*d.z);
		else if (kk==0)			_Ez = - (phi[ii][jj][kk+1]-phi[ii][jj][kk])     / d.z;
		else if (kk==N.z-1)		_Ez = - (phi[ii][jj][kk]-phi[ii][jj][kk-1])     / d.z;

		EE[ii][jj][kk] = sqrt(pow(_Ex,2) + pow(_Ey,2) + pow(_Ez,2));
		Ex[ii][jj][kk] = _Ex;
		Ey[ii][jj][kk] = _Ey;
		Ez[ii][jj][kk] = _Ez;
	};
    if(nn==-2){
		// Saves the ambient electric fields
        IO::write(Ex,(char*)"Ex3dAmb.dat");
        IO::write(Ey,(char*)"Ey3dAmb.dat");
        IO::write(Ez,(char*)"Ez3dAmb.dat");
    }else if(nn==-1){						
	    // Saves the most recent/final electric fields
		IO::write(Ex,(char*)"Ex3d.dat");
        IO::write(Ey,(char*)"Ey3d.dat");
        IO::write(Ez,(char*)"Ez3d.dat");
	}else{
        // Saves the electric fields on a Var::step3d basis
		IO::write(Ex,_strEx3D);
        IO::write(Ey,_strEy3D);
        IO::write(Ez,_strEz3D);
    }

	return EE;
};

CMatrix1D foo::eFieldFlux(const CMatrix3D& phi, const ResGrid& d, const SizeGrid& N)
{
	int nn = 1;
	CMatrix1D eFlux(8);
	for(int ii=0+nn ; ii<N.x-nn ; ii++) for(int jj=0+nn ; jj<N.y-nn ; jj++) for(int kk=0+nn ; kk<N.z-nn ; kk++)
	{
		if (ii==0+nn		&& jj==0+nn			&& kk==0+nn		)
		{
			eFlux[1] += -Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z/4;
			eFlux[2] += -Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x/4;
			eFlux[3] += -Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y/4;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/8;
		}
		if (ii==N.x-1-nn	&& jj==0+nn			&& kk==0+nn		)
		{
			eFlux[4] +=  Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z/4;
			eFlux[2] += -Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x/4;
			eFlux[3] += -Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y/4;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/8;
		}
		if (ii==0+nn		&& jj==N.y-1-nn	&& kk==0+nn		)
		{
			eFlux[1] += -Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z/4;
			eFlux[5] +=  Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x/4;
			eFlux[3] += -Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y/4;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/8;
		}
		if (ii==0+nn		&& jj==0+nn			&& kk==N.z-1-nn)
		{
			eFlux[1] += -Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z/4;
			eFlux[2] += -Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x/4;
			eFlux[6] +=  Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y/4;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/8;
		}
		if (ii==N.x-1-nn	&& jj==N.y-1-nn	&& kk==0+nn		)
		{
			eFlux[4] +=  Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z/4;
			eFlux[5] +=  Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x/4;
			eFlux[3] += -Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y/4;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/8;
		}
		if (ii==0+nn		&& jj==N.y-1-nn	&& kk==N.z-1-nn)
		{
			eFlux[1] += -Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z/4;
			eFlux[5] +=  Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x/4;
			eFlux[6] +=  Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y/4;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/8;
		}
		if (ii==N.x-1-nn	&& jj==0+nn		&& kk==N.z-1-nn	)
		{
			eFlux[4] +=  Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z/4;
			eFlux[2] += -Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x/4;
			eFlux[6] +=  Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y/4;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/8;
		}
		if (ii==N.x-1-nn	&& jj==N.y-1-nn	&& kk==N.z-1-nn)
		{
			eFlux[4] +=  Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z/4;
			eFlux[5] +=  Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x/4;
			eFlux[6] +=  Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y/4;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/8;
		}
		if ( 0+nn< ii		&& ii <N.x-1-nn	&& jj==0+nn		&& kk==0+nn		)
		{
			eFlux[2] += -Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x/2;
			eFlux[3] += -Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y/2;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/4;
		}
		if ( 0+nn< ii		&& ii <N.x-1-nn	&& jj==N.y-1-nn	&& kk==0+nn		)
		{
			eFlux[5] +=  Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x/2;
			eFlux[3] += -Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y/2;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/4;
		}
		if ( 0+nn< ii		&& ii <N.x-1-nn	&& jj==0+nn		&& kk==N.z-1-nn)
		{
			eFlux[2] += -Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x/2;
			eFlux[6] +=  Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y/2;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/4;
		}
		if ( 0+nn< ii		&& ii <N.x-1-nn	&& jj==N.y-1-nn	&& kk==N.z-1-nn)
		{
			eFlux[5] +=  Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x/2;
			eFlux[6] +=  Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y/2;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/4;
		}
		if (ii==0+nn		&&  0+nn< jj		&& jj <N.y-1-nn	&& kk==0+nn		)
		{
			eFlux[1] += -Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z/2;
			eFlux[3] += -Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y/2;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/4;
		}
		if (ii==N.x-1-nn	&&  0+nn< jj		&& jj <N.y-1-nn	&& kk==0+nn		)
		{
			eFlux[4] +=  Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z/2;
			eFlux[3] += -Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y/2;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/4;
		}
		if (ii==0+nn		&&  0+nn< jj		&& jj <N.y-1-nn	&& kk==N.z-1-nn)
		{
			eFlux[1] += -Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z/2;
			eFlux[6] +=  Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y/2;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/4;
		}
		if (ii==N.x-1-nn	&&  0+nn< jj		&& jj <N.y-1-nn	&& kk==N.z-1-nn)
		{
			eFlux[4] +=  Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z/2;
			eFlux[6] +=  Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y/2;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/4;
		}
		if (ii==0+nn		&&  jj==0+nn		&&  0+nn< kk	&& kk <N.z-1-nn)
		{
			eFlux[1] += -Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z/2;
			eFlux[2] += -Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x/2;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/4;
		}
		if (ii==N.x-1-nn	&&  jj==0+nn		&&  0+nn< kk	&& kk <N.z-1-nn)
		{
			eFlux[4] +=  Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z/2;
			eFlux[2] += -Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x/2;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/4;
		}
		if (ii==0+nn		&&  jj==N.y-1-nn	&&  0+nn< kk	&& kk <N.z-1-nn)
		{
			eFlux[1] += -Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z/2;
			eFlux[5] +=  Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x/2;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/4;
		}
		if (ii==N.x-1-nn	&&  jj==N.y-1-nn	&&  0+nn< kk	&& kk <N.z-1-nn)
		{
			eFlux[4] +=  Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z/2;
			eFlux[5] +=  Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x/2;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/4;
		}

		if ( 0+nn< ii		&& ii <N.x-1-nn	&&  0+nn< jj		&& jj <N.y-1-nn	&& kk==0+nn		)
		{
			eFlux[3] += -Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/2;
		}
		if ( 0+nn< ii		&& ii <N.x-1-nn	&&  0+nn< jj		&& jj <N.y-1-nn	&& kk==N.z-1-nn)
		{
			eFlux[6] +=  Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/2;
		}
		if ( 0+nn< ii		&& ii <N.x-1-nn	&& jj==0+nn			&&  0+nn< kk	&& kk <N.z-1-nn)
		{
			eFlux[2] += -Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/2;
		}
		if ( 0+nn< ii		&& ii <N.x-1-nn	&& jj==N.y-1-nn     &&  0+nn< kk	&& kk <N.z-1-nn)
		{
			eFlux[5] +=  Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/2;
		}
		if (ii==0+nn		&&  0+nn< jj		&& jj <N.y-1-nn	&&  0+nn< kk	&& kk <N.z-1-nn)
		{
			eFlux[1] += -Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/2;
		}
		if (ii==N.x-1-nn	&&  0+nn< jj		&& jj <N.y-1-nn	&&  0+nn< kk	&& kk <N.z-1-nn)
		{
			eFlux[4] +=  Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z;
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z/2;
		}
		if ( 0+nn< ii		&& ii <N.x-1-nn     &&  0+nn< jj	&& jj <N.y-1-nn	&&  0+nn< kk		&& kk <N.z-1-nn)
			eFlux[7] += rhoijk(ii,jj,kk,phi,d,N)*d.x*d.y*d.z;
		/*
		 if(ii == 0)		eFlux[0] += -Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z;
		 if(ii == N.x-1)	eFlux[0] +=  Eijk(ii,jj,kk,phi,d,N)[1]*d.y*d.z;
		 if(jj == 0)		eFlux[0] += -Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x;
		 if(jj == N.y-1)	eFlux[0] +=  Eijk(ii,jj,kk,phi,d,N)[2]*d.z*d.x;
		 if(kk == 0)		eFlux[0] += -Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y;
		 if(kk == N.z-1)	eFlux[0] +=  Eijk(ii,jj,kk,phi,d,N)[3]*d.x*d.y;
		 */
	};
	eFlux[1]*= eps0; eFlux[2]*= eps0; eFlux[3]*= eps0;
	eFlux[4]*= eps0; eFlux[5]*= eps0;	eFlux[6]*= eps0;
	eFlux[0] = eFlux[1] + eFlux[2] + eFlux[3] + eFlux[4] + eFlux[5] + eFlux[6];
	return eFlux;
}
/**************************************************************************************/

/**************************************************************************************/
/* For Charge																		  */
/**************************************************************************************/
double foo::rhoijk(int ii, int jj, int kk, const CMatrix3D& phi, const ResGrid& d, const SizeGrid& N)
{

	double d2x=0;
	double d2y=0;
	double d2z=0;
	if (ii == 0)			d2x = (phi(ii+2, jj , kk )-2*phi(ii+1, jj , kk )+phi( ii , jj , kk ))/(2*pow(d.x,2));
	if (0<ii && ii<N.x-1) 	d2x = (phi(ii+1, jj , kk )-2*phi( ii , jj , kk )+phi(ii-1, jj , kk ))/   pow(d.x,2) ;
	if (ii == N.x-1)		d2x = (phi( ii , jj , kk )-2*phi(ii-1, jj , kk )+phi(ii-2, jj , kk ))/(2*pow(d.x,2));

	if (jj == 0)			d2y = (phi( ii ,jj+2, kk )-2*phi( ii ,jj+1, kk )+phi( ii , jj , kk ))/(2*pow(d.y,2));
	if (0<jj && jj<N.y-1) 	d2y = (phi( ii ,jj+1, kk )-2*phi( ii , jj , kk )+phi( ii ,jj-1, kk ))/   pow(d.y,2) ;
	if (jj == N.y-1)		d2y = (phi( ii , jj , kk )-2*phi( ii ,jj-1, kk )+phi( ii ,jj-2, kk ))/(2*pow(d.y,2));

	if (kk == 0)			d2z = (phi( ii , jj ,kk+2)-2*phi( ii , jj ,kk+1)+phi( ii , jj , kk ))/(2*pow(d.z,2));
	if (0<kk && kk<N.z-1) 	d2z = (phi( ii , jj ,kk+1)-2*phi( ii , jj , kk )+phi( ii , jj ,kk-1))/   pow(d.z,2) ;
	if (kk == N.z-1)		d2z = (phi( ii , jj , kk )-2*phi( ii , jj ,kk-1)+phi( ii , jj ,kk-2))/(2*pow(d.z,2));

	if(ii!=0 && ii!=N.x-1 && jj!=0 && jj!=N.y-1 && kk!=0 && kk!=N.z-1)
		return  -eps0*(d2x+d2y+d2z);
	else
		return 0;
}

CMatrix3D foo::Globalrho(const CMatrix3D& phi, const ResGrid& d, const SizeGrid& N)
{
	CMatrix3D	_rho(N.x,N.y,N.z);
	for(int ii=0 ; ii<N.x ; ii++) for(int jj=0 ; jj<N.y ; jj++) for(int kk=0 ; kk<N.z ; kk++)
		_rho[ii][jj][kk]=rhoijk(ii,jj,kk,phi,d,N);
	return _rho;
}

double foo::ChannelCharge(const CMatrix3D& _rho, const CMatrix3D& UUn, const ResGrid& d, const SizeGrid& N)
{
	double QQ=0;
	for(int ii=0 ; ii<N.x ; ii++) for(int jj=0 ; jj<N.y ; jj++) for(int kk=0 ; kk<N.z ; kk++)
		//	for(int ii=1 ; ii<N.x-1 ; ii++) for(int jj=1 ; jj<N.y-1 ; jj++) for(int kk=1 ; kk<N.z-1 ; kk++)
		if(UUn[ii][jj][kk] != 0)
			QQ += _rho[ii][jj][kk]*d.x*d.y*d.z;
	return QQ;
};

double foo::ChannelChargePositive(const CMatrix3D& _rho, const CMatrix3D& UUn, const ResGrid& d, const SizeGrid& N)
{
	double QQ_plus=0;
	double rro;

	for(int ii=0 ; ii<N.x ; ii++) for(int jj=0 ; jj<N.y ; jj++) for(int kk=0 ; kk<N.z ; kk++)
		//	for(int ii=1 ; ii<N.x-1 ; ii++) for(int jj=1 ; jj<N.y-1 ; jj++) for(int kk=1 ; kk<N.z-1 ; kk++)
	{
		rro = _rho[ii][jj][kk];
		if(UUn[ii][jj][kk] != 0 && rro>=0) QQ_plus += rro*d.x*d.y*d.z;
	}
		return QQ_plus;
};

double foo::ChannelChargeNegative(const CMatrix3D& _rho, const CMatrix3D& UUn, const ResGrid& d, const SizeGrid& N)
{
	double QQ_minus=0;
	double rro;

	for(int ii=0 ; ii<N.x ; ii++) for(int jj=0 ; jj<N.y ; jj++) for(int kk=0 ; kk<N.z ; kk++)
		//	for(int ii=1 ; ii<N.x-1 ; ii++) for(int jj=1 ; jj<N.y-1 ; jj++) for(int kk=1 ; kk<N.z-1 ; kk++)
	{
		rro = _rho[ii][jj][kk];
		if(UUn[ii][jj][kk] == 1 && rro<=0) QQ_minus += rro*d.x*d.y*d.z;
	};
	return QQ_minus;
};

double foo::TotalCharge(const CMatrix3D& _rho, const ResGrid& d, const SizeGrid& N)
{
	double QQ=0;
	double ddV = d.x*d.y*d.z;
	for(int ii=0 ; ii<N.x ; ii++) for(int jj=0 ; jj<N.y ; jj++) for(int kk=0 ; kk<N.z ; kk++)
		//	for(int ii=1 ; ii<N.x-1 ; ii++) for(int jj=1 ; jj<N.y-1 ; jj++) for(int kk=1 ; kk<N.z-1 ; kk++)
		QQ+= _rho[ii][jj][kk]*ddV;
	return QQ;
};

CMatrix1D foo::ChannelLinearDensity(const CMatrix3D& _rho, const CMatrix3D& UUn, const ResGrid& d, const SizeGrid& N)
{
	cout<<"*** Derivation valid for a single link in z-direction. ***\n";

	CMatrix1D	rrhol(N.z);
	int			iic(0), jjc(0), kkcStart, kkcEnd;

	for(int ii=0 ; ii<N.x ; ii++) for(int jj=0 ; jj<N.y ; jj++) for(int kk=0 ; kk<N.z ; kk++)
		//	for(int ii=1 ; ii<N.x-1 ; ii++) for(int jj=1 ; jj<N.y-1 ; jj++) for(int kk=1 ; kk<N.z-1 ; kk++)
		if(UUn[ii][jj][kk] == 1)
		{
			iic=ii;
			jjc=jj;
			if(kk>0 && UUn[ii][jj][kk-1] == 0)	kkcStart=kk;
			if(kk<N.z-1 && UUn[ii][jj][kk+1] == 0)	kkcEnd=kk;
		};

	// solution 1 : integrate on a slice of the simulation box //
	/*
		for(int kk=0 ; kk<N.z ; kk++)
	 {
			rrhol(kk) = 0;
			for(int ii=0 ; ii<N.x ; ii++) for(int jj=0 ; jj<N.y ; jj++)
				rrhol(kk) += _rho[ii][jj][kk]*d.x*d.y;
	 }
	 */
	// solution 2 : assume all charge is concentrated in the channel //
	for(int kk=0 ; kk<N.z ; kk++) rrhol(kk) = _rho[iic][jjc][kk]*d.x*d.y;

	return rrhol;
}
/**************************************************************************************/

/**************************************************************************************/
/* For Dipole Moment																  */
/**************************************************************************************/
Vector foo::DipoleMoment(double& CCarriedCharge, const CMatrix3D& phi, const CMatrix3D& UUn, const SizeDomain& LL, const ResGrid& d, const SizeGrid& N)
{
	Vector pp;
	double XXc(LL.x/2), YYc(LL.y/2), ZZc(0);	// define origine of the simulation domain
	double ddq		= 0;						// elementary charge
	pp.x			= 0;
	pp.y			= 0;
	pp.z			= 0;
	CCarriedCharge	= 0;

	for(int ii=0 ; ii<N.x ; ii++) for(int jj=0 ; jj<N.y ; jj++) for(int kk=0 ; kk<N.z ; kk++)
		//	for(int ii=1 ; ii<N.x-1 ; ii++) for(int jj=1 ; jj<N.y-1 ; jj++) for(int kk=1 ; kk<N.z-1 ; kk++)
    {
		if(UUn[ii][jj][kk] == 1)
		{
			ddq = foo::rhoijk(ii,jj,kk, phi,d,N)*d.x*d.y*d.z;
			pp.x += ddq*(ii*d.x-XXc);
			pp.y += ddq*(jj*d.y-YYc);
			pp.z += ddq*(kk*d.z-ZZc);

			if(ddq >=0) CCarriedCharge += ddq;
        }}
	return pp;
}
/**************************************************************************************/

/**************************************************************************************/
/* For E, phi, rho																	  */
/**************************************************************************************/
bool foo::isfinite(const CMatrix3D& MM, const SizeGrid& N)
{
	bool flag = true;
	for(int ii=0 ; ii<N.x ; ii++) for(int jj=0 ; jj<N.y ; jj++) for(int kk=0 ; kk<N.z ; kk++)
		if (isnan(MM(ii,jj,kk)) || isinf(MM(ii,jj,kk)))
		{
			cout<<"!!! Error inproper value.\n Break at point: ["<<ii<<" "<<jj<<" "<<kk<<"].\n"<<endl;
			flag = false;
			break;
		};
	return flag;
}
/**************************************************************************************/
