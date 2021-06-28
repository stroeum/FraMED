/*
 *  UsefulFunctions.cpp
 *  Created by Jeremy Riousset on 10/25/07.
 */

#include "UsefulFunctions.h"

/**************************************************************************************/
/* For Electric Field																  */
/**************************************************************************************/
CMatrix1D foo::Eijk(int ii, int jj, int kk,const CMatrix3D& _phi,const ResGrid& _d, const SizeGrid& _N)
{
	CMatrix1D EE(4);
	// x component //
	if (ii!=0 && ii!=_N.x-1)	EE[1] = - (_phi[ii+1][jj][kk]-_phi[ii-1][jj][kk]) / (2*_d.x);
	else if (ii==0)				EE[1] = - (_phi[ii+1][jj][kk]-_phi[ii][jj][kk]) / _d.x;
	else if (ii==_N.x-1)		EE[1] = - (_phi[ii][jj][kk]-_phi[ii-1][jj][kk]) / _d.x;

	// y component //
	if (jj!=0 && jj!=_N.y-1)	EE[2] = - (_phi[ii][jj+1][kk]-_phi[ii][jj-1][kk]) / (2*_d.y);
	else if (jj==0)				EE[2] = - (_phi[ii][jj+1][kk]-_phi[ii][jj][kk]) / _d.y;
	else if (jj==_N.y-1)		EE[2] = - (_phi[ii][jj][kk]-_phi[ii][jj-1][kk]) / _d.y;

	// z component //
	if (kk!=0 && kk!=_N.z-1)	EE[3] = - (_phi[ii][jj][kk+1]-_phi[ii][jj][kk-1]) / (2*_d.z);
	else if (kk==0)				EE[3] = - (_phi[ii][jj][kk+1]-_phi[ii][jj][kk]) / _d.z;
	else if (kk==_N.z-1)		EE[3] = - (_phi[ii][jj][kk]-_phi[ii][jj][kk-1]) / _d.z;

	// magnitude //
    EE[0] = sqrt(pow(EE[1],2)+pow(EE[2],2)+pow(EE[3],2));
	return EE;
};

CMatrix3D foo::GlobalE(const CMatrix3D& _phi,const ResGrid& _d, const SizeGrid& _N, const int nn)
{
	double		EEx(0), EEy(0), EEz(0);
	CMatrix3D	EE(_N.x,_N.y,_N.z);
	CMatrix3D	Ex(_N.x,_N.y,_N.z);
	CMatrix3D	Ey(_N.x,_N.y,_N.z);
	CMatrix3D	Ez(_N.x,_N.y,_N.z);
    char		strEx3D[50];
    char		strEy3D[50];
    char		strEz3D[50];
    sprintf(strEx3D,"Ex3d%d.dat", nn);
    sprintf(strEy3D,"Ey3d%d.dat", nn);
    sprintf(strEz3D,"Ez3d%d.dat", nn);

	for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
	{
		// x component //
		if (ii!=0 && ii!=_N.x-1)	EEx = - (_phi[ii+1][jj][kk]-_phi[ii-1][jj][kk]) / (2*_d.x);
		else if (ii==0)				EEx = - (_phi[ii+1][jj][kk]-_phi[ii][jj][kk]) / _d.x;
		else if (ii==_N.x-1)		EEx = - (_phi[ii][jj][kk]-_phi[ii-1][jj][kk]) / _d.x;

		// y component //
		if (jj!=0 && jj!=_N.y-1)	EEy = - (_phi[ii][jj+1][kk]-_phi[ii][jj-1][kk]) / (2*_d.y);
		else if (jj==0)				EEy = - (_phi[ii][jj+1][kk]-_phi[ii][jj][kk]) / _d.y;
		else if (jj==_N.y-1)		EEy = - (_phi[ii][jj][kk]-_phi[ii][jj-1][kk]) / _d.y;

		// z component //
		if (kk!=0 && kk!=_N.z-1)	EEz = - (_phi[ii][jj][kk+1]-_phi[ii][jj][kk-1]) / (2*_d.z);
		else if (kk==0)				EEz = - (_phi[ii][jj][kk+1]-_phi[ii][jj][kk]) / _d.z;
		else if (kk==_N.z-1)		EEz = - (_phi[ii][jj][kk]-_phi[ii][jj][kk-1]) / _d.z;

		EE[ii][jj][kk] = sqrt(pow(EEx,2) + pow(EEy,2) + pow(EEz,2));
		Ex[ii][jj][kk] = EEx;
		Ey[ii][jj][kk] = EEy;
		Ez[ii][jj][kk] = EEz;
	};
    if(nn==-1) {
        IO::write(Ex,(char*)"Ex3d.dat");
        IO::write(Ey,(char*)"Ey3d.dat");
        IO::write(Ez,(char*)"Ez3d.dat");
    }else{
        IO::write(Ex,strEx3D);
        IO::write(Ey,strEy3D);
        IO::write(Ez,strEz3D);
    }

	return EE;
};

CMatrix1D foo::eFieldFlux(const CMatrix3D& _phi, const ResGrid& _d, const SizeGrid& _N)
{
	int nn = 1;
	CMatrix1D eeFlux(8);
	for(int ii=0+nn ; ii<_N.x-nn ; ii++) for(int jj=0+nn ; jj<_N.y-nn ; jj++) for(int kk=0+nn ; kk<_N.z-nn ; kk++)
	{
		if (ii==0+nn		&& jj==0+nn			&& kk==0+nn		)
		{
			eeFlux[1] += -Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z/4;
			eeFlux[2] += -Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x/4;
			eeFlux[3] += -Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y/4;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/8;
		}
		if (ii==_N.x-1-nn	&& jj==0+nn			&& kk==0+nn		)
		{
			eeFlux[4] +=  Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z/4;
			eeFlux[2] += -Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x/4;
			eeFlux[3] += -Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y/4;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/8;
		}
		if (ii==0+nn		&& jj==_N.y-1-nn	&& kk==0+nn		)
		{
			eeFlux[1] += -Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z/4;
			eeFlux[5] +=  Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x/4;
			eeFlux[3] += -Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y/4;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/8;
		}
		if (ii==0+nn		&& jj==0+nn			&& kk==_N.z-1-nn)
		{
			eeFlux[1] += -Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z/4;
			eeFlux[2] += -Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x/4;
			eeFlux[6] +=  Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y/4;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/8;
		}
		if (ii==_N.x-1-nn	&& jj==_N.y-1-nn	&& kk==0+nn		)
		{
			eeFlux[4] +=  Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z/4;
			eeFlux[5] +=  Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x/4;
			eeFlux[3] += -Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y/4;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/8;
		}
		if (ii==0+nn		&& jj==_N.y-1-nn	&& kk==_N.z-1-nn)
		{
			eeFlux[1] += -Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z/4;
			eeFlux[5] +=  Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x/4;
			eeFlux[6] +=  Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y/4;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/8;
		}
		if (ii==_N.x-1-nn	&& jj==0+nn		&& kk==_N.z-1-nn	)
		{
			eeFlux[4] +=  Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z/4;
			eeFlux[2] += -Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x/4;
			eeFlux[6] +=  Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y/4;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/8;
		}
		if (ii==_N.x-1-nn	&& jj==_N.y-1-nn	&& kk==_N.z-1-nn)
		{
			eeFlux[4] +=  Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z/4;
			eeFlux[5] +=  Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x/4;
			eeFlux[6] +=  Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y/4;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/8;
		}
		if ( 0+nn< ii		&& ii <_N.x-1-nn	&& jj==0+nn			&& kk==0+nn		)
		{
			eeFlux[2] += -Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x/2;
			eeFlux[3] += -Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y/2;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/4;
		}
		if ( 0+nn< ii		&& ii <_N.x-1-nn	&& jj==_N.y-1-nn	&& kk==0+nn		)
		{
			eeFlux[5] +=  Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x/2;
			eeFlux[3] += -Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y/2;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/4;
		}
		if ( 0+nn< ii		&& ii <_N.x-1-nn	&& jj==0+nn			&& kk==_N.z-1-nn)
		{
			eeFlux[2] += -Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x/2;
			eeFlux[6] +=  Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y/2;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/4;
		}
		if ( 0+nn< ii		&& ii <_N.x-1-nn	&& jj==_N.y-1-nn	&& kk==_N.z-1-nn)
		{
			eeFlux[5] +=  Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x/2;
			eeFlux[6] +=  Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y/2;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/4;
		}
		if (ii==0+nn		&&  0+nn< jj		&& jj <_N.y-1-nn	&& kk==0+nn		)
		{
			eeFlux[1] += -Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z/2;
			eeFlux[3] += -Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y/2;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/4;
		}
		if (ii==_N.x-1-nn	&&  0+nn< jj		&& jj <_N.y-1-nn	&& kk==0+nn		)
		{
			eeFlux[4] +=  Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z/2;
			eeFlux[3] += -Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y/2;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/4;
		}
		if (ii==0+nn		&&  0+nn< jj		&& jj <_N.y-1-nn	&& kk==_N.z-1-nn)
		{
			eeFlux[1] += -Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z/2;
			eeFlux[6] +=  Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y/2;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/4;
		}
		if (ii==_N.x-1-nn	&&  0+nn< jj		&& jj <_N.y-1-nn	&& kk==_N.z-1-nn)
		{
			eeFlux[4] +=  Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z/2;
			eeFlux[6] +=  Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y/2;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/4;
		}
		if (ii==0+nn		&&  jj==0+nn		&&  0+nn< kk		&& kk <_N.z-1-nn)
		{
			eeFlux[1] += -Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z/2;
			eeFlux[2] += -Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x/2;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/4;
		}
		if (ii==_N.x-1-nn	&&  jj==0+nn		&&  0+nn< kk		&& kk <_N.z-1-nn)
		{
			eeFlux[4] +=  Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z/2;
			eeFlux[2] += -Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x/2;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/4;
		}
		if (ii==0+nn		&&  jj==_N.y-1-nn	&&  0+nn< kk		&& kk <_N.z-1-nn)
		{
			eeFlux[1] += -Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z/2;
			eeFlux[5] +=  Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x/2;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/4;
		}
		if (ii==_N.x-1-nn	&&  jj==_N.y-1-nn	&&  0+nn< kk		&& kk <_N.z-1-nn)
		{
			eeFlux[4] +=  Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z/2;
			eeFlux[5] +=  Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x/2;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/4;
		}

		if ( 0+nn< ii		&& ii <_N.x-1-nn	&&  0+nn< jj		&& jj <_N.y-1-nn	&& kk==0+nn		)
		{
			eeFlux[3] += -Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/2;
		}
		if ( 0+nn< ii		&& ii <_N.x-1-nn	&&  0+nn< jj		&& jj <_N.y-1-nn	&& kk==_N.z-1-nn)
		{
			eeFlux[6] +=  Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/2;
		}
		if ( 0+nn< ii		&& ii <_N.x-1-nn	&& jj==0+nn			&&  0+nn< kk		&& kk <_N.z-1-nn)
		{
			eeFlux[2] += -Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/2;
		}
		if ( 0+nn< ii		&& ii <_N.x-1-nn	&& jj==_N.y-1-nn	&&  0+nn< kk		&& kk <_N.z-1-nn)
		{
			eeFlux[5] +=  Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/2;
		}
		if (ii==0+nn		&&  0+nn< jj		&& jj <_N.y-1-nn	&&  0+nn< kk		&& kk <_N.z-1-nn)
		{
			eeFlux[1] += -Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/2;
		}
		if (ii==_N.x-1-nn	&&  0+nn< jj		&& jj <_N.y-1-nn	&&  0+nn< kk		&& kk <_N.z-1-nn)
		{
			eeFlux[4] +=  Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z;
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z/2;
		}
		if ( 0+nn< ii		&& ii <_N.x-1-nn	&&  0+nn< jj		&& jj <_N.y-1-nn	&&  0+nn< kk		&& kk <_N.z-1-nn)
			eeFlux[7] += rhoijk(ii,jj,kk,_phi,_d,_N)*_d.x*_d.y*_d.z;
		/*
		 if(ii == 0)		eeFlux[0] += -Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z;
		 if(ii == _N.x-1)	eeFlux[0] +=  Eijk(ii,jj,kk,_phi,_d,_N)[1]*_d.y*_d.z;
		 if(jj == 0)		eeFlux[0] += -Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x;
		 if(jj == _N.y-1)	eeFlux[0] +=  Eijk(ii,jj,kk,_phi,_d,_N)[2]*_d.z*_d.x;
		 if(kk == 0)		eeFlux[0] += -Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y;
		 if(kk == _N.z-1)	eeFlux[0] +=  Eijk(ii,jj,kk,_phi,_d,_N)[3]*_d.x*_d.y;
		 */
	};
	eeFlux[1]*= eps0; eeFlux[2]*= eps0; eeFlux[3]*= eps0;
	eeFlux[4]*= eps0; eeFlux[5]*= eps0;	eeFlux[6]*= eps0;
	eeFlux[0] = eeFlux[1] + eeFlux[2] + eeFlux[3] + eeFlux[4] + eeFlux[5] + eeFlux[6];
	return eeFlux;
}
/**************************************************************************************/

/**************************************************************************************/
/* For Charge																		  */
/**************************************************************************************/
double foo::rhoijk(int ii, int jj, int kk, const CMatrix3D& _phi, const ResGrid& _d, const SizeGrid& _N)
{

	double d2x=0;
	double d2y=0;
	double d2z=0;
	if (ii == 0)			d2x = (_phi(ii+2, jj , kk )-2*_phi(ii+1, jj , kk )+_phi( ii , jj , kk ))/(2*pow(_d.x,2));
	if (0<ii && ii<_N.x-1) 	d2x = (_phi(ii+1, jj , kk )-2*_phi( ii , jj , kk )+_phi(ii-1, jj , kk ))/   pow(_d.x,2) ;
	if (ii == _N.x-1)		d2x = (_phi( ii , jj , kk )-2*_phi(ii-1, jj , kk )+_phi(ii-2, jj , kk ))/(2*pow(_d.x,2));

	if (jj == 0)			d2y = (_phi( ii ,jj+2, kk )-2*_phi( ii ,jj+1, kk )+_phi( ii , jj , kk ))/(2*pow(_d.y,2));
	if (0<jj && jj<_N.y-1) 	d2y = (_phi( ii ,jj+1, kk )-2*_phi( ii , jj , kk )+_phi( ii ,jj-1, kk ))/   pow(_d.y,2) ;
	if (jj == _N.y-1)		d2y = (_phi( ii , jj , kk )-2*_phi( ii ,jj-1, kk )+_phi( ii ,jj-2, kk ))/(2*pow(_d.y,2));

	if (kk == 0)			d2z = (_phi( ii , jj ,kk+2)-2*_phi( ii , jj ,kk+1)+_phi( ii , jj , kk ))/(2*pow(_d.z,2));
	if (0<kk && kk<_N.z-1) 	d2z = (_phi( ii , jj ,kk+1)-2*_phi( ii , jj , kk )+_phi( ii , jj ,kk-1))/   pow(_d.z,2) ;
	if (kk == _N.z-1)		d2z = (_phi( ii , jj , kk )-2*_phi( ii , jj ,kk-1)+_phi( ii , jj ,kk-2))/(2*pow(_d.z,2));

	if(ii!=0 && ii!=_N.x-1 && jj!=0 && jj!=_N.y-1 && kk!=0 && kk!=_N.z-1)
		return  -eps0*(d2x+d2y+d2z);
	else
		return 0;
}

CMatrix3D foo::Globalrho(const CMatrix3D& _phi, const ResGrid& _d, const SizeGrid& _N)
{
	CMatrix3D	_rho(_N.x,_N.y,_N.z);
	for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
		_rho[ii][jj][kk]=rhoijk(ii,jj,kk,_phi,_d,_N);
	return _rho;
}

double foo::ChannelCharge(const CMatrix3D& _rho, const CMatrix3D& UUn, const ResGrid& _d, const SizeGrid& _N)
{
	double QQ=0;
	for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
		//	for(int ii=1 ; ii<_N.x-1 ; ii++) for(int jj=1 ; jj<_N.y-1 ; jj++) for(int kk=1 ; kk<_N.z-1 ; kk++)
		if(UUn[ii][jj][kk] != 0)
			QQ += _rho[ii][jj][kk]*_d.x*_d.y*_d.z;
	return QQ;
};

double foo::ChannelChargePositive(const CMatrix3D& _rho, const CMatrix3D& UUn, const ResGrid& _d, const SizeGrid& _N)
{
	double QQ_plus=0;
	double rro;

	for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
		//	for(int ii=1 ; ii<_N.x-1 ; ii++) for(int jj=1 ; jj<_N.y-1 ; jj++) for(int kk=1 ; kk<_N.z-1 ; kk++)
	{
		rro = _rho[ii][jj][kk];
		if(UUn[ii][jj][kk] != 0 && rro>=0) QQ_plus += rro*_d.x*_d.y*_d.z;
	}
		return QQ_plus;
};

double foo::ChannelChargeNegative(const CMatrix3D& _rho, const CMatrix3D& UUn, const ResGrid& _d, const SizeGrid& _N)
{
	double QQ_minus=0;
	double rro;

	for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
		//	for(int ii=1 ; ii<_N.x-1 ; ii++) for(int jj=1 ; jj<_N.y-1 ; jj++) for(int kk=1 ; kk<_N.z-1 ; kk++)
	{
		rro = _rho[ii][jj][kk];
		if(UUn[ii][jj][kk] == 1 && rro<=0) QQ_minus += rro*_d.x*_d.y*_d.z;
	};
	return QQ_minus;
};

double foo::TotalCharge(const CMatrix3D& _rho, const ResGrid& _d, const SizeGrid& _N)
{
	double QQ=0;
	double ddV = _d.x*_d.y*_d.z;
	for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
		//	for(int ii=1 ; ii<_N.x-1 ; ii++) for(int jj=1 ; jj<_N.y-1 ; jj++) for(int kk=1 ; kk<_N.z-1 ; kk++)
		QQ+= _rho[ii][jj][kk]*ddV;
	return QQ;
};

CMatrix1D foo::ChannelLinearDensity(const CMatrix3D& _rho, const CMatrix3D& UUn, const ResGrid& _d, const SizeGrid& _N)
{
	cout<<"*** Derivation valid for a single link in z-direction. ***\n";

	CMatrix1D	rrhol(_N.z);
	int			iic(0), jjc(0), kkcStart, kkcEnd;

	for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
		//	for(int ii=1 ; ii<_N.x-1 ; ii++) for(int jj=1 ; jj<_N.y-1 ; jj++) for(int kk=1 ; kk<_N.z-1 ; kk++)
		if(UUn[ii][jj][kk] == 1)
		{
			iic=ii;
			jjc=jj;
			if(kk>0 && UUn[ii][jj][kk-1] == 0)	kkcStart=kk;
			if(kk<_N.z-1 && UUn[ii][jj][kk+1] == 0)	kkcEnd=kk;
		};

	// solution 1 : integrate on a slice of the simulation box //
	/*
		for(int kk=0 ; kk<_N.z ; kk++)
	 {
			rrhol(kk) = 0;
			for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++)
				rrhol(kk) += _rho[ii][jj][kk]*_d.x*_d.y;
	 }
	 */
	// solution 2 : assume all charge is concentrated in the channel //
	for(int kk=0 ; kk<_N.z ; kk++) rrhol(kk) = _rho[iic][jjc][kk]*_d.x*_d.y;

	return rrhol;
}
/**************************************************************************************/

/**************************************************************************************/
/* For Dipole Moment																  */
/**************************************************************************************/
Vector foo::DipoleMoment(double& CCarriedCharge, const CMatrix3D& _phi, const CMatrix3D& UUn, const SizeDomain& LL, const ResGrid& _d, const SizeGrid& _N)
{
	Vector pp;
	double XXc(LL.x/2), YYc(LL.y/2), ZZc(0);	// define origine of the simulation domain
	double ddq		= 0;						// elementary charge
	pp.x			= 0;
	pp.y			= 0;
	pp.z			= 0;
	CCarriedCharge	= 0;

	for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
		//	for(int ii=1 ; ii<_N.x-1 ; ii++) for(int jj=1 ; jj<_N.y-1 ; jj++) for(int kk=1 ; kk<_N.z-1 ; kk++)
    {
		if(UUn[ii][jj][kk] == 1)
		{
			ddq = foo::rhoijk(ii,jj,kk, _phi,_d,_N)*_d.x*_d.y*_d.z;
			pp.x += ddq*(ii*_d.x-XXc);
			pp.y += ddq*(jj*_d.y-YYc);
			pp.z += ddq*(kk*_d.z-ZZc);

			if(ddq >=0) CCarriedCharge += ddq;
        }}
	return pp;
}
/**************************************************************************************/

/**************************************************************************************/
/* For E, phi, rho																	  */
/**************************************************************************************/
bool foo::isfinite(const CMatrix3D& MM, const SizeGrid& _N)
{
	bool flag = true;
	for(int ii=0 ; ii<_N.x ; ii++) for(int jj=0 ; jj<_N.y ; jj++) for(int kk=0 ; kk<_N.z ; kk++)
		if (isnan(MM(ii,jj,kk)) || isinf(MM(ii,jj,kk)))
		{
			cout<<"!!! Error inproper value.\n Break at point: ["<<ii<<" "<<jj<<" "<<kk<<"].\n"<<endl;
			flag = false;
			break;
		};
	return flag;
}
/**************************************************************************************/
