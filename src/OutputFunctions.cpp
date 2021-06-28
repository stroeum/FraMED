/* OutputFunctions.cpp */

#include "OutputFunctions.h"

/**************************************************************************************/
/* List of doubles (channel potential)												  */
/**************************************************************************************/
ListDouble read(string ss)
{
	char *fname = &ss[0];
	ListDouble LL;
	double xx;
	ifstream inFile;
	
	/**********************************************************************************/
	/* Open the file to be writen.													  */
	/**********************************************************************************/

	/* printf("Opening file '%s'\n",relativePath); */
	inFile.open(fname);
	if(!inFile)
		cerr<<"cannot open file";
	//nrerror("cannot open file"); 
	
	/**********************************************************************************/
	/* Write datas.																	  */
	/**********************************************************************************/
	while(inFile>>xx)
		LL.push_back(xx);
	inFile.close();
	return LL;
}
/**************************************************************************************/

/**************************************************************************************/
/* List of doubles (channel potential)												  */
/**************************************************************************************/
void write(ListDouble& LL, string ss)
{
	char *fname = &ss[0];
	
	ListDouble::iterator it;
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
		fprintf(file,"%f\n", *it);
	fclose(file);
}
/**************************************************************************************/

/**************************************************************************************/
/* List of Vectors (e.g., Dipole moments)											  */
/**************************************************************************************/
void write(ListVector& LL, string ss)
{
	char *fname = &ss[0];
	
	ListVector::iterator it;
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
		fprintf(file,"%f %f %f\n", it->x,it->y,it->z);
	}
	fclose(file);
}
/**************************************************************************************/

/**************************************************************************************/
/* List of CMatrix1D (e.g., Ez on the central vertical axis)						  */
/**************************************************************************************/
void write(ListCMatrix1D& LL, string ss)
{
	char *fname = &ss[0];
	
	ListCMatrix1D::iterator it;
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
		for (int ii=0 ; ii<it->getNbElem() ; ii++)
			fprintf(file,"%f ",it->getElem(ii));
		fprintf(file,"\n");
	}
	fclose(file);
}
/**************************************************************************************/

/**************************************************************************************/
/* For Electric Field																  */
/**************************************************************************************/
CMatrix1D Eijk(int ii, int jj, int kk,const CMatrix3D& pphi,const StepsSizes& dd, const BoxSteps& NN)
{
	CMatrix1D EE(4);
	// x component //
	if (ii!=0 && ii!=NN.x-1)	EE[1] = - (pphi[ii+1][jj][kk]-pphi[ii-1][jj][kk]) / (2*dd.x);
	else if (ii==0)				EE[1] = - (pphi[ii+1][jj][kk]-pphi[ii][jj][kk]) / dd.x;
	else if (ii==NN.x-1)		EE[1] = - (pphi[ii][jj][kk]-pphi[ii-1][jj][kk]) / dd.x; 
				
	// y component //
	if (jj!=0 && jj!=NN.y-1)	EE[2] = - (pphi[ii][jj+1][kk]-pphi[ii][jj-1][kk]) / (2*dd.y);
	else if (jj==0)				EE[2] = - (pphi[ii][jj+1][kk]-pphi[ii][jj][kk]) / dd.y;
	else if (jj==NN.y-1)		EE[2] = - (pphi[ii][jj][kk]-pphi[ii][jj-1][kk]) / dd.y;
	
	// z component //
	if (kk!=0 && kk!=NN.z-1)	EE[3] = - (pphi[ii][jj][kk+1]-pphi[ii][jj][kk-1]) / (2*dd.z);
	else if (kk==0)				EE[3] = - (pphi[ii][jj][kk+1]-pphi[ii][jj][kk]) / dd.z;
	else if (kk==NN.z-1)		EE[3] = - (pphi[ii][jj][kk]-pphi[ii][jj][kk-1]) / dd.z;
				
	// magnitude //
    EE[0] = sqrt(pow(EE[1],2)+pow(EE[2],2)+pow(EE[3],2)); 	
	return EE;
};

CMatrix3D GlobalE(const CMatrix3D& pphi,const StepsSizes& dd, const BoxSteps& NN)
{
	double		EEx(0), EEy(0), EEz(0);
	CMatrix3D	EE(NN.x,NN.y,NN.z);
	CMatrix3D	Ex(NN.x,NN.y,NN.z);
	CMatrix3D	Ey(NN.x,NN.y,NN.z);
	CMatrix3D	Ez(NN.x,NN.y,NN.z);
	
	for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
	{
		// x component //
		if (ii!=0 && ii!=NN.x-1)	EEx = - (pphi[ii+1][jj][kk]-pphi[ii-1][jj][kk]) / (2*dd.x);
		else if (ii==0)				EEx = - (pphi[ii+1][jj][kk]-pphi[ii][jj][kk]) / dd.x;
		else if (ii==NN.x-1)		EEx = - (pphi[ii][jj][kk]-pphi[ii-1][jj][kk]) / dd.x; 
		
		// y component //
		if (jj!=0 && jj!=NN.y-1)	EEy = - (pphi[ii][jj+1][kk]-pphi[ii][jj-1][kk]) / (2*dd.y);
		else if (jj==0)				EEy = - (pphi[ii][jj+1][kk]-pphi[ii][jj][kk]) / dd.y;
		else if (jj==NN.y-1)		EEy = - (pphi[ii][jj][kk]-pphi[ii][jj-1][kk]) / dd.y;
		
		// z component //
		if (kk!=0 && kk!=NN.z-1)	EEz = - (pphi[ii][jj][kk+1]-pphi[ii][jj][kk-1]) / (2*dd.z);
		else if (kk==0)				EEz = - (pphi[ii][jj][kk+1]-pphi[ii][jj][kk]) / dd.z;
		else if (kk==NN.z-1)		EEz = - (pphi[ii][jj][kk]-pphi[ii][jj][kk-1]) / dd.z;
		
		EE[ii][jj][kk] = sqrt(pow(EEx,2) + pow(EEy,2) + pow(EEz,2));
		Ex[ii][jj][kk] = EEx;
		Ey[ii][jj][kk] = EEy;
		Ez[ii][jj][kk] = EEz;
	};
	Ex.fwrite("results/Ex.dat");
	Ey.fwrite("results/Ey.dat");
	Ez.fwrite("results/Ez.dat");

	return EE;
};

CMatrix1D eFieldFlux(const CMatrix3D& pphi, const StepsSizes& dd, const BoxSteps& NN)
{
	int nn = 1;
	CMatrix1D eeFlux(8);
	for(int ii=0+nn ; ii<NN.x-nn ; ii++) for(int jj=0+nn ; jj<NN.y-nn ; jj++) for(int kk=0+nn ; kk<NN.z-nn ; kk++)
	{
		if (ii==0+nn		&& jj==0+nn			&& kk==0+nn		)
		{
			eeFlux[1] += -Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z/4;
			eeFlux[2] += -Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x/4;
			eeFlux[3] += -Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y/4;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/8;
		}
		if (ii==NN.x-1-nn	&& jj==0+nn			&& kk==0+nn		)
		{
			eeFlux[4] +=  Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z/4;
			eeFlux[2] += -Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x/4;
			eeFlux[3] += -Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y/4;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/8;
		}
		if (ii==0+nn		&& jj==NN.y-1-nn	&& kk==0+nn		)
		{
			eeFlux[1] += -Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z/4;
			eeFlux[5] +=  Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x/4;
			eeFlux[3] += -Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y/4;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/8;
		}
		if (ii==0+nn		&& jj==0+nn			&& kk==NN.z-1-nn)
		{
			eeFlux[1] += -Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z/4;
			eeFlux[2] += -Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x/4;
			eeFlux[6] +=  Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y/4;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/8;
		}
		if (ii==NN.x-1-nn	&& jj==NN.y-1-nn	&& kk==0+nn		)
		{
			eeFlux[4] +=  Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z/4;
			eeFlux[5] +=  Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x/4;
			eeFlux[3] += -Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y/4;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/8;
		}
		if (ii==0+nn		&& jj==NN.y-1-nn	&& kk==NN.z-1-nn)
		{
			eeFlux[1] += -Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z/4;
			eeFlux[5] +=  Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x/4;
			eeFlux[6] +=  Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y/4;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/8;
		}
		if (ii==NN.x-1-nn	&& jj==0+nn		&& kk==NN.z-1-nn	)
		{
			eeFlux[4] +=  Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z/4;
			eeFlux[2] += -Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x/4;
			eeFlux[6] +=  Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y/4;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/8;
		}
		if (ii==NN.x-1-nn	&& jj==NN.y-1-nn	&& kk==NN.z-1-nn)
		{
			eeFlux[4] +=  Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z/4;
			eeFlux[5] +=  Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x/4;
			eeFlux[6] +=  Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y/4;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/8;
		}
		if ( 0+nn< ii		&& ii <NN.x-1-nn	&& jj==0+nn			&& kk==0+nn		)
		{
			eeFlux[2] += -Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x/2;
			eeFlux[3] += -Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y/2;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/4;
		}
		if ( 0+nn< ii		&& ii <NN.x-1-nn	&& jj==NN.y-1-nn	&& kk==0+nn		)
		{
			eeFlux[5] +=  Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x/2;
			eeFlux[3] += -Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y/2;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/4;
		}
		if ( 0+nn< ii		&& ii <NN.x-1-nn	&& jj==0+nn			&& kk==NN.z-1-nn)
		{
			eeFlux[2] += -Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x/2;
			eeFlux[6] +=  Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y/2;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/4;
		}
		if ( 0+nn< ii		&& ii <NN.x-1-nn	&& jj==NN.y-1-nn	&& kk==NN.z-1-nn)
		{
			eeFlux[5] +=  Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x/2;
			eeFlux[6] +=  Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y/2; 
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/4;			
		}
		if (ii==0+nn		&&  0+nn< jj		&& jj <NN.y-1-nn	&& kk==0+nn		)
		{
			eeFlux[1] += -Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z/2;
			eeFlux[3] += -Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y/2;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/4;
		}
		if (ii==NN.x-1-nn	&&  0+nn< jj		&& jj <NN.y-1-nn	&& kk==0+nn		)
		{
			eeFlux[4] +=  Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z/2;
			eeFlux[3] += -Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y/2; 
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/4;			
		}
		if (ii==0+nn		&&  0+nn< jj		&& jj <NN.y-1-nn	&& kk==NN.z-1-nn)
		{
			eeFlux[1] += -Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z/2;
			eeFlux[6] +=  Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y/2;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/4;
		}
		if (ii==NN.x-1-nn	&&  0+nn< jj		&& jj <NN.y-1-nn	&& kk==NN.z-1-nn)
		{
			eeFlux[4] +=  Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z/2;
			eeFlux[6] +=  Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y/2;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/4;
		}		
		if (ii==0+nn		&&  jj==0+nn		&&  0+nn< kk		&& kk <NN.z-1-nn)
		{
			eeFlux[1] += -Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z/2;
			eeFlux[2] += -Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x/2;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/4;
		}
		if (ii==NN.x-1-nn	&&  jj==0+nn		&&  0+nn< kk		&& kk <NN.z-1-nn)
		{
			eeFlux[4] +=  Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z/2;
			eeFlux[2] += -Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x/2;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/4;
		}
		if (ii==0+nn		&&  jj==NN.y-1-nn	&&  0+nn< kk		&& kk <NN.z-1-nn)
		{
			eeFlux[1] += -Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z/2;
			eeFlux[5] +=  Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x/2;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/4;
		}
		if (ii==NN.x-1-nn	&&  jj==NN.y-1-nn	&&  0+nn< kk		&& kk <NN.z-1-nn)
		{
			eeFlux[4] +=  Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z/2;
			eeFlux[5] +=  Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x/2; 
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/4;
		}
		
		if ( 0+nn< ii		&& ii <NN.x-1-nn	&&  0+nn< jj		&& jj <NN.y-1-nn	&& kk==0+nn		)
		{	
			eeFlux[3] += -Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/2;
		}
		if ( 0+nn< ii		&& ii <NN.x-1-nn	&&  0+nn< jj		&& jj <NN.y-1-nn	&& kk==NN.z-1-nn)
		{	
			eeFlux[6] +=  Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/2;
		}
		if ( 0+nn< ii		&& ii <NN.x-1-nn	&& jj==0+nn			&&  0+nn< kk		&& kk <NN.z-1-nn)
		{	
			eeFlux[2] += -Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/2;
		}
		if ( 0+nn< ii		&& ii <NN.x-1-nn	&& jj==NN.y-1-nn	&&  0+nn< kk		&& kk <NN.z-1-nn)
		{	
			eeFlux[5] +=  Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/2;
		}
		if (ii==0+nn		&&  0+nn< jj		&& jj <NN.y-1-nn	&&  0+nn< kk		&& kk <NN.z-1-nn)
		{	
			eeFlux[1] += -Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/2;
		}
		if (ii==NN.x-1-nn	&&  0+nn< jj		&& jj <NN.y-1-nn	&&  0+nn< kk		&& kk <NN.z-1-nn)
		{	
			eeFlux[4] +=  Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z;
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z/2;
		}
		if ( 0+nn< ii		&& ii <NN.x-1-nn	&&  0+nn< jj		&& jj <NN.y-1-nn	&&  0+nn< kk		&& kk <NN.z-1-nn)
			eeFlux[7] += rhoijk(ii,jj,kk,pphi,dd,NN)*dd.x*dd.y*dd.z;		
		/*
		 if(ii == 0)		eeFlux[0] += -Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z;
		 if(ii == NN.x-1)	eeFlux[0] +=  Eijk(ii,jj,kk,pphi,dd,NN)[1]*dd.y*dd.z;
		 if(jj == 0)		eeFlux[0] += -Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x;
		 if(jj == NN.y-1)	eeFlux[0] +=  Eijk(ii,jj,kk,pphi,dd,NN)[2]*dd.z*dd.x;
		 if(kk == 0)		eeFlux[0] += -Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y;
		 if(kk == NN.z-1)	eeFlux[0] +=  Eijk(ii,jj,kk,pphi,dd,NN)[3]*dd.x*dd.y;
		 */
	};
	eeFlux[1]*= PMC.eps0; eeFlux[2]*= PMC.eps0; eeFlux[3]*= PMC.eps0;
	eeFlux[4]*= PMC.eps0; eeFlux[5]*= PMC.eps0;	eeFlux[6]*= PMC.eps0;
	eeFlux[0] = eeFlux[1] + eeFlux[2] + eeFlux[3] + eeFlux[4] + eeFlux[5] + eeFlux[6];
	return eeFlux;
}
/**************************************************************************************/

/**************************************************************************************/
/* For Charge																		  */
/**************************************************************************************/

double rhoijk(int ii, int jj, int kk, const CMatrix3D& pphi, const StepsSizes& dd, const BoxSteps& NN)
{
	
	double d2x=0;
	double d2y=0;
	double d2z=0;
	if (ii == 0)			d2x = (pphi(ii+2, jj , kk )-2*pphi(ii+1, jj , kk )+pphi( ii , jj , kk ))/(2*pow(dd.x,2));
	if (0<ii && ii<NN.x-1) 	d2x = (pphi(ii+1, jj , kk )-2*pphi( ii , jj , kk )+pphi(ii-1, jj , kk ))/   pow(dd.x,2) ;
	if (ii == NN.x-1)		d2x = (pphi( ii , jj , kk )-2*pphi(ii-1, jj , kk )+pphi(ii-2, jj , kk ))/(2*pow(dd.x,2));

	if (jj == 0)			d2y = (pphi( ii ,jj+2, kk )-2*pphi( ii ,jj+1, kk )+pphi( ii , jj , kk ))/(2*pow(dd.y,2));
	if (0<jj && jj<NN.y-1) 	d2y = (pphi( ii ,jj+1, kk )-2*pphi( ii , jj , kk )+pphi( ii ,jj-1, kk ))/   pow(dd.y,2) ;
	if (jj == NN.y-1)		d2y = (pphi( ii , jj , kk )-2*pphi( ii ,jj-1, kk )+pphi( ii ,jj-2, kk ))/(2*pow(dd.y,2));

	if (kk == 0)			d2z = (pphi( ii , jj ,kk+2)-2*pphi( ii , jj ,kk+1)+pphi( ii , jj , kk ))/(2*pow(dd.z,2));
	if (0<kk && kk<NN.z-1) 	d2z = (pphi( ii , jj ,kk+1)-2*pphi( ii , jj , kk )+pphi( ii , jj ,kk-1))/   pow(dd.z,2) ;
	if (kk == NN.z-1)		d2z = (pphi( ii , jj , kk )-2*pphi( ii , jj ,kk-1)+pphi( ii , jj ,kk-2))/(2*pow(dd.z,2));

	if(ii!=0 && ii!=NN.x-1 && jj!=0 && jj!=NN.y-1 && kk!=0 && kk!=NN.z-1)
		return  -PMC.eps0*(d2x+d2y+d2z);
	else
		return 0.;
}

CMatrix3D Globalrho(const CMatrix3D& pphi, const StepsSizes& dd, const BoxSteps& NN)
{
	CMatrix3D	rrho(NN.x,NN.y,NN.z);
	for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
		rrho[ii][jj][kk]=rhoijk(ii,jj,kk,pphi,dd,NN);
	return rrho;
}

double ChannelCharge(const CMatrix3D& rrho, const CMatrix3D& UUn, const StepsSizes& dd, const BoxSteps& NN)
{
	double QQ=0;
	for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
//	for(int ii=1 ; ii<NN.x-1 ; ii++) for(int jj=1 ; jj<NN.y-1 ; jj++) for(int kk=1 ; kk<NN.z-1 ; kk++)
		if(UUn[ii][jj][kk] != 0)
			QQ += rrho[ii][jj][kk]*dd.x*dd.y*dd.z;
	return QQ;
};

double ChannelChargePositive(const CMatrix3D& rrho, const CMatrix3D& UUn, const StepsSizes& dd, const BoxSteps& NN)
{
	double QQ_plus=0;
	double rro;

	for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
//	for(int ii=1 ; ii<NN.x-1 ; ii++) for(int jj=1 ; jj<NN.y-1 ; jj++) for(int kk=1 ; kk<NN.z-1 ; kk++)
	{
		rro = rrho[ii][jj][kk];
		if(UUn[ii][jj][kk] != 0 && rro>=0) QQ_plus += rro*dd.x*dd.y*dd.z;
	}
	return QQ_plus;
};

double ChannelChargeNegative(const CMatrix3D& rrho, const CMatrix3D& UUn, const StepsSizes& dd, const BoxSteps& NN)
{
	double QQ_minus=0;
	double rro;

	for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)	
//	for(int ii=1 ; ii<NN.x-1 ; ii++) for(int jj=1 ; jj<NN.y-1 ; jj++) for(int kk=1 ; kk<NN.z-1 ; kk++)
	{
		rro = rrho[ii][jj][kk];
		if(UUn[ii][jj][kk] == 1 && rro<=0) QQ_minus += rro*dd.x*dd.y*dd.z;
	};
	return QQ_minus;
};

double TotalCharge(const CMatrix3D& rrho, const StepsSizes& dd, const BoxSteps& NN)
{
	double QQ=0;
	double ddV = dd.x*dd.y*dd.z;
	for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
//	for(int ii=1 ; ii<NN.x-1 ; ii++) for(int jj=1 ; jj<NN.y-1 ; jj++) for(int kk=1 ; kk<NN.z-1 ; kk++)
		QQ+= rrho[ii][jj][kk]*ddV;
	return QQ;
};

CMatrix1D ChannelLinearDensity(const CMatrix3D& rrho, const CMatrix3D& UUn, const StepsSizes& dd, const BoxSteps& NN)
{
	cout<<"*** Derivation valid for a single link in z-direction. ***\n";
	
	CMatrix1D	rrhol(NN.z);
	int			iic(0), jjc(0), kkcStart, kkcEnd;

	for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
//	for(int ii=1 ; ii<NN.x-1 ; ii++) for(int jj=1 ; jj<NN.y-1 ; jj++) for(int kk=1 ; kk<NN.z-1 ; kk++)
		if(UUn[ii][jj][kk] == 1)
		{
			iic=ii;
			jjc=jj;
			if(kk>0 && UUn[ii][jj][kk-1] == 0)	kkcStart=kk;
			if(kk<NN.z-1 && UUn[ii][jj][kk+1] == 0)	kkcEnd=kk;					
		};
		
	// solution 1 : integrate on a slice of the simulation box //
	/*	
	for(int kk=0 ; kk<NN.z ; kk++)
	{
		rrhol(kk) = 0;
		for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++)
				rrhol(kk) += rrho[ii][jj][kk]*dd.x*dd.y;
	}
	*/	
	// solution 2 : assume all charge is concentrated in the channel //
	for(int kk=0 ; kk<NN.z ; kk++) rrhol(kk) = rrho[iic][jjc][kk]*dd.x*dd.y;
	
	return rrhol;
}
/**************************************************************************************/

/**************************************************************************************/
/* For Dipole Moment																  */
/**************************************************************************************/
Vector DipoleMoment(double& CCarriedCharge, const CMatrix3D& pphi, const CMatrix3D& UUn, const BoxLengths& LL, const StepsSizes& dd, const BoxSteps& NN)
{
	Vector pp;
	double XXc(LL.x/2), YYc(LL.y/2), ZZc(0);	// define origine of the simulation domain
	double ddq		= 0;						// elementary charge
	pp.x			= 0;
	pp.y			= 0;
	pp.z			= 0;
	CCarriedCharge	= 0;

	for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
//	for(int ii=1 ; ii<NN.x-1 ; ii++) for(int jj=1 ; jj<NN.y-1 ; jj++) for(int kk=1 ; kk<NN.z-1 ; kk++)
		if(UUn[ii][jj][kk] == 1)
		{
			ddq = rhoijk(ii,jj,kk, pphi,dd,NN)*dd.x*dd.y*dd.z;
			pp.x += ddq*(ii*dd.x-XXc);
			pp.y += ddq*(jj*dd.y-YYc);
			pp.z += ddq*(kk*dd.z-ZZc);
			
			if(ddq >=0) CCarriedCharge += ddq;
		}
	return pp;
}
/**************************************************************************************/

/**************************************************************************************/
/* Misc.																			  */
/**************************************************************************************/
void SwitchValues(double& x1, double& x2)
{
	double tmp(x1);
	x1 = x2;
	x2 = tmp;
}
void SwitchValues(int& n1, int& n2)
{
	int tmp(n1);
	n1 = n2;
	n2 = tmp;
}

bool isfinite(const CMatrix3D& MM, const BoxSteps& NN)
{
	bool flag = true;
	for(int ii=0 ; ii<NN.x ; ii++) for(int jj=0 ; jj<NN.y ; jj++) for(int kk=0 ; kk<NN.z ; kk++)
		if (isnan(MM(ii,jj,kk)) || isinf(MM(ii,jj,kk)))
		{
			cout<<"!!! Error inproper value.\n Break at point: ["<<ii<<" "<<jj<<" "<<kk<<"].\n"<<endl;
			flag = false;
			break;
		};
	return flag;
}
/**************************************************************************************/