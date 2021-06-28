#include <stdlib.h>
#include <memory>
#include <iostream>
#include <string.h>
#include <cstdio>
#include "Matrix.h"
using namespace std;

/**************************************************************************************/
/* Virtual class of matrices, to be derived but no instance should be used.			  */
/**************************************************************************************/
const int CMatrix::DOUBLE_BYTE_SIZE = sizeof(double);

void CMatrix::freeMemIfRequired(void)
{
	if (pElems != NULL)	{free (pElems);}
	pElems = NULL;
} // Free Memory

bool CMatrix::init()
{
	freeMemIfRequired();
	int iMemSize = iNbElem* DOUBLE_BYTE_SIZE;
	if(!(pElems = (double*)malloc(iMemSize))) {return false;}
	memset(pElems, 0, iMemSize);
	return true;
} // Init

CMatrix& CMatrix::operator=(const CMatrix& matrix)
{
	iNbElem = matrix.iNbElem;
//	if (!init()) {cout<<"Not Enough Memory"<<endl;};
	if (!init()) {throw ENotEnoughMemory();};

	int iMemSize = iNbElem * DOUBLE_BYTE_SIZE;
	memcpy(pElems, matrix.pElems, iMemSize);
	
	return *this;
} // Operator =

CMatrix& CMatrix::operator+=(const CMatrix& matrix)
{
//	if (iNbElem != matrix.iNbElem) {cout<<"Impossible Addition"<<endl;}
	if (iNbElem != matrix.iNbElem) {throw EAdditionImpossible();}
	for (int iiElem = 0; iiElem < iNbElem; iiElem++) {pElems[iiElem] += matrix.pElems[iiElem];}
	return *this;
} // Operator+=

CMatrix& CMatrix::operator-=(const CMatrix& matrix)
{
	//	if (iNbElem != matrix.iNbElem) {cout<<"Impossible Addition"<<endl;}
	if (iNbElem != matrix.iNbElem) {throw EAdditionImpossible();}
	for (int iiElem = 0; iiElem < iNbElem; iiElem++) {pElems[iiElem] -= matrix.pElems[iiElem];}
	return *this;
} // Operator-=

CMatrix::CMatrix (int iiNbElem)
{
	iNbElem = iiNbElem;
	pElems  = NULL;
//	if (!init()) {cout<<"Not Enough Memory"<<endl;}
	if (!init()) {throw ENotEnoughMemory();}
} // Default constructor

CMatrix::CMatrix(const CMatrix& matrix) 
{
	pElems  = NULL;
	iNbElem = matrix.iNbElem;
	(*this) = matrix;
} // Constructor by copy

int CMatrix::getNbElem()
{
	return iNbElem;
}

double CMatrix::getElem(int ii)
{
	return pElems[ii];
}

bool CMatrix::operator==(const CMatrix& matrix) const
{
	if (iNbElem != matrix.iNbElem)
		return false;
	
	if (memcmp(pElems, matrix.pElems, iNbElem) != 0)
		return false;
	
	return true;
}  // Operator ==

FILE* CMatrix::Open(char *fname,char *flags)
{
	FILE *file;
	char *relativePath;
#ifdef __MWERKS__
	/* Parse the file name for the Macintosh */
	int hierarchy;
	char MacPath[256],*Mname,*Uname;
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
	if(!(file=fopen(relativePath,flags)))
		cout<<"cannot open file";
	//nrerror("cannot open file"); 
	return file;
}

void CMatrix::MinMax(double& MMin, double& MMax)
{
	if ( iNbElem <= 0 ) 
		throw; // Pb, max sur une matrice vide
	else
	{
		MMin=pElems[0];
		MMax=pElems[0];
		for(int ii=0; ii<iNbElem; ii++)
		{
			if(pElems[ii]<=MMin) {MMin=pElems[ii];};
			if(pElems[ii]>=MMax) {MMax=pElems[ii];};
		}
	}
}
/**************************************************************************************/

/**************************************************************************************/
/* A 3 dimension matrix containing doubles.											  */
/**************************************************************************************/
void CMatrix3D::freeMemIfRequired(void)
{
	CMatrix::freeMemIfRequired();
	if (pMatrix3d != NULL)
	{
		free (pMatrix2d);
		free (pMatrix3d);
	}
	pMatrix3d = NULL;	
} // freeMemIfRequired

bool CMatrix3D::init()
{
        if (!CMatrix::init())
		return false;
	
	if(!(pMatrix3d = (double***)malloc(xSize * sizeof(double**)))) 
	{
		CMatrix::freeMemIfRequired();
		return false;
	}
	if(!(pMatrix2d = (double**)malloc(xSize * ySize * sizeof(double* )))) 
	{
		CMatrix::freeMemIfRequired();
		free (pMatrix3d);
		pMatrix3d = NULL;
		return false;
	}
	
	for(int x=0 ; x<xSize ; x++)
	{
		pMatrix3d[x] = pMatrix2d + x*ySize;
		for(int y=0 ; y < ySize ; y++)
		{
			pMatrix3d[x][y] = pElems + zSize * ( y + x * ySize );
		}
	}
	return true;
} // Init

bool CMatrix3D::init(int xxSize, int yySize, int zzSize)
{
	iNbElem = xxSize*yySize*zzSize;
	xSize = xxSize;
	ySize = yySize;
	zSize = zzSize;
	return init();
} // bool CMatrix3D::init

CMatrix3D::CMatrix3D(int xxSize, int yySize, int zzSize) : CMatrix(xxSize * yySize * zzSize)
{
	xSize = xxSize;
	ySize = yySize;
	zSize = zzSize;
        pMatrix2d = NULL;
	pMatrix3d = NULL;
//	if (!init()) {cout<<"Not Enough Memory";}
	if (!init()) {throw ENotEnoughMemory();}
} // Default constructor

CMatrix3D::CMatrix3D(const CMatrix3D& matrix) : CMatrix(matrix.xSize * matrix.ySize * matrix.zSize)
{
	xSize	  = matrix.xSize;
	ySize	  = matrix.ySize;
	zSize     = matrix.zSize;
	pMatrix3d = NULL;
	*this = matrix;
}// Copy constructor

CMatrix3D& CMatrix3D::operator=(const CMatrix3D& matrix)
{
	xSize = matrix.xSize;
	ySize = matrix.ySize;
	zSize = matrix.zSize;
	CMatrix::operator=(matrix);
	return *this;
} // operator=

CMatrix3D& CMatrix3D::operator+=(const CMatrix3D& matrix)
{
	if (xSize != matrix.xSize 
		|| ySize != matrix.ySize 
		|| zSize != matrix.zSize)
		throw EAdditionImpossible();
//	{cout<<"Impossible Addition";};
	
	CMatrix::operator +=(matrix);
	return *this;
} // CMatrix3D& CMatrix3D::operator+=

CMatrix3D& CMatrix3D::operator-=(const CMatrix3D& matrix)
{
	if (xSize != matrix.xSize 
		|| ySize != matrix.ySize 
		|| zSize != matrix.zSize)
		throw EAdditionImpossible();
	//	{cout<<"Impossible Addition";};
	
	CMatrix::operator -=(matrix);
	return *this;
} // CMatrix3D& CMatrix3D::operator-=

CMatrix3D CMatrix3D::operator+(const CMatrix3D& matrix) const
{
	CMatrix3D result(*this);
	result += matrix;
	return result;
} // CMatrix3D& CMatrix3D::operator+

CMatrix3D CMatrix3D::operator-(const CMatrix3D& matrix) const
{
	CMatrix3D result(*this);
	result -= matrix;
	return result;
} // CMatrix3D& CMatrix3D::operator-

ostream& operator<<(ostream& os, const CMatrix3D& matrix)
{
	for (int z = 0; z < matrix.zSize; z++)
	{
		os << "z=" << z << endl;
		for (int x = 0; x < matrix.xSize; x++)
		{
			for (int y = 0; y < matrix.ySize; y++)
			{
				os <<setw(sizeof(double)+1)  << matrix[x][y][z] << '\t';
			}
			os << endl;
		}
		os << endl;
	}
	return os;
} // Operator<<

void CMatrix3D::fwrite(string s)
{
	char * c = &s[0];
	FILE * ffile = Open (c, "w");
	for (int kk=0 ; kk<zSize ; kk++)
	{
//		fprintf(ffile,"z = %d\n",kk);
		for (int ii=0 ; ii<xSize ; ii++)
		{
			for (int jj=0 ; jj<ySize ; jj++) fprintf(ffile,"%f ", pMatrix3d[ii][jj][kk]);
			fprintf(ffile,"\n");
		};
		fprintf(ffile,"\n");
	};
	fclose(ffile);
}// fwrite

/**************************************************************************************/

/**************************************************************************************/
/* A 2 dimension matrix containing doubles.											  */
/**************************************************************************************/
void CMatrix2D::freeMemIfRequired(void)
{
	if (pMatrix2d != NULL)
	{
		free (pMatrix2d);
	}
	pMatrix2d = NULL;
	CMatrix::freeMemIfRequired();
} // CMatrix2D::freeMem

bool CMatrix2D::init()
{
	if (!CMatrix::init()) {return false;}
	
	if(!(pMatrix2d = (double**)malloc(xSize * sizeof(double*)))) 
	{
		CMatrix::freeMemIfRequired();
		return false;
	}
	
	for(int x=0 ; x<xSize ; x++) {pMatrix2d[x] = pElems + x*ySize;}
	
	return true;
} // bool CMatrix2D::init

bool CMatrix2D::init(int xxSize, int yySize)
{
	iNbElem = xxSize*yySize;
	xSize = xxSize;
	ySize = yySize;
	return init();
} // bool CMatrix2D::init

CMatrix2D::CMatrix2D(int xxSize, int yySize) : CMatrix(xxSize * yySize)
{
	xSize	  = xxSize; 
	ySize	  = yySize;
	pMatrix2d = NULL;
//	if (!init()) {cout<<"Not Enough Memory"<<endl;}
	if (!init()) {throw ENotEnoughMemory();}
} // CMatrix2D::CMatrix2D

CMatrix2D::CMatrix2D(const CMatrix2D& matrix)
{
	xSize     = matrix.xSize;
	ySize	  = matrix.ySize;
	pMatrix2d = NULL;
	*this = matrix;
}// CMatrix2D::CMatrix2D

CMatrix2D& CMatrix2D::operator=(const CMatrix2D& matrix)
{
	xSize = matrix.xSize;
	ySize = matrix.ySize;
	CMatrix::operator=(matrix);
	
	return *this;
} // CMatrix2D& CMatrix2D::operator=

CMatrix2D CMatrix2D::operator+(const CMatrix2D& matrix) const
{
	CMatrix2D result(*this);
	result += matrix;
	return result;
} // CMatrix2D& CMatrix2D::operator+

CMatrix2D CMatrix2D::operator-(const CMatrix2D& matrix) const
{
	CMatrix2D result(*this);
	result -= matrix;
	return result;
} // CMatrix2D& CMatrix2D::operator-

CMatrix2D& CMatrix2D::operator+=(const CMatrix2D& matrix)
{
//	if (xSize != matrix.xSize || ySize != matrix.ySize ) {cout<<"Impossible Addition"<<endl;}
	if (xSize != matrix.xSize || ySize != matrix.ySize ) {throw EAdditionImpossible();}
	
	CMatrix::operator +=(matrix);
	return *this;
} // CMatrix2D& CMatrix2D::operator+=

CMatrix2D& CMatrix2D::operator-=(const CMatrix2D& matrix)
{
	//	if (xSize != matrix.xSize || ySize != matrix.ySize ) {cout<<"Impossible Addition"<<endl;}
	if (xSize != matrix.xSize || ySize != matrix.ySize ) {throw EAdditionImpossible();}
	
	CMatrix::operator -=(matrix);
	return *this;
} // CMatrix2D& CMatrix2D::operator-=

ostream& operator<<(ostream& os, const CMatrix2D& matrix)
{
	for (int x = 0; x < matrix.xSize; x++)
	{
		for (int y = 0; y < matrix.ySize; y++)
		{
			os << matrix(x, y) << '\t';
		}
		os << endl;
	}
	os << endl;
	return os;
} // std::ostream& operator<<

void CMatrix2D::fwrite(string s)
{
	char * c = &s[0];
	FILE * ffile = Open (c, "w");
	for (int jj=0 ; jj<ySize ; jj++)
	{
		for (int ii=0 ; ii<xSize ; ii++) fprintf(ffile,"%f ", pMatrix2d[ii][jj]);
		fprintf(ffile,"\n");
	};
	fclose(ffile);
}// fwrite
/**************************************************************************************/

/**************************************************************************************/
/* A 1 dimension matrix containing doubles.											  */
/**************************************************************************************/
bool CMatrix1D::init(int xSize)
{
	iNbElem = xSize;
	return CMatrix::init();
} // bool CMatrix1D::init

CMatrix1D& CMatrix1D::operator=(const CMatrix1D& matrix)
{
	CMatrix::operator=(matrix);
	return *this;
} // CMatrix1D& CMatrix1D::operator=

CMatrix1D CMatrix1D::operator+(const CMatrix1D& matrix) const
{
	CMatrix1D result(*this);
	result += matrix;
	return result;
} // CMatrix1D& CMatrix1D::operator+

CMatrix1D CMatrix1D::operator-(const CMatrix1D& matrix) const
{
	CMatrix1D result(*this);
	result -= matrix;
	return result;
} // CMatrix1D& CMatrix1D::operator-

CMatrix1D& CMatrix1D::operator+=(const CMatrix1D& matrix)
{
	CMatrix::operator +=(matrix);
	return *this;
} // CMatrix1D& CMatrix1D::operator+=

CMatrix1D& CMatrix1D::operator-=(const CMatrix1D& matrix)
{
	CMatrix::operator -=(matrix);
	return *this;
} // CMatrix1D& CMatrix1D::operator-=

ostream& operator<<(ostream& os, const CMatrix1D& matrix)
{
	for (int x = 0; x < matrix.iNbElem; x++)
	{
		os << matrix(x) << '\t';
	}
	os << endl;
	return os;
} // std::ostream& operator<<

void CMatrix1D::fwrite(string s)
{
	char * c = &s[0];
	FILE * ffile = Open (c, "w");
	for (int ii=0 ; ii<iNbElem ; ii++)
		fprintf(ffile,"%f\n",pElems[ii]);
	fclose(ffile);
}// fwrite
/**************************************************************************************/