/* @file matrix.h this file presents a base class for matrixes and the operations available on it. */

#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <iomanip>
//#include <memory>
//#include <exception>
using namespace std;

/**************************************************************************************/
/* Virtual class of matrices, to be derived but no instance should be used.			  */
/**************************************************************************************/
class CMatrix
{
protected:
	static const int DOUBLE_BYTE_SIZE;				// Size in bytes of a double.
	virtual void freeMemIfRequired();				// Free the memory occupied by _pElems if necessary.
	virtual bool init();							// Initialize the dimension of the matrix and its content (all values are set to 0) according to the value of <i>_iNbItem</i>.
													// @return true if the memory could be allocated, false in the other case.
	CMatrix& operator=(const CMatrix& matrix);		// Overloading affectation operator
													// @throw ENotEnoughMemory if there is not enough memory to complete the copy.
	CMatrix& operator+=(const CMatrix& matrix);		// Overloading self addition operator.
													// throw EAdditionImpossible is the addition is impossible
	CMatrix& operator-=(const CMatrix& matrix);		// Overloading self difference operator.
													// throw EAdditionImpossible is the addition is impossible
	int iNbElem;									// The number of elements in <i>_pElems</i>.
	double* pElems;									// An array of <i>iNbElem</i> doubles which are the content of the matrix.

public:
	class EMatrix : public exception{};				// Base class for exceptions raised by methods of the CMatrix class.
	class ENotEnoughMemory : public EMatrix{};		// Exception raised if memory is missing to complete an operation.
	class EAdditionImpossible : public EMatrix{};	// Exception raised if the addition of 2 matrix is impossible.

	CMatrix(int iNbElem=0);							// Default constructor.
													// @throw ENotEnoughMemory if there is not enough memory to complete the creation.
	CMatrix( const CMatrix& matrix);				// constructor by copy.
													// @throw ENotEnoughMemory if there is not enough memory to complete the creation.
	virtual ~CMatrix(void) {freeMemIfRequired();};	// Destructor.

	int getNbElem();								// Get the number of elements in the matrix
	double getElem(int ii);							// Get the ii-th elements of the matrix
	bool operator==(const CMatrix& matrix) const;	// Compare two matrices, return TRUE if identical
	FILE* Open(char *fname,char *flags);			// Open file for writing
	void MinMax(double& MMin, double& MMax);		// Return min and max values of the matrix
};

/**************************************************************************************/

/**************************************************************************************/
/* A 3 dimension matrix containing doubles.											  */
/**************************************************************************************/
class CMatrix3D : public CMatrix
{

protected:
	int xSize;										// Width of the matrix.
	int ySize;										// Height of the matrix.
	int zSize;										// The depth of the matrix.

	double ** pMatrix2d;
	virtual void freeMemIfRequired(void);			// Free the memory occupied by _pMatrix3d if necessary.
	virtual bool init();							// Initialize the dimension of the matrix and its content (all values are set to 0).
													// @return true if the memory could be allocated, false in the other case.
public:
	double *** pMatrix3d;							// A shortcul to access the elements stored in the
													// <i>pElems</i> variable
													//		* inherited from the base class.
													//		* The element [x, y, z] is available at the adress : _pMatrix3d[x][y][z]
	bool init(int xxSize, int yySize, int zzSize);	// Dimensionnate the matrix
	CMatrix3D( int xSize = 0,  int ySize = 0,  int zSize = 0);
													// Default constructor.
													// @throw ENotEnoughMemory if there is not enough memory to complete the creation.
	CMatrix3D(const CMatrix3D& matrix);				// Copy constructor.
													// @throw ENotEnoughMemory if there is not enough memory to complete the creation.
	virtual ~CMatrix3D(void) {freeMemIfRequired();};

	inline double** operator[](int x) __attribute__((always_inline)) {return pMatrix3d[x];};
													// Overloading the bracket operator so that an element of the matrix can be accessed using [x][y][z]
	inline double** operator[](int x) const __attribute__((always_inline)) {return pMatrix3d[x];};
													// Overloading the bracket operator so that an element of the matrix can be accessed using [x][y][z]
	inline double& operator()(int x, int y, int z) __attribute__((always_inline)) {return pMatrix3d[x][y][z];};
													// Overloading the parenthesis operator so that an element of the matrix can be accessed using (x, y, z).
	inline const double& operator()(int x, int y, int z) const __attribute__((always_inline)) {return pMatrix3d[x][y][z];};
													// Overloading the parenthesis operator so that an element of the matrix can be accessed using (x, y, z).

	CMatrix3D& operator=(const CMatrix3D& matrix);	// Overloading affectation operator
													// @throw ENotEnoughMemory if there is not enough memory to complete the copy.
	CMatrix3D& operator+=(const CMatrix3D& matrix);	// Overloading self addition operator.
													// @throw EAdditionImpossible is the addition is impossible.
	CMatrix3D& operator-=(const CMatrix3D& matrix);	// Overloading self difference operator.
													// @throw EAdditionImpossible is the addition is impossible.
	CMatrix3D operator+(const CMatrix3D& matrix) const;
													// Overloading addition operator.
													//	* @throw EAdditionImpossible is the addition is impossible.
													//  * @throw ENotEnoughMemory if there is not enough memory to complete the creation.
	CMatrix3D operator-(const CMatrix3D& matrix) const;
													// Overloading difference operator.
													//	* @throw EAdditionImpossible is the addition is impossible.
													//  * @throw ENotEnoughMemory if there is not enough memory to complete the creation.
	friend ostream& operator<<(ostream& os, const CMatrix3D& matrix);
													// Overloading display operator.
	void fwrite(char *);							// Write Matrix in a file
}; // class CMatrix3D
/**************************************************************************************/

/**************************************************************************************/
/* A 2 dimension matrix containing doubles.											  */
/**************************************************************************************/
class CMatrix2D : public CMatrix
{
protected:
	int xSize;										// Width of the matrix.
	int ySize;										// Height of the matrix.
	double ** pMatrix2d;							// A shortcul to access the elements stored in the <i>pElems</i> variable
													// * inherited from the base class.
													// * The element [x, y] is available at the adress: pMatrix2d[x][y]
	virtual void freeMemIfRequired(void);			// Free the memory occupied by _pMatrix2d if necessary.
	virtual bool init();							// Initialize the dimension of the matrix and its content (all values are set to 0).
													// @return true if the memory could be allocated, false in the other case.

public:
	bool init(int xSize, int ySize);				// Initialize the matrix.
	CMatrix2D( int xSize = 0,  int ySize = 0);		// Default constructor.
													// @throw ENotEnoughMemory if there is not enough memory to complete the creation.
	CMatrix2D(const CMatrix2D& matrix);				// Copy constructor.
													// @throw ENotEnoughMemory if there is not enough memory to complete the creation.
	virtual ~CMatrix2D(void) {freeMemIfRequired();};// Destructor

	double* operator[](int x) __attribute__((always_inline)) {return pMatrix2d[x];};
													// Overloading the bracket operator so that an element of the matrix can be accessed using [x][y]
	double* operator[](int x) const __attribute__((always_inline)) {return pMatrix2d[x];};
													// Overloading the parenthesis operator so that an element of the matrix can be accessed using (x, y).
	double& operator()(int x, int y) __attribute__((always_inline))	{return pMatrix2d[x][y];};
													// Overloading the parenthesis operator so that an element of the matrix can be accessed using (x, y).
	const double& operator()(int x, int y) const __attribute__((always_inline))	{return pMatrix2d[x][y];};
													// Overloading the parenthesis operator so that an element of the matrix can be accessed using (x, y).

	CMatrix2D& operator=(const CMatrix2D& matrix);	// Overloading affectation operator
													// @throw ENotEnoughMemory if there is not enough memory to complete the copy.
	CMatrix2D operator+(const CMatrix2D& matrix) const;
													// Overloading addition operator.
													// * @throw EAdditionImpossible is the addition is impossible.
													// * @throw ENotEnoughMemory if there is not enough	memory to complete the creation.
	CMatrix2D operator-(const CMatrix2D& matrix) const;
													// Overloading difference operator.
													// * @throw EAdditionImpossible is the addition is impossible.
													// * @throw ENotEnoughMemory if there is not enough	memory to complete the creation.
	CMatrix2D& operator+=(const CMatrix2D& matrix); // Overloading self addition operator.
													// @throw EAdditionImpossible is the addition is impossible.
	CMatrix2D& operator-=(const CMatrix2D& matrix); // Overloading self difference operator.
													// @throw EAdditionImpossible is the addition is impossible.
	friend ostream& operator<<(ostream& os, const CMatrix2D& matrix);
													// Overloading display operator.
	void fwrite(char *);							// Write Matrix in a file
}; // class CMatrix2D
/**************************************************************************************/

/**************************************************************************************/
/* A 1 dimension matrix containing doubles.											  */
/**************************************************************************************/
class CMatrix1D : public CMatrix
{
public:
	CMatrix1D(int xSize = 0) : CMatrix(xSize) {};	// Default constructor

	CMatrix1D(const CMatrix1D& matrix) : CMatrix(matrix) {};
													// Copy constructor

	virtual ~CMatrix1D(void) {freeMemIfRequired();} // CMatrix::~CMatrix

	bool init(int xSize);							// Initialise la matrice

	double& operator[](int x) __attribute__((always_inline)) {return pElems[x];};
													// Overloading the bracket operator so that an element of the matrix can be accessed using [x]
	double& operator[](int x) const __attribute__((always_inline)) {return pElems[x];};
													// Overloading the bracket operator so that an element of the matrix can be accessed using [x]
	double& operator()(int x) __attribute__((always_inline)) {return pElems[x];};
													// Overloading the parenthesis operator so that an element of the matrix can be accessed using (x).
	const double& operator()(int x) const __attribute__((always_inline)) {return pElems[x];};
													// Overloading the parenthesis operator so that an element of the matrix can be accessed using (x).

	CMatrix1D& operator=(const CMatrix1D& matrix);	// Overloading affectation operator
													// @throw ENotEnoughMemory if there is not enough memory to complete the copy.

	CMatrix1D operator+(const CMatrix1D& matrix) const;
													// Overloading addition operator.
													// * @throw EAdditionImpossible is the addition is impossible.
													// * @throw ENotEnoughMemory if there is not enough memory to complete the creation.
	CMatrix1D operator-(const CMatrix1D& matrix) const;
													// Overloading difference operator.
													// * @throw EAdditionImpossible is the addition is impossible.
													// * @throw ENotEnoughMemory if there is not enough memory to complete the creation.

	CMatrix1D& operator+=(const CMatrix1D& matrix); // Overloading self addition operator.
													// @throw EAdditionImpossible is the addition is impossible.
	CMatrix1D& operator-=(const CMatrix1D& matrix); // Overloading self difference operator.
													// @throw EAdditionImpossible is the addition is impossible.

	friend std::ostream& operator<<(std::ostream& os, const CMatrix1D& matrix);
													// Overloading display operator.
	void fwrite(char *);							// Write Matrix in a file
}; // class CMatrix1D
/**************************************************************************************/
#endif // MATRIX_H
