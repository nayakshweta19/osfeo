/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.20 $
// $Date: 2009/12/23 22:54:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/matrix/Vector.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class implementation for Vector.
//
// What: "@(#) Vector.C, revA"

#include "Vector.h"
#include "Matrix.h"
#include "ID.h"
#include <Tensor.h>
#include <iostream>
using std::nothrow;

#include <stdlib.h>
#include <math.h>

double Vector::VECTOR_NOT_VALID_ENTRY =0.0;

// Vector():
//	Standard constructor, sets size = 0;

Vector::Vector()
: sz(0), theData(0), fromFree(0)
{

}

// Vector(int size):
//	Constructor used to allocate a vector of size size.

Vector::Vector(int size)
: sz(size), theData(0), fromFree(0)
{
#ifdef _G3DEBUG
  if (sz < 0) {
    opserr << "Vector::Vector(int) - size " << size << " specified < 0\n";
    sz = 1;
  }
#endif

  // get some space for the vector
  //  theData = (double *)malloc(size*sizeof(double));
  if (size > 0) {
    theData = new (nothrow) double [size];

    if (theData == 0) {
      opserr << "Vector::Vector(int) - out of memory creating vector of size " << size << endln;
      sz = 0; // set this should fatal error handler not kill process!!
    }
    
    // zero the components
    for (int i=0; i<sz; i++)
      theData[i] = 0.0;
  }
}


// Vector::Vector(double *data, int size)

Vector::Vector(double *data, int size)
: sz(size),theData(data),fromFree(1)
{
#ifdef _G3DEBUG
  if (sz <= 0) {
    opserr << "Vector::Vector(double *, size) - size " << size << " specified <= 0\n";
    sz = 0;
  }
#endif
}
 


// Vector(const Vector&):
//	Constructor to init a vector from another.

Vector::Vector(const Vector &other)
: sz(other.sz),theData(0),fromFree(0)
{
#ifdef _G3DEBUG
  if (sz < 0) {
    opserr << "Vector::Vector(int) - size " << sz << " specified <= 0\n";
    sz = 1;
  }
#endif

  //  theData = (double *)malloc(other.sz*sizeof(double));    
  theData = new (nothrow) double [other.sz];    
  
  if (theData == 0) {
    opserr << "Vector::Vector(int) - out of memory creating vector of size " << sz << endln;
    exit(-1);
  }

  // copy the component data
  for (int i=0; i<sz; i++)
    theData[i] = other.theData[i];
}	


// ~Vector():
// 	destructor, deletes the [] data

Vector::~Vector()
{
  if (sz != 0 && fromFree == 0) 
    delete [] theData;
  //  free((void *)theData);
}


int 
Vector::setData(double *newData, int size){
  if (sz != 0 && fromFree == 0) 
    delete [] theData;      
  sz = size;
  theData = newData;
  fromFree = 1;

  if (sz <= 0) {
    opserr << " Vector::Vector(double *, size) - size specified: " << size << " <= 0\n";
    sz = 0;
  }

  return 0;
}



int 
Vector::resize(int newSize){

  // first check that newSize is valid
  if (newSize < 0) {
    opserr << "Vector::resize) - size specified " << newSize << " <= 0\n";
    return -1;
  } 
  
  // otherwise if newSize is gretaer than oldSize free old space and get new space
  else if (newSize > sz) {

    // delete the old array
    if (sz != 0 && fromFree == 0) 
	delete [] theData;
    sz = 0;
    fromFree = 0;
    
    // create new memory
    // theData = (double *)malloc(newSize*sizeof(double));    
    theData = new (nothrow) double[newSize];
    if (theData == 0) {
      opserr << "Vector::resize() - out of memory for size " << newSize << endln;
      sz = 0;
      return -2;
    }
    sz = newSize;
  }  

  // just set the size to be newSize .. penalty of holding onto additional
  // memory .. but then save the free() and malloc() calls
  else 
    sz = newSize;    

  return 0;
}


// Assemble(Vector &x, ID &y, double fact ):
//	Method to assemble into object the Vector V using the ID l.
//	If ID(x) does not exist program writes error message if
//	VECTOR_CHECK defined, otherwise ignores it and goes on.

int 
Vector::Assemble(const Vector &V, const ID &l, double fact )
{
  int result = 0;
  int pos;
  for (int i=0; i<l.Size(); i++) {
    pos = l(i);
    
    if (pos < 0)
      ;
    else if ((pos < sz) && (i < V.Size()))
      // assemble into vector
      theData[pos] += V.theData[i] *fact;
    else {
      result = -1;
      if (pos < sz)
	opserr << "Vector::Assemble() " << pos << " out of range [1, " << sz-1 << "]\n";
      else
	opserr << "Vector::Assemble() " << pos << " out of range [1, "<< V.Size()-1 << "]\n";
    }
  }
  return result;
}
    

int
Vector::Normalize(void) 
{
  double length = 0.0;
  for (int i=0; i<sz; i++)
    length += theData[i] * theData[i];
  length = sqrt(length);
  
  if (length == 0.0) 
    return -1;

  length = 1.0/length;
  for (int j=0; j<sz; j++)
    theData[j] *= length;

  return 0;
}

int
Vector::addVector(double thisFact, const Vector &other, double otherFact )
{
  // check if quick return
  if (otherFact == 0.0 && thisFact == 1.0)
    return 0; 

  // if sizes are compatable add
#ifdef _G3DEBUG
  if (sz != other.sz) {
    // else sizes are incompatable, do nothing but warning
    opserr <<  "WARNING Vector::addVector() - incompatable Vector sizes\n";
    return -1;
  }
#endif


  if (thisFact == 1.0) {

    // want: this += other * otherFact
    double *dataPtr = theData;
    double *otherDataPtr = other.theData;
    if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) 
	*dataPtr++ += *otherDataPtr++;
    } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) 
	*dataPtr++ -= *otherDataPtr++;
    } else 
      for (int i=0; i<sz; i++) 
	*dataPtr++ += *otherDataPtr++ * otherFact;
  } 

  else if (thisFact == 0.0) {

    // want: this = other * otherFact
    double *dataPtr = theData;
    double *otherDataPtr = other.theData;
    if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) 
	*dataPtr++ = *otherDataPtr++;
    } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) 
	*dataPtr++ = -(*otherDataPtr++);
    } else 
      for (int i=0; i<sz; i++) 
	*dataPtr++ = *otherDataPtr++ * otherFact;
  }

  else {

    // want: this = this * thisFact + other * otherFact
    double *dataPtr = theData;
    double *otherDataPtr = other.theData;
    if (otherFact == 1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) {
	double value = *dataPtr * thisFact + *otherDataPtr++;
	*dataPtr++ = value;
      }
    } else if (otherFact == -1.0) { // no point doing a multiplication if otherFact == 1.0
      for (int i=0; i<sz; i++) {
	double value = *dataPtr * thisFact - *otherDataPtr++;
	*dataPtr++ = value;
      }
    } else 
      for (int i=0; i<sz; i++) {
	double value = *dataPtr * thisFact + *otherDataPtr++ * otherFact;
	*dataPtr++ = value;
      }
  } 

  // successfull
  return 0;
}
	    
	
int
Vector::addMatrixVector(double thisFact, const Matrix &m, const Vector &v, double otherFact )
{
  // see if quick return
  if (thisFact == 1.0 && otherFact == 0.0)
    return 0;

  // check the sizes are compatable
#ifdef _G3DEBUG
  // check the sizes are compatable
  if ((sz != m.noRows()) && (m.noCols() != v.sz)) {
    // otherwise incompatable sizes
    opserr << "Vector::addMatrixVector() - incompatable sizes\n";
    return -1;    
  }
#endif

  if (thisFact == 1.0) {

    // want: this += m * v * otherFact
    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++;
	for (int j=0; j<sz; j++)
	  theData[j] += *matrixDataPtr++ * otherData;
      }
    } 
    else if (otherFact == -1.0) { // no point doing multiplication if otherFact = -1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++;
	for (int j=0; j<sz; j++)
	  theData[j] -= *matrixDataPtr++ * otherData;
      }
    } 
    else { // have to do the multiplication
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++ * otherFact;
	for (int j=0; j<sz; j++)
	  theData[j] += *matrixDataPtr++ * otherData;
      }
    }
  }

  else if (thisFact == 0.0) {
    
    // want: this = m * v * otherFact
    for (int i=0; i<sz; i++)
      theData[i] = 0.0;

    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++;
	for (int j=0; j<sz; j++)
	  theData[j] += *matrixDataPtr++ * otherData;
      }
    } 
    else if (otherFact == -1.0) { // no point doing multiplication if otherFact = -1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++;
	for (int j=0; j<sz; j++)
	  theData[j] -= *matrixDataPtr++ * otherData;
      }
    } else {
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++ * otherFact;
	for (int j=0; j<sz; j++)
	  theData[j] += *matrixDataPtr++ * otherData;
      }
    }
  }

  else {

    // want: this = this * thisFact + m * v * otherFact
    for (int i=0; i<sz; i++)
      theData[i] *= thisFact;

    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++;
	for (int j=0; j<sz; j++)
	  theData[j] += *matrixDataPtr++ * otherData;
      }
    } else if (otherFact == -1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++;
	for (int j=0; j<sz; j++)
	  theData[j] -= *matrixDataPtr++ * otherData;
      }
    } else {
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtr = v.theData;
      for (int i=0; i<otherSize; i++) {
	double otherData = *otherDataPtr++ * otherFact;
	for (int j=0; j<sz; j++)
	  theData[j] += *matrixDataPtr++ * otherData;
      }
    }
  }
  
  // successfull
  return 0;
}



int
Vector::addMatrixTransposeVector(double thisFact, 
				 const Matrix &m, 
				 const Vector &v, 
				 double otherFact )
{
  // see if quick return
  if (otherFact == 0.0 && thisFact == 1.0)
    return 0;

#ifdef _G3DEBUG
  // check the sizes are compatable
  if ((sz != m.noRows()) && (m.noRows() != v.sz)) {
    // otherwise incompatable sizes
    opserr << "Vector::addMatrixTransposeVector() - incompatable sizes\n";
    return -1;    
  }
#endif

  if (thisFact == 1.0) {

    // want: this += m^t * v * otherFact
    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	theData[i] += sum;
      }
    } else if (otherFact == -1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	theData[i] -= sum;
      }
    } else {
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	theData[i] += sum * otherFact;
      }
    }
  }

  else if (thisFact == 0.0) {

    // want: this = m^t * v * otherFact
    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	theData[i] = sum;
      }
    } else if (otherFact == -1.0) { // no point doing multiplication if otherFact = -1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	theData[i] = -sum;
      }
    } else {
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	theData[i] = sum * otherFact;
      }
    }
  } 

  else {

    // want: this = this * thisFact + m^t * v * otherFact
    if (otherFact == 1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	double value = theData[i] * thisFact + sum;
	theData[i] = value;
      }
    } else if (otherFact == -1.0) { // no point doing multiplication if otherFact = 1.0
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	double value = theData[i] * thisFact - sum;
	theData[i] = value;
      }
    } else {
      int otherSize = v.sz;
      double *matrixDataPtr = m.data;
      double *otherDataPtrA = v.theData;
      for (int i=0; i<sz; i++) {
	double *otherDataPtr = otherDataPtrA;
	double sum = 0.0;
	for (int j=0; j<otherSize; j++)
	  sum += *matrixDataPtr++ * *otherDataPtr++;
	double value = theData[i] * thisFact + sum * otherFact;
	theData[i] = value;
      }
    }
}

  return 0;
}
	
	



// double Norm();
//	Method to return the norm of  vector. (non-const as may save norm for later)

double
Vector::Norm(void) const
{
  double value = 0;
  for (int i=0; i<sz; i++) {
    double data = theData[i];
    value += data*data;
  }
  return sqrt(value);
}

double
Vector::pNorm(int p) const
{
  double value = 0;
  
  if (p>0) {
    for (int i=0; i<sz; i++) {
      double data = fabs(theData[i]);
      value += pow(data,p);
    }
    return pow(value,1.0/p);
  } else {
    for (int i=0; i<sz; i++) {
      double data = fabs(theData[i]);
      value = (data>value) ? data : value;
    }
    return value;
  }
}


double &
Vector::operator[](int x) 
{
#ifdef _G3DEBUG
  // check if it is inside range [0,sz-1]
  if (x < 0) {
    opserr << "Vector::() - x " << x << " outside range 0 - << " << sz-1 << endln;
    return VECTOR_NOT_VALID_ENTRY;
  }
#endif
  
  if (x >= sz) {
    double *dataNew = new (nothrow) double[x+1];
    for (int i=0; i<sz; i++)
      dataNew[i] = theData[i];
    for (int j=sz; j<x; j++)
      dataNew[j] = 0.0;
    
    if (fromFree == 1)
      if (theData != 0)
	delete [] theData;

    theData = dataNew;
    sz = x+1;
  }

  return theData[x];
}

double Vector::operator[](int x) const
{
#ifdef _G3DEBUG
  // check if it is inside range [0,sz-1]
  if (x < 0 || x >= sz) {
    opserr << "Vector::() - x " << x << " outside range 0 - " <<  sz-1 << endln;
    return VECTOR_NOT_VALID_ENTRY;
  }
#endif

  return theData[x];
}


// operator()(const ID &rows) const
//	Method to return a vector whose components are the components of the
//	current vector located in positions given by the ID rows.


Vector 
Vector::operator()(const ID &rows) const
{
  // create a new Vector to be returned
  Vector result(rows.Size());

  // check if obtained VEctor of correct size
  if (result.Size() != rows.Size()) {
    opserr << "Vector::()(ID) - new Vector could not be constructed\n";
    return result;
  }

  // copy the appropraite contents from current to result     
  int pos;
  for (int i=0; i<rows.Size(); i++) {
    pos = rows(i);
    if (pos <0 || pos >= sz) {
      opserr << "Vector::()(ID) - invalid location " << pos << " outside range [0, " << sz-1 << "]\n";
    } else
      result(i) = (*this)(pos);
  }
  return result;
}


// Vector &operator=(const Vector  &V):
//	the assignment operator, This is assigned to be a copy of V. if sizes
//	are not compatable this.theData [] is deleted. The data pointers will not
//	point to the same area in mem after the assignment.
//

Vector &
Vector::operator=(const Vector &V) 
{
  // first check we are not trying v = v
  if (this != &V) {

	  /*
#ifdef _G3DEBUG
    // check size compatability, if different warning
    if (sz != V.sz) 
      opserr << "Vector::operator=() - vectors of differing sizes\n");
#endif
	  */

      if (sz != V.sz)  {

#ifdef _G3DEBUG
	  opserr << "Vector::operator=() - vectors of differing sizes\n";
#endif

	  // Check that we are not deleting an empty Vector
	  if (this->theData != 0) delete [] this->theData;

	  this->sz = V.sz;
	  
	  // Check that we are not creating an empty Vector
	  theData = (sz != 0) ? new (nothrow) double[sz] : 0;
      }


      //	 copy the data
      for (int i=0; i<sz; i++)
	  theData[i] = V.theData[i];
  }

  return *this;
}




Vector &
Vector::operator=(const Tensor &V) 
{
  int rank = V.rank();
  if (rank != 2) {
    opserr << "Vector::operator=() - tensor must be of rank 2\n";
    return *this;
  }
  int dim = V.dim(1);
  if (dim != V.dim(2)) {
    opserr << "Vector::operator=() - tensor must have square dimensions\n";
    return *this;
  }
  
  if (dim != 2 || dim != 3 || dim != 1) {
    opserr << "Vector::operator=() - tensor must be of dimension 2 or 3\n";
    return *this;
  }      
  
  if (dim == 1) {
    if (sz != 1) {
      opserr << "Vector::operator=() - Vector size must be 1\n"; 
      return *this;
    }
    theData[0] = V.cval(1,1);
  } else if (dim == 2) {
    if (sz != 3) {
      opserr << "Vector::operator=() - Vector size must be 3\n"; 
      return *this;
    }
    theData[0] = V.cval(1,1);
    theData[1] = V.cval(2,2);
    theData[2] = V.cval(1,2);
  } else {
    if (sz != 6) {
      opserr << "Vector::operator=() - Vector size must be 6\n"; 
      return *this;
    }      
    theData[0] = V.cval(1,1);
    theData[1] = V.cval(2,2);
    theData[2] = V.cval(3,3);
    theData[3] = V.cval(1,2);
    theData[4] = V.cval(1,3);
    theData[5] = V.cval(2,3);
  }
  return *this;
}



// Vector &operator+=(double fact):
//	The += operator adds fact to each element of the vector, data[i] = data[i]+fact.

Vector &Vector::operator+=(double fact)
{
  if (fact != 0.0)
    for (int i=0; i<sz; i++)
      theData[i] += fact;
  return *this;    
}



// Vector &operator-=(double fact)
//	The -= operator subtracts fact from each element of the vector, data[i] = data[i]-fact.

Vector &Vector::operator-=(double fact)
{
  if (fact != 0.0)
    for (int i=0; i<sz; i++)
      theData[i] -= fact;
  return *this;    
}



// Vector &operator*=(double fact):
//	The -= operator subtracts fact from each element of the vector, data[i] = data[i]-fact.

Vector &Vector::operator*=(double fact)
{
  for (int i=0; i<sz; i++)
    theData[i] *= fact;
  return *this;
}



// Vector &operator/=(double fact):
//	The /= operator divides each element of the vector by fact, theData[i] = theData[i]/fact.
//	Program exits if divide-by-zero would occur with warning message.

Vector &Vector::operator/=(double fact)
{
  if (fact == 0.0) { // instead of divide-by-zero error set to VECTOR_VERY_LARGE_VALUE
    for (int i=0; i<sz; i++)
      theData[i] = VECTOR_VERY_LARGE_VALUE;
  } else {
    for (int i=0; i<sz; i++)
      theData[i] /= fact;
  }
  
  return *this;
}




// Vector operator+(double fact):
//	The + operator returns a Vector of the same size as current, whose components
//	are return(i) = theData[i]+fact;

Vector 
Vector::operator+(double fact) const
{
  Vector result(*this);
  if (result.Size() != sz) 
    opserr << "Vector::operator+(double) - ran out of memory for new Vector\n";

  result += fact;
  return result;
}



// Vector operator-(double fact):
//	The + operator returns a Vector of the same size as current, whose components
//	are return(i) = theData[i]-fact;

Vector 
Vector::operator-(double fact) const
{
    Vector result(*this);
    if (result.Size() != sz) 
      opserr << "Vector::operator-(double) - ran out of memory for new Vector\n";

    result -= fact;
    return result;
}



// Vector operator*(double fact):
//	The + operator returns a Vector of the same size as current, whose components
//	are return(i) = theData[i]*fact;

Vector 
Vector::operator*(double fact) const
{
    Vector result(*this);
    if (result.Size() != sz) 
      opserr << "Vector::operator*(double) - ran out of memory for new Vector\n";

    result *= fact;
    return result;
}


// Vector operator/(double fact):
//	The + operator returns a Vector of the same size as current, whose components
//	are return(i) = theData[i]/fact; Exits if divide-by-zero error.

Vector 
Vector::operator/(double fact) const
{
    if (fact == 0.0) 
      opserr << "Vector::operator/(double fact) - divide-by-zero error coming\n";

    Vector result(*this);
    if (result.Size() != sz) 
      opserr << "Vector::operator/(double) - ran out of memory for new Vector\n";

    result /= fact;
    return result;
}



// Vector &operator+=(const Vector &V):
//	The += operator adds V's data to data, data[i]+=V(i). A check to see if
//	vectors are of same size is performed if VECTOR_CHECK is defined.

Vector &
Vector::operator+=(const Vector &other)
{
#ifdef _G3DEBUG
  if (sz != other.sz) {
    opserr << "WARNING Vector::operator+=(Vector):Vectors not of same sizes: " << sz << " != " << other.sz << endln;
    return *this;
  }    
#endif

  for (int i=0; i<sz; i++)
    theData[i] += other.theData[i];
  return *this;	    
}



// Vector &operator-=(const Vector &V):
//	The -= operator subtracts V's data from  data, data[i]+=V(i). A check 
//   	to see if vectors are of same size is performed if VECTOR_CHECK is defined.

Vector &
Vector::operator-=(const Vector &other)
{
#ifdef _G3DEBUG
  if (sz != other.sz) {
    opserr << "WARNING Vector::operator+=(Vector):Vectors not of same sizes: " << sz << " != " << other.sz << endln;
    return *this;
  }
#endif
  
  for (int i=0; i<sz; i++)
    theData[i] -= other.theData[i];
  return *this;    
}



// Vector operator+(const Vector &V):
//	The + operator checks the two vectors are of the same size if VECTOR_CHECK is defined.
// 	Then returns a Vector whose components are the vector sum of current and V's data.

Vector 
Vector::operator+(const Vector &b) const
{
#ifdef _G3DEBUG
  if (sz != b.sz) {
    opserr << "WARNING Vector::operator+=(Vector):Vectors not of same sizes: " << sz << " != " << b.sz << endln;
    return *this;
  }
#endif

    Vector result(*this);

    // check new Vector of correct size
  if (result.Size() != sz) {
    opserr << "Vector::operator-(Vector): new Vector not of correct size \n";
    return result;
  }
  result += b;
  return result;
}


// Vector operator-(const Vector &V):
//	The - operator checks the two vectors are of the same size and then returns a Vector
//	whose components are the vector difference of current and V's data.

Vector 
Vector::operator-(const Vector &b) const
{
#ifdef _G3DEBUG
  if (sz != b.sz) {
    opserr << "WARNING Vector::operator+=(Vector):Vectors not of same sizes: " << sz << " != " << b.sz << endln;
    return *this;
  }
#endif

  Vector result(*this);

  // check new Vector of correct size
  if (result.Size() != sz) {
    opserr << "Vector::operator-(Vector): new Vector not of correct size \n";
    return result;
  }

  result -= b;
  return result;
}



// double operator^(const Vector &V) const;
//	Method to perform (Vector)transposed * vector.
double
Vector::operator^(const Vector &V) const
{
#ifdef _G3DEBUG
  if (sz != V.sz) {
    opserr << "WARNING Vector::operator+=(Vector):Vectors not of same sizes: " << sz << " != " << V.sz << endln;
    return 0.0;
  }
#endif

  double result = 0.0;
  double *dataThis = theData;
  double *dataV = V.theData;
    for (int i=0; i<sz; i++)
      result += *dataThis++ * *dataV++;

    return result;
}


// Vector operator/(const Matrix &M) const;    
//	Method to return inv(M)*this

Vector
Vector::operator/(const Matrix &M) const
{
  Vector res(M.noRows());
    
  if (M.noRows() != M.noCols()) { // if not square do least squares solution
    Matrix A(M^M);
    A.Solve(*this, res);    
  }
  else {
    M.Solve(*this, res);
  }
  return res;
}
    
	
// Vector operator==(const Vector &V):
//	The == operator checks the two vectors are of the same size if VECTOR_CHECK is defined.
// 	Then returns 1 if all the components of the two vectors are equal and 0 otherwise.

int 
Vector::operator==(const Vector &V) const
{
  if (sz != V.sz)
    return 0;

  double *dataThis = theData;
  double *dataV = V.theData;

  for (int i=0; i<sz; i++)
    if (*dataThis++ != *dataV++)
      return 0;

  return 1;
}

int 
Vector::operator==(double value) const
{
  double *dataThis = theData;

  for (int i=0; i<sz; i++)
    if (*dataThis++ != value)
      return 0;

  return 1;
}


// Vector operator!=(const Vector &V):
//	The != operator checks the two vectors are of the same size if VECTOR_CHECK is defined.
// 	Then returns 1 if any of the components of the two vectors are unequal and 0 otherwise.

int
Vector::operator!=(const Vector &V) const
{
  if (sz != V.sz)
    return 1;

  double *dataThis = theData;
  double *dataV = V.theData;

  for (int i=0; i<sz; i++)
    if (*dataThis++ != *dataV++)
      return 1;

  return 0;
}

int 
Vector::operator!=(double value) const
{
  double *dataThis = theData;

  for (int i=0; i<sz; i++)
    if (*dataThis++ != value)
      return 1;

  return 0;
}


// friend OPS_Stream &operator<<(OPS_Stream &s, const Vector &V)
//	A function is defined to allow user to print the vectors using OPS_Streams.

OPS_Stream &operator<<(OPS_Stream &s, const Vector &V)
{
  /*
  for (int i=0; i<V.Size(); i++) 
      s << V(i) << " ";

  return s << endln;
  */
  return s.write(V.theData, V.sz);

}

// friend istream &operator>>(istream &s, Vector &V)
//	A function is defined to allow user to input the data into a Vector which has already
//	been constructed with data, i.e. Vector(int) or Vector(const Vector &) constructors.

/*
istream &operator>>(istream &s, Vector &V)
{
  for (int i=0; i<V.Size(); i++) 
    s >> V(i);
  
    return s;
}
*/


Vector operator*(double a, const Vector &V)
{
  return V * a;
}


int
Vector::Assemble(const Vector &V, int init_pos, double fact) 
{
  int res = 0;
  int cur_pos   = init_pos;  
  int final_pos = init_pos + V.sz - 1;
  
  if ((init_pos >= 0) && (final_pos < sz))
  {
     for (int j=0; j<V.sz; j++) 
        (*this)(cur_pos++) += V(j)*fact;
  }
  else 
  {
     opserr << "WARNING: Vector::Assemble(const Vector &V, int init_pos, double fact): ";
     opserr << "position outside bounds \n";
     res = -1;
  }

  return res;
}




int
Vector::Extract(const Vector &V, int init_pos, double fact) 
{
  int res = 0;
  int cur_pos   = init_pos;  
  int final_pos = init_pos + sz - 1;
  
  if ((init_pos >= 0) && (final_pos < V.sz))
  {
     for (int j=0; j<sz; j++) 
        (*this)(j) = V(cur_pos++)*fact;
  }
  else 
  {
     opserr << "WARNING: Vector::Assemble(const Vector &V, int init_pos, double fact): ";
     opserr << "position outside bounds \n";
     res = -1;
  }

  return res;
}


Matrix Vector::operator%(const Vector &V) const
{
  // if sizes are compatable add
#ifdef _G3DEBUG
  if (sz != V.sz) {
    // else sizes are incompatable, do nothing but warning
    opserr << "WARNING Vector::tensor multiplication operator % - incompatable Vector sizes\n";
    return -1;
  }
#endif

  //we want result=a*b'

  Matrix result(sz, sz);
  for (int i = 0; i<sz; i++)
    for (int j = 0; j<sz; j++)
      result(i, j) = theData[i] * V.theData[j];

  return result;

}

//added by neallee@tju.edu.cn
void
Vector::computePrincipalValues(Vector &answer)
{
  //
  // This function computes sorted principal values of a strain vector.
  // Engineering notation is used.
  //

  if (sz == 1) {
    answer.resize(1);
    answer.Zero();
    answer.setData(theData, 1);
  }
  else if (sz == 3) {
    answer.resize(2);
    answer.Zero();
    // 2D problem
    double ast, dst, D;

    ast = theData[0] + theData[1];
    dst = theData[0] - theData[1];
    D = sqrt(dst * dst + theData[2] * theData[2]);
    answer(1) = 0.5 * (ast + D);
    answer(2) = 0.5 * (ast - D);
  }
  else {
    answer.resize(3);
    answer.Zero();
    // 3D problem
    double I1 = 0.0, I2 = 0.0, I3 = 0.0, s1, s2, s3;
    Matrix s;
    int nonzeroFlag = 0;

    /*this->convertToFullForm(s);
    for (int i = 1; i <= size; i++) {
    if (fabs(s.at(i)) > 1.e-20) {
    nonzeroFlag = 1;
    }
    }*/
    if (sz == 6) {
      I1 = theData[0] + theData[1] + theData[2];
      I2 = theData[0] * theData[1] + theData[1] * theData[2] + theData[2] * theData[0] -
        0.25 * (theData[3] * theData[3] + theData[4] * theData[4] + theData[5] * theData[5]);
      I3 = theData[0] * theData[1] * theData[2] +
        0.25 * (theData[3] * theData[4] * theData[5] - theData[0] * theData[3] * theData[3] -
        theData[1] * theData[4] * theData[4] - theData[2] * theData[5] * theData[5]);

      /*
      * Call cubic3r to ensure, that all three real eigenvalues will be found, because we have symmetric tensor.
      * This allows to overcome various rounding errors when solving general cubic equation.
      */
      int i;
      cubic3r((double)-1., I1, -I2, I3, &s1, &s2, &s3, &i);

      if (i > 0) {
        answer(0) = s1;
      }

      if (i > 1) {
        answer(1) = s2;
      }

      if (i > 2) {
        answer(2) = s3;
      }

      // sort results
      for (int i = 1; i < 3; i++) {
        for (int j = 1; j < 3; j++) {
          if (answer(j) > answer(j-1)) {
            std::swap(answer(j), answer(j-1));
          }
        }
      }
    }
    else {
      opserr << "Vector::computePrincipalValues(), 3D mat principle value error!!" << endln;
    }
  }
}

#define CUBIC_ZERO 1.0e-100
#define sgn(x) ((x) < 0 ? -1.0 : 1.0)
void
Vector::cubic3r(double a, double b, double c, double d, double *r1, double *r2, double *r3, int *num)
//
//   solves cubic equation for real roots
//
//   input:
//     a,b,c,d - coefficients of equation in form:
//     ax^3 + bx^2 + cx + d = 0
//
//   output:
//     r1,r2,r3 - roots (only first num roots is valid)
//     num      - number of roots resolved
//
{
  double aa, r, q, p, D, phi;
  double help;

  if (fabs(a) < CUBIC_ZERO) {
    if (fabs(b) < CUBIC_ZERO) {
      if (fabs(c) < CUBIC_ZERO) {
        *num = 0;
        return;
      }
      else {
        *r1 = -d / c;
        *num = 1;
        return;
      }
    }
    else { //fabs(b) < CUBIC_ZERO
      if ((D = c * c - 4.0 * b * d) < 0.0) {
        *num = 0;
        return;
      }
      else {
        //*r1 = (-c + sqrt(D)) / 2.0 / b;
        //*r2 = (-c - sqrt(D)) / 2.0 / b;
        if (fabs(c) < CUBIC_ZERO) {
          help = -d / b;
          if (help >= 0.) {
            *r1 = sqrt(help);
            *r2 = -sqrt(help);
            *num = 2;
            return;
          }
          else {
            *num = 0;
            return;
          }
        }
        else {
          help = -(c + sgn(c) * sqrt(D)) / 2.0;
          *r1 = help / b;
          *r2 = d / help;
          *num = 2;
          return;
        }
      }
    }
  }
  else {
    aa = a;
    a = b / aa;
    b = c / aa;
    c = d / aa;

    q = (a * a - 3.0 * b) / 9.0;
    r = (2.0 * a * a * a - 9.0 * a * b + 27.0 * c) / 54.0;

    D = q * q * q - r * r;
    //  if (D > 0.) {
    // three real roots
    help = r / sqrt(q * q * q);
    if (fabs(help) > 1.0) {
      help = sgn(help);                // prevent rounding errors
    }

    phi = acos(help);
    p = sqrt(q);

    *r1 = -2.0 *p *cos(phi / 3.0) - a / 3.0;
    *r2 = -2.0 *p *cos((phi + 2.0 *PI) / 3.0) - a / 3.0;
    *r3 = -2.0 *p *cos((phi - 2.0 *PI) / 3.0) - a / 3.0;
    *num = 3;
    /*  } else {
    *
    * help = fabs(r) + sqrt(-D);
    * A = -sgn(r)*pow(D, 1./3.);
    * if (fabs(A) > CUBIC_ZERO)
    * B = q/A;
    * else
    * B = 0.0;
    *
    *****r1 = (A+B) - a/3.0;
    *****num = 1;
    * }
    */
    return;
  }
}

void
Vector::computePrincipalValDir(Vector &answer, Matrix &dir) const
{
  //
  // This function cumputes Principal values & directions corresponding to receiver.
  //
  // Return Values:
  //
  // matrix dir -> principal directions of strains or stresses
  // array sp -> principal strains or stresses
  //

  Matrix ss;
  Vector sp;
  int nval;
  int nonzeroFlag = 0;
  int size = sz;


  if (sz == 1) { // _1dMat
    answer = *this;
    dir.resize(1, 1);
    dir(0, 0) = 1.0;
    return;
  }
  else if (sz == 3) { // _PlaneStress
    // 2D problem
    ss.resize(2, 2);
    answer.resize(2);

    for (int i = 0; i < size; i++) {
      if (fabs(theData[i]) > 1.e-20) {
        nonzeroFlag = 1;
      }
    }

    if (nonzeroFlag == 0) {
      answer.Zero();
      ss.Zero();
      return;
    }

    ss(0, 0) = theData[0];
    ss(1, 1) = theData[1];
    ss(0, 1) = ss(1, 0) = theData[2];
  }
  else {
    // 3D problem
    ss.resize(3, 3);
    //Vector s;

    answer.resize(3);

    if (sz == 6) {
      //this->convertToFullForm(s);
      for (int i = 0; i < size; i++) {
        if (fabs(theData[i]) > 1.e-20) {
          nonzeroFlag = 1;
        }
      }

      if (nonzeroFlag == 0) {
        answer.Zero();
        ss.Zero();
        return;
      }

      ss(0, 0) = theData[0];
      ss(1, 1) = theData[1];
      ss(2, 2) = theData[2];
      ss(0, 1) = ss(1, 0) = theData[5];
      ss(0, 2) = ss(2, 0) = theData[4];
      ss(1, 2) = ss(2, 1) = theData[3];
    }
    else {
      opserr << "Vector::computePrincipalValDir(), 3D mat principle value error!!" << endln;
    }
  }

  ss.jaco_(answer, dir, 10);

  // sort results
  nval = 3;
  if (sz == 3) {
    nval = 2;
  }

  for (int ii = 0; ii < nval-1; ii++) {
    for (int jj = 0; jj < nval-1; jj++) {
      if (answer(jj + 1) > answer(jj)) {
        // swap eigen values and eigen vectors
        std::swap(answer(jj + 1), answer(jj));
        for (int kk = 0; kk < nval; kk++) {
          std::swap(dir(kk, jj + 1), dir(kk, jj));
        }
      }
    }
  }
}

void
Vector::computeDeviatoricVolumetricSplit(Vector &dev, double &vol) const
{

  if (sz == 1) { //_1dMat
    // 1D model
    opserr << "Vector::computeDeviatoricVolumetricSplit: No Split for 1D!" << endln;
    //  dev.resize(1); dev.at(1) = 0.0;
    // vol = this->at (1);
  }
  else if (sz == 3) {
    // plane stress problem
    opserr << "Vector::computeDeviatoricVolumetricSplit: No Split for plane stress!" << endln;
    //  dev = *this;
    // vol = (this->at(1)+this->at(2))/2.0;
    //    dev.at(1) -= vol;
    // dev.at(2) -= vol;
  }
  else {
    // 3d, plane strain or axisymmetric problem
    dev = *this;
    vol = (theData[0] + theData[1] + theData[2]) / 3.0;
    dev(1) -= vol;
    dev(2) -= vol;
    dev(3) -= vol;
  }
}

double
Vector::computeSecondCoordinate() const
{
  //
  // This function computes the second Haigh-Westergaard coordinate
  // from the deviatoric stress state
  //
  ;
  return sqrt(2. * computeSecondInvariant());
}

double
Vector::computeSecondInvariant() const
{
  //
  // This function computes the second invariant of
  // of the deviatoric stress vector
  //
  if (sz == 1) { //_1dMat
    // 1d problem
    return .5 * theData[0] * theData[0];
  }
  else if (sz == 3) { //_PlaneStress
    // 2d problem: plane stress
    return .5 * (theData[0] * theData[0] + theData[1] * theData[1]) + theData[2] * theData[2];
  }
  else if (sz == 4) { //_PlaneStrain
    //  plane strain or axisymmetry
    return .5 * (theData[0] * theData[0] + theData[1] * theData[1] + theData[2] * theData[2]) +
      theData[3] * theData[3];
  }
  else {
    // 3d problem
    return .5 * (theData[0] * theData[0] + theData[1] * theData[1] + theData[2] * theData[2]) +
      theData[3] * theData[3] + theData[4] * theData[4] + theData[5] * theData[5];
  }
}

double
Vector::computeStrainNorm() const
{
  //
  // This function computes the tensorial Norm of the strain in engineering notation
  //
  if (sz == 1) { //_1dMat
    // 1d problem
    return sqrt(theData[0] * theData[0]);
  }
  else if (sz == 3) { //_PlaneStress
    // 2d problem: plane stress
    return sqrt(.5 * (2. * theData[0] * theData[0] + 2 * theData[1] * theData[1] + theData[2] * theData[2]));
  }
  else if (sz == 4) { //_PlaneStrain
    //  plane strain or axisymmetry
    return sqrt(.5 * (2. * theData[0] * theData[0] + 2. * theData[1] * theData[1] + 2. * theData[2] * theData[2] +
      theData[3] * theData[3]));
  }
  else {
    // 3d problem
    return sqrt(.5 * (2. * theData[0] * theData[0] + 2. * theData[1] * theData[1] + 2. * theData[2] * theData[2] +
      theData[3] * theData[3] + theData[4] * theData[4] + theData[5] * theData[5]));
  }
}


double
Vector::computeThirdCoordinate() const
{
  //
  // This function computes the third Haigh-Westergaard coordinate
  // from the deviatoric stress state
  //
  double c1 = 0.0;
  if (computeSecondInvariant() == 0.) {
    c1 = 0.0;
  }
  else {
    c1 = (3. * sqrt(3.) / 2.) * computeThirdInvariant() / (pow(computeSecondInvariant(), (3. / 2.)));
  }

  if (c1 > 1.0) {
    c1 = 1.0;
  }

  if (c1 < -1.0) {
    c1 = -1.0;
  }

  return 1. / 3. * acos(c1);
}

double
Vector::computeThirdInvariant() const
{
  //
  // This function computes the third invariant
  // of the deviatoric stress vector.
  //
  if (sz == 1) { //_1dMat
    // 1d problem
    return (1. / 3.) * theData[0] * theData[0] * theData[0];
  }
  else if (sz == 3) { //_PlaneStress
    // 2d problem: plane stress
    return (1. / 3.) * (theData[0] * theData[0] * theData[0] + 3. * theData[1] * theData[2] * theData[2]
      + 3. * theData[0] * theData[2] * theData[2] + theData[1] * theData[1] * theData[1]);
  }
  else if (sz == 4) { //_PlaneStrain
    // plane strain or axisymmetry
    return (1. / 3.) * (theData[0] * theData[0] * theData[0] + 3. * theData[0] * theData[3] * theData[3] +
      3. * theData[1] * theData[3] * theData[3] + theData[1] * theData[1] * theData[1] +
      theData[2] * theData[2] * theData[2]);
  }
  else {
    // 3d problem
    return (1. / 3.) * (theData[0] * theData[0] * theData[0] + 3. * theData[0] * theData[5] * theData[5] +
      3. * theData[0] * theData[4] * theData[4] + 6. * theData[3] * theData[5] * theData[4] +
      3. * theData[1] * theData[5] * theData[5] + 3 * theData[2] * theData[4] * theData[4] +
      theData[1] * theData[1] * theData[1] + 3. * theData[1] * theData[3] * theData[3] +
      3. * theData[2] * theData[3] * theData[3] + theData[2] * theData[2] * theData[2]);
  }
}
