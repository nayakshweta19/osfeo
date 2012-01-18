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

// $Revision: 1.4 $
// $Date: 2006/02/28 19:19:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/ritzSOE/FullGenRitzSOE.cpp,v $


// File: ~/system_of_eqn/linearSOE/fullGEN/FullGenRitzSOE.C
//
// Written: fmk 
// Created: Febuary 1997
// Revision: A
//
// Description: This file contains the implementation for FullGenRitzSOE


#include <FullGenRitzSOE.h>
#include <FullGenRitzSOESolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <f2c.h>
#include <math.h>


#include <Channel.h>
#include <FEM_ObjectBroker.h>

#define RitzSOE_TAGS_FullGenRitzSOE    25

FullGenRitzSOE::FullGenRitzSOE(FullGenRitzSOESolver &theSolvr)
:RitzSOE(theSolvr, RitzSOE_TAGS_FullGenRitzSOE),
size(0), A(0), B(0), M(0), X(0), vectX(0), vectB(0), 
Asize(0), Bsize(0), Msize(0),
factored(false)
{
	theSolvr.setRitzSOE(*this);
}


FullGenRitzSOE::FullGenRitzSOE(int N, FullGenRitzSOESolver &theSolvr)
:RitzSOE(theSolvr, RitzSOE_TAGS_FullGenRitzSOE),
size(0), A(0), B(0), M(0), X(0), vectX(0), vectB(0),
Asize(0), Bsize(0), Msize(0),
factored(false)
{
	size = N;

	A = new double[size*size];

	if (A == 0) {
		opserr << "WARNING :FullGenRitzSOE::FullGenRitzSOE :";
		opserr << " ran out of memory for A (size,size) (";
		opserr << size <<", " << size << ") \n";
		size = 0; 
	} else {
		// zero the matrix
		Asize = size*size;
		for (int i=0; i<Asize; i++)
			A[i] = 0;

		B = new double[size];
		X = new double[size];

		if (B == 0 || X == 0) {
			opserr << "WARNING :FullGenRitzSOE::FullGenRitzSOE :";
			opserr << " ran out of memory for vectors (size) (";
			opserr << size << ") \n";
			size = 0; Bsize = 0;
		} else {
			Bsize = size;
			// zero the vectors
			for (int j=0; j<size; j++) {
				B[j] = 0;
				X[j] = 0;
			}
		}
	}

	vectX = new Vector(X,size);
	vectB = new Vector(B,size);    

	theSolvr.setRitzSOE(*this);

	// invoke setSize() on the Solver        
	if (theSolvr.setSize() < 0) {
		opserr << "WARNING :FullGenRitzSOE::FullGenRitzSOE :";
		opserr << " solver failed setSize() in constructor\n";
	}    

}


FullGenRitzSOE::~FullGenRitzSOE()
{
	if (A != 0) delete [] A;
	if (B != 0) delete [] B;
	if (M != 0) delete [] M;
	if (X != 0) delete [] X;
	if (vectX != 0) delete vectX;    
	if (vectB != 0) delete vectB;        
}


int
FullGenRitzSOE::getNumEqn(void) const
{
	return size;
}

int 
FullGenRitzSOE::setSize(Graph &theGraph)
{
	int result = 0;
	int oldSize = size;
	size = theGraph.getNumVertex();

	if (size*size > Asize) { // we have to get another space for A

		if (A != 0) 
			delete [] A;

		A = new double[size*size];

		if (A == 0) {
			opserr << "WARNING FullGenRitzSOE::FullGenRitzSOE :";
			opserr << " ran out of memory for A (size,size) (";
			opserr << size <<", " << size << ") \n";
			size = 0; Asize = 0;
			result =  -1;
		} else
			Asize = size*size;
	}

	// zero the matrix
	for (int i=0; i<Asize; i++)
		A[i] = 0;

	factored = false;

	// we have to get another space for M
	if (size*size > Msize) {
		if (M != 0) 
			delete [] M;

		M = new double[size*size];
		if (M == 0) {
			opserr << "WARNING FullGenEigenSOE::setSize() - "
				<< "ran out of memory for M (size,size) ("
				<< size << ", " << size << ")\n";
			Msize = 0; size = 0;
			result= -1;
		}
		else
			Msize = size*size;
	}

	// zero the matrix
	for (int i=0; i<Msize; i++)
		M[i] = 0;

	if (size > Bsize) { // we have to get space for the vectors

		// delete the old	
		if (B != 0) delete [] B;
		if (X != 0) delete [] X;

		// create the new
		B = new double[size];
		X = new double[size];

		if (B == 0 || X == 0) {
			opserr << "WARNING FullGenRitzSOE::FullGenRitzSOE :";
			opserr << " ran out of memory for vectors (size) (";
			opserr << size << ") \n";
			size = 0; Bsize = 0;
			result =  -1;
		}
		else
			Bsize = size;
	}

	// zero the vectors
	for (int j=0; j<Bsize; j++) {
		B[j] = 0;
		X[j] = 0;
	}

	// create new Vectors
	if (size != oldSize) {
		if (vectX != 0)
			delete vectX;

		if (vectB != 0)
			delete vectB;

		vectX = new Vector(X,Bsize);
		vectB = new Vector(B,Bsize);	

	}

	// invoke setSize() on the Solver    
	RitzSOESolver *theSolvr = this->getSolver();
	int solverOK = theSolvr->setSize();
	if (solverOK < 0) {
		opserr << "WARNING:FullGenRitzSOE::setSize :";
		opserr << " solver failed setSize()\n";
		return solverOK;
	}    

	return result;
}

int 
FullGenRitzSOE::addA(const Matrix &m, const ID &id, double fact)
{
	// check for a quick return 
	if (fact == 0.0)  return 0;

	int idSize = id.Size();

	// check that m and id are of similar size
	if (idSize != m.noRows() && idSize != m.noCols()) {
		opserr << "FullGenRitzSOE::addA()	- Matrix and ID not of similar sizes\n";
		return -1;
	}

	if (fact == 1.0) { // do not need to multiply 
		for (int i=0; i<idSize; i++) {
			int col = id(i);
			if (col < size && col >= 0) {
				double *startColiPtr = A + col*size;
				for (int j=0; j<idSize; j++) {
					int row = id(j);
					if (row <size && row >= 0) {
						double *APtr = startColiPtr + row;
						*APtr += m(j,i);
					}
				}  // for j
			} 
		}  // for i
	} else {
		for (int i=0; i<idSize; i++) {
			int col = id(i);
			if (col < size && col >= 0) {
				double *startColiPtr = A + col*size;
				for (int j=0; j<idSize; j++) {
					int row = id(j);
					if (row <size && row >= 0) {
						double *APtr = startColiPtr + row;
						*APtr += m(j,i) * fact;
					}
				}  // for j
			} 
		}  // for i
	}    
	return 0;
}




int 
FullGenRitzSOE::addB(const Vector &v, const ID &id, double fact)
{
	// check for a quick return 
	if (fact == 0.0)  return 0;

	int idSize = id.Size();    
	// check that m and id are of similar size
	if (idSize != v.Size() ) {
		opserr << "FullGenRitzSOE::addB()	- Vector and ID not of similar sizes\n";
		return -1;
	}    

	if (fact == 1.0) { // do not need to multiply if fact == 1.0
		for (int i=0; i<idSize; i++) {
			int pos = id(i);
			if (pos <size && pos >= 0)
				B[pos] += v(i);
		}
	} else if (fact == -1.0) { // do not need to multiply if fact == -1.0
		for (int i=0; i<idSize; i++) {
			int pos = id(i);
			if (pos <size && pos >= 0)
				B[pos] -= v(i);
		}
	} else {
		for (int i=0; i<idSize; i++) {
			int pos = id(i);
			if (pos <size && pos >= 0)
				B[pos] += v(i) * fact;
		}
	}	
	return 0;
}



int
FullGenRitzSOE::setB(const Vector &v, double fact)
{
	// check for a quick return 
	if (fact == 0.0)  return 0;


	if (v.Size() != size) {
		opserr << "WARNING BandGenLinSOE::setB() -";
		opserr << " incomptable sizes " << size << " and " << v.Size() << endln;
		return -1;
	}

	if (fact == 1.0) { // do not need to multiply if fact == 1.0
		for (int i=0; i<size; i++) {
			B[i] = v(i);
		}
	} else if (fact == -1.0) {
		for (int i=0; i<size; i++) {
			B[i] = -v(i);
		}
	} else {
		for (int i=0; i<size; i++) {
			B[i] = v(i) * fact;
		}
	}	
	return 0;
}

void 
FullGenRitzSOE::zeroA(void)
{
	double *Aptr = A;
	int theSize = size*size;
	for (int i=0; i<theSize; i++)
		*Aptr++ = 0;

	factored = false;
}

void 
FullGenRitzSOE::zeroB(void)
{
	double *Bptr = B;
	for (int i=0; i<size; i++)
		*Bptr++ = 0;
}

void 
FullGenRitzSOE::setX(int loc, double value)
{
	if (loc < size && loc >=0)
		X[loc] = value;
}

void 
FullGenRitzSOE::setX(const Vector &x)
{
	if (x.Size() == size && vectX != 0)
		*vectX = x;
}

const Vector &
FullGenRitzSOE::getX(void)
{
	if (vectX == 0) {
		opserr << "FATAL FullGenRitzSOE::getX - vectX == 0";
		exit(-1);
	}
	return *vectX;
}

const Vector &
FullGenRitzSOE::getB(void)
{
	if (vectB == 0) {
		opserr << "FATAL FullGenRitzSOE::getB - vectB == 0";
		exit(-1);
	}        
	return *vectB;
}

const Matrix&
FullGenRitzSOE::getK(void)
{
	matA = new Matrix(A, Bsize, Bsize);	
	if (matA == 0) {
		opserr << "FATAL FullGenRitzSOE::getK - matA == 0";
		exit(-1);
	}        
	return *matA;
}

const Matrix&
FullGenRitzSOE::getM(void)
{
	matA = new Matrix(M, Bsize, Bsize);	
	if (matA == 0) {
		opserr << "FATAL FullGenRitzSOE::getM - matA == 0";
		exit(-1);
		}        
	return *matA;
}

double 
FullGenRitzSOE::normRHS(void)
{
	double norm =0.0;
	for (int i=0; i<size; i++) {
		double Yi = B[i];
		norm += Yi*Yi;
	}
	return sqrt(norm);

}    


int
FullGenRitzSOE::setFullGenSolver(FullGenRitzSOESolver &newSolver)
{
	newSolver.setRitzSOE(*this);

	if (size != 0) {
		int solverOK = newSolver.setSize();
		if (solverOK < 0) {
			opserr << "WARNING:FullGenRitzSOE::setSolver :";
			opserr << "the new solver could not setSeize() - staying with old\n";
			return -1;
		}
	}

	return this->RitzSOE::setSolver(newSolver);
}



int FullGenRitzSOE::addM( const Matrix &m, const ID &id, double fact )
{
	    // check for quick return 
    if (fact == 0.0)
        return 0;

    // check that m and id are of similar size
    int idSize = id.Size();    
    if (idSize != m.noRows() && idSize != m.noCols()) {
        opserr << "FullGenEigenSOE::addM() - Matrix and ID not of similar sizes\n";
        return -1;
    }

    if (fact == 1.0) { // do not need to multiply 
        for (int i=0; i<idSize; i++) {
            int col = id(i);
            if (col < size && col >= 0) {
                double *startColiPtr = M + col*size;
                for (int j=0; j<idSize; j++) {
                    int row = id(j);
                    if (row <size && row >= 0) {
                        double *MPtr = startColiPtr + row;
                        *MPtr += m(j,i);
                    }
                }  // for j
            } 
        }  // for i
    } else {
        for (int i=0; i<idSize; i++) {
            int col = id(i);
            if (col < size && col >= 0) {
                double *startColiPtr = M + col*size;
                for (int j=0; j<idSize; j++) {
                    int row = id(j);
                    if (row <size && row >= 0) {
                        double *MPtr = startColiPtr + row;
                        *MPtr += m(j,i)*fact;
                    }
                }  // for j
            } 
        }  // for i
    }

    return 0;
}

void FullGenRitzSOE::zeroM( void )
{
	double *Mptr = M;
	for (int i=0; i<size; i++)
		*Mptr++ = 0;

	factored = false;
}

int
FullGenRitzSOE::sendSelf(int commitTag, Channel &theChannel)

{
	// nothing to do
	return 0;
}

int
FullGenRitzSOE::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	// nothing to do
	return 0;
}



