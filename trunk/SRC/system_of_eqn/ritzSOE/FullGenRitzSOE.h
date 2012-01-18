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

// $Revision: 1.2 $
// $Date: 2001/12/07 00:17:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/ritzSOE/FullGenRitzSOE.h,v $


#ifndef FullGenRitzSOE_h
#define FullGenRitzSOE_h

// File: ~/system_of_eqn/ritzSOE/FullGenRitzSOE.h
//
// Written: fmk 
// Created: 02/97
// Revision: A
//
// Description: This file contains the class definition for FullGenRitzSOE
// FullGenRitzSOE is a subclass of LinearSOE. It stores all the components
// of the linear system of equation in 1d arrays.
//
// What: "@(#) FullGenRitzSOE.h, revA"

#include <RitzSOE.h>
#include <Vector.h>

class FullGenRitzSOESolver;

class FullGenRitzSOE : public RitzSOE
{
public:
	FullGenRitzSOE(FullGenRitzSOESolver &theSolver);        
	FullGenRitzSOE(int N, FullGenRitzSOESolver &theSolver);        

	~FullGenRitzSOE();

	int getNumEqn(void) const;
	int setSize(Graph &theGraph);
	int addA(const Matrix &, const ID &, double fact = 1.0);
	int addB(const Vector &, const ID &, double fact = 1.0);    
	int setB(const Vector &, double fact = 1.0);  
	int addM(const Matrix &, const ID &, double fact = 1.0);

	void zeroA(void);
	void zeroB(void);
	void zeroM(void);

	const Vector &getX(void);
	const Vector &getB(void);
	const Matrix &getK(void);
	const Matrix &getM(void);
	double normRHS(void);

	void setX(int loc, double value);        
	void setX(const Vector &x);        

	int setFullGenSolver(FullGenRitzSOESolver &newSolver);    

	friend class FullGenRitzLapackSolver;    

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

protected:

private:
	int size;    
	double *A, *B, *X, *M;
	Vector *vectX;
	Vector *vectB;
	Matrix *matA;
	int Asize, Bsize, Msize;
	bool factored;
};


#endif

