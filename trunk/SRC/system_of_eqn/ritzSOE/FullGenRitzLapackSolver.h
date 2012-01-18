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

// $Revision: 1.1.1.1 $
// $Date: 2000/09/15 08:23:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/ritzSOE/FullGenRitzLapackSolver.h,v $


// File: ~/system_of_eqn/ritzSOE/FullGenRitzLapackSolver.h
//
// Written: fmk 
// Created: Tue Sep 26 16:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// FullGenRitzLapackSolver. It solves the FullGenLinSOE object by calling
// Lapack routines.
//
// What: "@(#) FullGenRitzLapackSolver.h, revA"

#ifndef FullGenRitzLapackSolver_h
#define FullGenRitzLapackSolver_h

#include <FullGenRitzSOESolver.h>

class FullGenRitzLapackSolver : public FullGenRitzSOESolver
{
public:
	FullGenRitzLapackSolver();    
	~FullGenRitzLapackSolver();

	int solve(void) {return this->solve(theSOE->size);};
	int solve(int numRitz);
	int setSize(void);
	const Vector& getRitzvector(int mode);
	double getRitzvalue(int mode);

	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

protected:

private:
	void sort(int length, double *x, int *id);

	int *iPiv;
	int sizeIpiv;
	int numRitz;

	double *eigenvalue;
	double *eigenvector;
	int *sortingID;
	Vector *eigenV;
};

#endif

