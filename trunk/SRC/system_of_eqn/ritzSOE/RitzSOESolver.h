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
// $Date: 2001/02/17 06:32:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/RitzSOESolver.h,v $


// File: ~/system_of_eqn/linearSOE/RitzSOESolver.h
//
// Written: fmk 
// Created: Tue Sep 26 16:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for RitzSOESolver.
// RitzSOESolver is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes.  Instances of RitzSOESolver are used
// to solve a system of equations.
//
// What: "@(#) RitzSOESolver.h, revA"

#ifndef RitzSOESolver_h
#define RitzSOESolver_h

#include <Solver.h>
#include <MovableObject.h>
#include <Vector.h>

class RitzSOE;

class RitzSOESolver : public MovableObject
{
public:
	RitzSOESolver(int classTag);    
	virtual ~RitzSOESolver();

	virtual int solve(void) = 0;
	virtual int solve(int numRitz) = 0;
	//    virtual int setLinearSOE(LinearSOE &theSOE) =0;
	virtual const Vector &getRitzvector(int numRitz) = 0;
	virtual double getRitzvalue(int numRitz) = 0; 

	virtual int setSize(void) = 0;

	virtual double getDeterminant(void) {return 1.0;};

protected:

private:

};

#endif

