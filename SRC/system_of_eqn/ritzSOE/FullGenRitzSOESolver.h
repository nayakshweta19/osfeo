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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/ritzSOE/FullGenRitzSOESolver.h,v $


#ifndef FullGenRitzSOESolver_h
#define FullGenRitzSOESolver_h

// File: ~/system_of_eqn/ritzSOE/FullGenRitzSOESolver.h
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for FullGenRitzSOESolver.
// FullGenRitzSOESolver is a concrete subclass of LinearSOE. It stores full
// unsymmetric linear system of equations using 1d arrays in Fortan style
//
// What: "@(#) FullGenRitzSOESolver.h, revA"

#include <RitzSOESolver.h>
#include <FullGenRitzSOE.h>

class FullGenRitzSOESolver : public RitzSOESolver
{
public:
	FullGenRitzSOESolver(int classTag);    
	virtual ~FullGenRitzSOESolver();
   
	int setRitzSOE(FullGenRitzSOE &theSOE);

	FullGenRitzSOE *theSOE;

private:

};

#endif

