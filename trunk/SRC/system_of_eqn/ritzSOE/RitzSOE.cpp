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

// $Revision: 1.3 $
// $Date: 2001/07/20 22:36:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/ritzSOE/RitzSOE.cpp,v $


// File: ~/system_of_eqn/RitzSOE.C
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the implementation of RitzSOE.
//
// What: "@(#) RitzSOE.C, revA"

#include<RitzSOE.h>
#include<RitzSOESolver.h>

RitzSOE::RitzSOE(RitzSOESolver &theRitzSOESolver, int classtag)
:MovableObject(classtag), theSolver(&theRitzSOESolver)
{

}

RitzSOE::~RitzSOE()
{
	delete theSolver;
}

int 
RitzSOE::solve(void)
{
	opserr << "ERROR RitzSOE::solve(void) - need to specify numRitzs\n";
	return -1;
}

int RitzSOE::solve( int numRitz )
{
	return theSolver->solve(numRitz);
}


double
RitzSOE::getDeterminant(void)
{
	return theSolver->getDeterminant();
}



int 
RitzSOE::setSolver(RitzSOESolver &newSolver)
{
	theSolver = &newSolver;
	return 0;
}

RitzSOESolver *
RitzSOE::getSolver(void)
{
	return theSolver;
}

const Vector& RitzSOE::getRitzvector(int numRitz)
{
	return theSolver->getRitzvector(numRitz);
}

double RitzSOE::getRitzvalue(int numRitz)
{
	return theSolver->getRitzvalue(numRitz);
}