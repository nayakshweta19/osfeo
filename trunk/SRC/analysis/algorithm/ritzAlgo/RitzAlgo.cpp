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
// $Date: 2003/02/14 23:00:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/eigenAlgo/RitzAlgo.cpp,v $


// File: ~/analysis/algorithm/eigenAlgo/RitzAlgo.C
//
// Written: Jun Peng
// Created: Mon Feb. 8, 1999
// Revision: A
//
// Description: This file contains the class definition of RitzAlgo.
// RitzAlgo is a class which performs a eigen solution algorithm
// to solve the Generalized eigen equations. It is not expected that 
// this class will have subclasses.
//
// This class is inheritanted from the base class of SolutionAlgorithm
// which was created by fmk (Frank).


#include <RitzAlgo.h>
#include <AnalysisModel.h>
#include <RitzAnalysis.h>
#include <RitzIntegrator.h>
#include <RitzSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Timer.h>

#define RitzALGORITHM_TAGS 3

RitzAlgo::RitzAlgo()
:RitzAlgorithm(RitzALGORITHM_TAGS)
{
	// do nothing here.
}

RitzAlgo::~RitzAlgo()
{
	// do nothing here.
}

int 
RitzAlgo::solveCurrentStep(int numRitz)
{
	AnalysisModel *theModel = this->getAnalysisModelPtr();
	RitzSOE *theSOE = this->getRitzSOEptr();
	RitzIntegrator *theIntegrator = this->getRitzIntegratorPtr();

	if ((theModel == 0) || (theIntegrator == 0) || (theSOE == 0)) {
		opserr << "WARNING RitzAlgo::solverCurrentStep() - ";
		opserr << "setLinks() has not been called. \n";
		return -1;
	}

	if (theIntegrator->formK() < 0) {
		opserr << "WARNING RitzAlgo::solverCurrentStep() - ";
		opserr << "the Integrator failed in formK().\n";
		return -2;
	}

	if (theIntegrator->formM() < 0) {
		opserr << "WARNING RitzAlgo::solverCurrentStep() - ";
		opserr << "the Integrator failed in formK().\n";
		return -3;
	}
	
	if (theIntegrator->formB() < 0) {
		opserr << "WARNING RitzAlgo::solverCurrentStep() - ";
		opserr << "the Integrator failed in formB().\n";
		return -5;
	}
	
	/*if (theSOE->solve(numRitz) < 0) {
		opserr << "Warning RitzAlgo::solveCurrentStep() - ";
		opserr << "the EigenSOE failed in solve().\n";
		return -4;
	}

	// now set the ritzvalues and ritzvectors in the model
	// Commented while the real ritz algor is developed and set the ritzvectors
	/*theModel->setNumEigenvectors(numRitz);
	Vector theRitzvalues(numRitz);
	for (int i=1; i<=numRitz; i++) {
		theRitzvalues[i-1] = theSOE->getRitzvalue(i);
		theModel->setEigenvector(i, theSOE->getRitzvector(i));
	}    
	theModel->setEigenvalues(theRitzvalues);
	*/
	return 0;
}

int 
RitzAlgo::sendSelf(int cTag, Channel &theChannel)
{
	return 0;
}

int 
RitzAlgo::recvSelf(int cTag, Channel &theChannel,
						FEM_ObjectBroker &theBroker)
{
	return 0;
}

void 
RitzAlgo::Print(OPS_Stream &s, int flag)
{
	s << "\t Ritz Algorithm \n";
}


