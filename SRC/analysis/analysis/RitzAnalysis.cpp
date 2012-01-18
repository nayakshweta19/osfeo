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
// $Date: 2005/08/31 17:39:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/RitzAnalysis.cpp,v $


// File: ~/analysis/analysis/eigenAnalysis/RitzAnalysis.C
//
// Written: Jun Peng
// Created: Wed Jan 27, 1999
// Revision: A
//
// Description: This file contains the class definition of RitzAnalysis.
// RitzAnalysis is a subclass of Analysis, it is used to perform the 
// eigen vlaue analysis on the FE_Model.
//
// This class is inheritanted from the base class of Analysis
// which was created by fmk (Frank).


#include <RitzAnalysis.h>
#include <RitzAlgorithm.h>
#include <AnalysisModel.h>
#include <RitzSOE.h>
#include <DOF_Numberer.h>
#include <ConstraintHandler.h>
#include <RitzIntegrator.h>
#include <Domain.h>
#include <Timer.h>
#include <FE_Element.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <DOF_GrpIter.h>
#include <Matrix.h>
#include <ID.h>
#include <Graph.h>


RitzAnalysis::RitzAnalysis(Domain &the_Domain,
							 ConstraintHandler &theHandler,
							 DOF_Numberer &theNumberer,
							 AnalysisModel &theModel,
							 RitzAlgorithm &theAlgo,
							 RitzSOE &theSOE,
							 RitzIntegrator &theRitzIntegrator)
							 :Analysis(the_Domain), theConstraintHandler(&theHandler),
							 theDOF_Numberer(&theNumberer), theAnalysisModel(&theModel),
							 theAlgorithm(&theAlgo), theSOE(&theSOE),
							 theIntegrator(&theRitzIntegrator), domainStamp(0)
{
	// first set up the links needed by the elements in the aggregation.
	theAnalysisModel->setLinks(the_Domain, *theConstraintHandler);
	theConstraintHandler->setLinks(the_Domain, theModel, theRitzIntegrator);
	theDOF_Numberer->setLinks(theModel);
	theIntegrator->setLinks(theModel, theSOE);
	theAlgorithm->setLinks(theModel, theRitzIntegrator, theSOE);
}


RitzAnalysis::~RitzAnalysis()
{
	// do nothing now.
	this->clearAll();
}

void
RitzAnalysis::clearAll(void)
{
	// invoke the destructor on all the objects in the aggregation
	delete theAnalysisModel;
	delete theConstraintHandler;
	delete theDOF_Numberer;
	delete theIntegrator;
	delete theAlgorithm;
	delete theSOE;
}    



int 
RitzAnalysis::analyze(int numRitz)
{
	int result = 0;
	Domain *the_Domain = this->getDomainPtr();

	// check for change in Domain since last step. As a change can
	// occur in a commit() in a domaindecomp with load balancing
	// this must now be inside the loop
	int stamp = the_Domain->hasDomainChanged();
	if (stamp != domainStamp) {
		domainStamp = stamp;
		result = this->domainChanged();
		if (result < 0) {
			opserr << "RitzAnalysis::analyze() - domainChanged failed\n";
			return -1;
		}	
	}

	result = theIntegrator->newStep();
	if (result < 0) {
		opserr << "RitzAnalysis::analyze() - integrator failed\n";
		return -2;
	}

	result = theAlgorithm->solveCurrentStep(numRitz);
	if (result < 0) {
		opserr << "RitzAnalysis::analyze() - algorithm failed\n";
		return -3;
	}

	return 0;
}

int 
RitzAnalysis::domainChanged()
{
	theAnalysisModel->clearAll();    
	theConstraintHandler->clearAll();      
	theConstraintHandler->handle();

	theDOF_Numberer->numberDOF();
	theConstraintHandler->doneNumberingDOF();

	Graph &theGraph = theAnalysisModel->getDOFGraph();
	theSOE->setSize(theGraph);

	theIntegrator->domainChanged();
	theAlgorithm->domainChanged();

	return 0;
}

int 
RitzAnalysis::setAlgorithm(RitzAlgorithm &theAlgo)
{
	opserr << "RitzAnalysis::setAlgorithm() - does nothing yet\n";
	return 0;
}

int 
RitzAnalysis::setIntegrator(RitzIntegrator &theIntegrator)
{
	opserr << "RitzAnalysis::setIntegrator() - does nothing yet\n";    
	return 0;
}

int 
RitzAnalysis::setRitzSOE(RitzSOE &theSOE)
{
	opserr << "RitzAnalysis::setRitzSOE() - does nothing yet\n";    
	return 0;
}

ConstraintHandler *
RitzAnalysis::getConstraintHandlerPtr() const
{
	return theConstraintHandler;
}

DOF_Numberer *
RitzAnalysis::getDOF_NumbererPtr() const
{
	return theDOF_Numberer;
}

AnalysisModel *
RitzAnalysis::getAnalysisModelPtr() const
{
	return theAnalysisModel;
}

RitzAlgorithm * 
RitzAnalysis::getRitzAlgorithm() const
{
	return theAlgorithm;
}

RitzSOE *
RitzAnalysis::getRitzSOE() const
{
	return theSOE;
}

RitzIntegrator * 
RitzAnalysis::getRitzIntegrator() const
{
	return theIntegrator;
}
