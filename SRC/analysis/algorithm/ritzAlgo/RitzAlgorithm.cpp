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
// $Date: 2000/09/15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/eigenAlgo/RitzAlgorithm.cpp,v $


// File: ~/analysis/algorithm/eigenAlgo/RitzAlgorithm.C
//
// Written: Jun Peng
// Created: Wed Jan 27, 1999
// Revision: A
//
// Description: This file contains the class definition of RitzAlgorithm.
// RitzAlgorithm is a class which performs a eigen solution algorithm
// to solve the equations. 
//
// This class is inheritanted from the base class of SolutionAlgorithm
// which was created by fmk (Frank).


#include <RitzAlgorithm.h>
#include <AnalysisModel.h>
#include <RitzIntegrator.h>
#include <RitzSOE.h>

RitzAlgorithm::RitzAlgorithm(int classTag)
:SolutionAlgorithm(classTag),
theModel(0), theIntegrator(0), theSOE(0)
{
	// need do nothing here.
}


RitzAlgorithm::~RitzAlgorithm()
{
	// do nothing here.
}

void 
RitzAlgorithm::setLinks(AnalysisModel &theNewModel,
						RitzIntegrator &theNewIntegrator,
						RitzSOE &theNewSOE)
{
	theModel = &theNewModel;
	theIntegrator = &theNewIntegrator;
	theSOE = &theNewSOE;
}

AnalysisModel *
RitzAlgorithm::getAnalysisModelPtr() const
{
	return theModel;
}

RitzIntegrator * 
RitzAlgorithm::getRitzIntegratorPtr() const
{
	return theIntegrator;
}

RitzSOE *
RitzAlgorithm::getRitzSOEptr() const
{
	return theSOE;
}

