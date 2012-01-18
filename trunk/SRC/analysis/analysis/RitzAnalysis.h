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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/RitzAnalysis.h,v $


// File: ~/analysis/analysis/eigenAnalysis/RitzAnalysis.h
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


#ifndef RitzAnalysis_h
#define RitzAnalysis_h

#include <Analysis.h>

class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class RitzAlgorithm;
class RitzIntegrator; 
class RitzSOE;

class RitzAnalysis : public Analysis
{
public:
	RitzAnalysis(Domain &theDomain,
		ConstraintHandler &theHandler,
		DOF_Numberer &theNumberer,
		AnalysisModel &theModel,
		RitzAlgorithm &theAlgo,
		RitzSOE &theSOE,
		RitzIntegrator &theIntegrator);

	virtual ~RitzAnalysis();

	virtual int analyze(int numRitz);
	void clearAll(void);	         
	virtual int domainChanged();

	virtual int setAlgorithm(RitzAlgorithm &theAlgo);
	virtual int setIntegrator(RitzIntegrator &theIntegrator);
	virtual int setRitzSOE(RitzSOE &theSOE);

protected:
	ConstraintHandler	*getConstraintHandlerPtr() const;
	DOF_Numberer	*getDOF_NumbererPtr() const;
	AnalysisModel	*getAnalysisModelPtr() const;
	RitzAlgorithm *getRitzAlgorithm() const;
	RitzSOE 		*getRitzSOE() const;
	RitzIntegrator *getRitzIntegrator() const;

private:
	ConstraintHandler 	*theConstraintHandler;
	DOF_Numberer	*theDOF_Numberer;
	AnalysisModel 	*theAnalysisModel;
	RitzAlgorithm	*theAlgorithm;
	RitzSOE			*theSOE;
	RitzIntegrator	*theIntegrator;
	int domainStamp;
}; 

#endif

