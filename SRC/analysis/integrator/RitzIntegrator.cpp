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
// $Date: 2005/12/19 22:43:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/RitzIntegrator.cpp,v $

// Written: Jun Peng
// Created: Wed Jan 27, 1999
// Revision: A
//
// Description: This file contains the class definition of RitzIntegrator.
// RitzIntegrator is an algorithmic class for setting up the finite element 
// equations for a eigen problem analysis.
//
// This class is inheritanted from the base class of Integrator which was
// created by fmk (Frank).


#include <RitzIntegrator.h>
#include <FE_Element.h>
#include <AnalysisModel.h>
#include <RitzSOE.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <DOF_GrpIter.h>

#define RitzINTEGRATOR_TAGS        2

RitzIntegrator::RitzIntegrator()
:Integrator(RitzINTEGRATOR_TAGS),
theSOE(0), theAnalysisModel(0)
{

}

RitzIntegrator::~RitzIntegrator()
{

}

void
RitzIntegrator::setLinks(AnalysisModel &theModel, RitzSOE &theSysOE)
{
	theAnalysisModel = &theModel;
	theSOE = &theSysOE;
}

int 
RitzIntegrator::formEleTangent(FE_Element *theEle)
{
	if (flagK == 0)
		return this->formEleTangK(theEle);
	else
		return this->formEleTangM(theEle);
}

int 
RitzIntegrator::formNodTangent(DOF_Group *theDof)
{
	return this->formNodTangM(theDof);
}

int 
RitzIntegrator::formEleResidual(FE_Element *theEle)
{
	// only elements residual needed
	theEle->zeroResidual();
	//theEle->addRInertia();
	theEle->addRIncInertiaToResidual();
	return 0;
}

int 
RitzIntegrator::formNodUnbalance(DOF_Group *theDof)
{
	// only nodes unbalance need be added
	theDof->zeroUnbalance();
	//theDof->addPInertia();
	theDof->addPIncInertiaToUnbalance();
	return 0;
}

int 
RitzIntegrator::formElementResidual()
{
	// loop through the FE_Elements and add the residual
	FE_Element *elePtr;
	
	int res = 0;    
	
	FE_EleIter &theEles2 = theAnalysisModel->getFEs();    
	while((elePtr = theEles2()) != 0) {
	//      opserr << "ELEPTR " << elePtr->getResidual(this);
	
	if (theSOE->addB(elePtr->getResidual(this),elePtr->getID()) <0) {
		opserr << "WARNING RitzIntegrator::formElementResidual -";
		opserr << " failed in addB for ID " << elePtr->getID();
		res = -2;
		}
	}
	
	return res;	 
}

int 
RitzIntegrator::formNodalUnbalance(void)
{
	// loop through the DOF_Groups and add the unbalance
	DOF_GrpIter &theDOFs = theAnalysisModel->getDOFs();
	DOF_Group *dofPtr;
	int res = 0;
	
	while ((dofPtr = theDOFs()) != 0) { 
	//      opserr << "NODPTR: " << dofPtr->getUnbalance(this);
	
	if (theSOE->addB(dofPtr->getUnbalance(this),dofPtr->getID()) <0) {
		opserr << "WARNING RitzIntegrator::formNodalUnbalance -";
		opserr << " failed in addB for ID " << dofPtr->getID();
		res = -2;
		}
	}
	
	return res;
}

int 
RitzIntegrator::formUnbalance(void)
{
	if (theAnalysisModel == 0 || theSOE == 0) {
	opserr << "WARNING RitzIntegrator::formUnbalance -";
	opserr << " no AnalysisModel or LinearSOE has been set\n";
	return -1;
	}
	
	theSOE->zeroB();
	
	if (this->formElementResidual() < 0) {
	opserr << "WARNING RitzIntegrator::formUnbalance ";
	opserr << " - this->formElementResidual failed\n";
	return -1;
	}
	
	if (this->formNodalUnbalance() < 0) {
	opserr << "WARNING RitzIntegrator::formUnbalance ";
	opserr << " - this->formNodalUnbalance failed\n";
	return -2;
	}    
	
	return 0;
}


int 
RitzIntegrator::newStep()
{
	return 0;
}

int 
RitzIntegrator::getLastResponse(Vector &result, const ID &id)
{
	return 0;
}

int
RitzIntegrator::formK()
{
	if (theAnalysisModel == 0 || theSOE == 0) {
		opserr << "WARNING RitzIntegrator::formK -";
		opserr << " no AnalysisModel or EigenSOE has been set\n";
		return -1;
	}

	// the loops to form and add the tangents are broken into two for 
	// efficiency when performing parallel computations

	// loop through the FE_Elements getting them to form the tangent
	// FE_EleIter &theEles1 = theAnalysisModel->getFEs();
	FE_Element *elePtr;

	flagK = 0;

	theSOE->zeroA();

	//while((elePtr = theEles1()) != 0) 
	//  elePtr->formTangent(this);

	// loop through the FE_Elements getting them to add the tangent    
	int result = 0;
	FE_EleIter &theEles2 = theAnalysisModel->getFEs();    
	while((elePtr = theEles2()) != 0) {

		if (theSOE->addA(elePtr->getTangent(this), elePtr->getID()) < 0) {
			opserr << "WARNING RitzIntegrator::formK -";
			opserr << " failed in addA for ID " << elePtr->getID();	    
			result = -2;
		}
	}

	return result;    
}


int
RitzIntegrator::formM()
{
	if (theAnalysisModel == 0 || theSOE == 0) {
		opserr << "WARNING RitzIntegrator::formM -";
		opserr << " no AnalysisModel or EigenSOE has been set\n";
		return -1;
	}

	// the loops to form and add the tangents are broken into two for 
	// efficiency when performing parallel computations

	// loop through the FE_Elements getting them to form the tangent
	// FE_EleIter &theEles1 = theAnalysisModel->getFEs();
	FE_Element *elePtr;

	flagK = 1;
	theSOE->zeroM();

	// while((elePtr = theEles1()) != 0) 
	//     elePtr->formTangent(this);

	// loop through the FE_Elements getting them to add the tangent    
	int result = 0;
	FE_EleIter &theEles2 = theAnalysisModel->getFEs();    
	while((elePtr = theEles2()) != 0) {     
		if (theSOE->addM(elePtr->getTangent(this), elePtr->getID()) < 0) {
			opserr << "WARNING RitzIntegrator::formM -";
			opserr << " failed in addM for ID " << elePtr->getID();	    
			result = -2;
		}
	}

	DOF_Group *dofPtr;
	DOF_GrpIter &theDofs = theAnalysisModel->getDOFs();    
	while((dofPtr = theDofs()) != 0) {
		//   	dofPtr->formTangent(this);
		if (theSOE->addM(dofPtr->getTangent(this),dofPtr->getID()) < 0) {
			opserr << "WARNING RitzIntegrator::formM -";
			opserr << " failed in addM for ID " << dofPtr->getID();	    
			result = -3;
		}
	}

	return result;    
}

int
RitzIntegrator::formB()
{
	if (theAnalysisModel == 0 || theSOE == 0) {
		opserr << "WARNING RitzIntegrator::formB -";
		opserr << " no AnalysisModel or EigenSOE has been set\n";
		return -1;
	}

	// the loops to form and add the tangents are broken into two for 
	// efficiency when performing parallel computations

	// loop through the FE_Elements getting them to form the tangent
	// FE_EleIter &theEles1 = theAnalysisModel->getFEs();
	FE_Element *elePtr;

	theSOE->zeroB();

	// while((elePtr = theEles1()) != 0) 
	//     elePtr->formTangent(this);

	// loop through the FE_Elements getting them to add the tangent    
	int result = 0;
	FE_EleIter &theEles2 = theAnalysisModel->getFEs();    
	while((elePtr = theEles2()) != 0) {     
		if (theSOE->addB(elePtr->getResidual(this), elePtr->getID()) < 0) {
			opserr << "WARNING RitzIntegrator::formB -";
			opserr << " failed in addB for ID " << elePtr->getID();	    
			result = -2;
		}
	}

	DOF_Group *dofPtr;
	DOF_GrpIter &theDofs = theAnalysisModel->getDOFs();    
	while((dofPtr = theDofs()) != 0) {
		//   	dofPtr->formTangent(this);
		if (theSOE->addB(dofPtr->getM_Force(theSOE->getX(),1.0),dofPtr->getID()) < 0) {
			opserr << "WARNING RitzIntegrator::formB -";
			opserr << " failed in addB for ID " << dofPtr->getID();	    
			result = -3;
		}
	}

	return result;    
}

int 
RitzIntegrator::formEleTangK(FE_Element *theEle)
{
	theEle->zeroTangent();
	theEle->addKtToTang(1.0);
	return 0;
}

int 
RitzIntegrator::formEleTangM(FE_Element *theEle)
{
	theEle->zeroTangent();
	theEle->addMtoTang(1.0);
	return 0;
}

int 
RitzIntegrator::formNodTangM(DOF_Group *theDof)
{
	theDof->zeroTangent();
	theDof->addMtoTang(1.0);
	return 0;
}

int 
RitzIntegrator::update(const Vector &deltaU)
{
	return 0;
}

RitzSOE *
RitzIntegrator::getRitzSOEPtr() const
{
	return theSOE;
}

AnalysisModel *
RitzIntegrator::getAnalysisModelPtr() const
{
	return theAnalysisModel;
}

int 
RitzIntegrator::sendSelf(int commitTag, Channel &theChannel)
{
	return 0;
}

int 
RitzIntegrator::recvSelf(int commitTag, Channel &theChannel,
						  FEM_ObjectBroker &theBroker)
{
	return 0;
}

void 
RitzIntegrator::Print(OPS_Stream &s, int flag)
{
	s << "\t RitzIntegrator: \n";
}


