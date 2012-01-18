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
// $Date: 2003/02/14 23:00:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/RitzIntegrator.h,v $


// File: ~/analysis/integrator/eigenIntegrator/RitzIntegrator.h
//
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


#ifndef RitzIntegrator_h
#define RitzIntegrator_h

#include <Integrator.h>

class RitzSOE;
class AnalysisModel;
class FE_Element;
class DOF_Group;
class Vector;

class RitzIntegrator : public Integrator
{
public:
	RitzIntegrator();
	virtual ~RitzIntegrator();

	virtual void setLinks(AnalysisModel &theModel,
		RitzSOE &theSOE);

	// methods to form the M and K matrices.
	virtual int formK();
	virtual int formM();
	// methods to form the B vector.
	virtual int formB();

	// methods to instruct the FE_Element and DOF_Group objects
	// how to determine their contribution to M and K
	virtual int formEleTangK(FE_Element *theEle);
	virtual int formEleTangM(FE_Element *theEle);
	virtual int formNodTangM(DOF_Group *theDof);
	virtual int update(const Vector &deltaU);

	virtual int formEleTangent(FE_Element *theEle);
	virtual int formNodTangent(DOF_Group *theDof);
	virtual int formEleResidual(FE_Element *theEle);
	virtual int formNodUnbalance(DOF_Group *theDof);

	virtual int newStep(void);

	virtual int formElementResidual(void);
	virtual int formNodalUnbalance(void);
	virtual int formUnbalance(void);

	virtual int getLastResponse(Vector &result, const ID &id);

	virtual int sendSelf(int commitTag, Channel &theChannel);
	virtual int recvSelf(int commitTag, Channel &theChannel,
		FEM_ObjectBroker &theBroker);     
	virtual void Print(OPS_Stream &s, int flag = 0);

protected:
	virtual RitzSOE *getRitzSOEPtr() const;
	virtual AnalysisModel *getAnalysisModelPtr() const;

private:
	RitzSOE *theSOE;
	AnalysisModel *theAnalysisModel;
	int flagK;

};

#endif




