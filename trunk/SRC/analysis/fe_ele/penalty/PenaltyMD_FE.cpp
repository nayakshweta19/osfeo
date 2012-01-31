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
// $Date: 2008-07-03 17:58:49 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/analysis/fe_ele/penalty/PenaltyMD_FE.cpp,v $
                                                                        
                                                                        
#include <PenaltyMD_FE.h>
#include <stdlib.h>

#include <Element.h>
#include <Domain.h>
#include <Node.h>
#include <DOF_Group.h>
#include <Integrator.h>
#include <IncrementalIntegrator.h>
#include <Subdomain.h>
#include <AnalysisModel.h>
#include <Matrix.h>
#include <Vector.h>
#include <Node.h>
#include <MD_Constraint.h>
#include <DOF_Group.h>
#include <iostream>

PenaltyMD_FE::PenaltyMD_FE(int tag, Domain &theDomain, 
			   MP_Constraint &TheMP, double Alpha)
			   :FE_Element(tag, 2,(TheMP.getConstrainedDOFs()).Size()+
			   (TheMP.getRetainedDOFs()).Size()),
			   theMP(&TheMP), theConstrainedNode(0) , theRetainedNode(0),
			   tang(0), resid(0), alpha(Alpha)
{
    
    int size;
    size = 9;

    tang = new Matrix(size,size);
    resid = new Vector(size);

    if (tang == 0 || resid == 0 || tang->noCols() != size || resid->Size() != size) {
	opserr << "FATAL PenaltyMP_FE::PenaltyMP_FE() - out of memory\n";
	exit(-1);
    }
	    
    theRetainedNode = theDomain.getNode(theMP->getNodeRetained());    
    theConstrainedNode = theDomain.getNode(theMP->getNodeConstrained());

    if (theRetainedNode == 0 || theConstrainedNode == 0) {
	opserr << "FATAL PenaltyMP_FE::PenaltyMP_FE() - Constrained or Retained";
	opserr << " Node does not exist in Domain\n";
	opserr << theMP->getNodeRetained() << " " << theMP->getNodeConstrained() << endln;
	exit(-1);
    }	

    // set up the dof groups tags
    DOF_Group *dofGrpPtr = 0;
 
    dofGrpPtr = theConstrainedNode->getDOF_GroupPtr();
    if (dofGrpPtr != 0) 
	myDOF_Groups(0) = dofGrpPtr->getTag();	        
    else
	opserr << "WARNING PenaltyMP_FE::PenaltyMP_FE() - node no Group yet?\n"; 
    
    
    if (theMP->isTimeVarying() == false) {
	this->determineTangent();
    }
}

PenaltyMD_FE::~PenaltyMD_FE()
{
    if (tang != 0)
	delete tang;
    if (resid != 0)
	delete resid;
}    

// void setID(int index, int value);
//	Method to set the correMPonding index of the ID to value.
int
PenaltyMD_FE::setID(void)
{
    int result = 0;

    // first determine the IDs in myID for those DOFs marked
    // as constrained DOFs, this is obtained from the DOF_Group
    // associated with the constrained node
    DOF_Group *theConstrainedNodesDOFs = theConstrainedNode->getDOF_GroupPtr();
    if (theConstrainedNodesDOFs == 0) {
	opserr << "WARNING PenaltyMP_FE::setID(void)";
	opserr << " - no DOF_Group with Constrained Node\n";
	return -2;
    }    

    const ID &constrainedDOFs = theMP->getConstrainedDOFs();
    const ID &theConstrainedNodesID = theConstrainedNodesDOFs->getID();    
    
    int size1 = 9;
    for (int i=0; i<size1; i++) {
	int constrained = i;
	if (constrained < 0 || 
	    constrained >= theConstrainedNode->getNumberDOF()) {
	    
	    opserr << "WARNING PenaltyMP_FE::setID(void) - unknown DOF ";
	    opserr << constrained << " at Node\n";
	    myID(i) = -1; // modify so nothing will be added to equations
	    result = -3;
	}    	
	else {
	    if (constrained >= theConstrainedNodesID.Size()) {
		opserr << "WARNING PenaltyMP_FE::setID(void) - ";
		opserr << " Nodes DOF_Group too small\n";
		myID(i) = -1; // modify so nothing will be added to equations
		result = -4;
	    }
	    else
		myID(i) = theConstrainedNodesID(constrained);
	}
    }
    
    myDOF_Groups(0) = theConstrainedNodesDOFs->getTag();

    return result;
}

const Matrix &
PenaltyMD_FE::getTangent(Integrator *theNewIntegrator)
{
    if (theMP->isTimeVarying() == true)
	this->determineTangent(); 

    //IncrementalIntegrator *theincrIntegr = (IncrementalIntegrator *) theNewIntegrator;
    //if (theincrIntegr->getStatus() == 1)
    //   cout<<"it worked"<<endl;   
    return *tang;
}

const Vector &
PenaltyMD_FE::getResidual(Integrator *theNewIntegrator)
{
    // zero residual, CD = 0
    return *resid;
}



const Vector &
PenaltyMD_FE::getTangForce(const Vector &disp, double fact)
{
    // does nothing , zero residual for CD = 0
    return *resid;
}

const Vector &
PenaltyMD_FE::getK_Force(const Vector &disp, double fact)
{
 //opserr << "WARNING PenaltyMP_FE::getK_Force() - not yet implemented\n";
 (*resid).Zero();
 return *resid;
}

const Vector &
PenaltyMD_FE::getC_Force(const Vector &disp, double fact)
{
 //opserr << "WARNING PenaltyMP_FE::getC_Force() - not yet implemented\n";
 (*resid).Zero();
 return *resid;
}

const Vector &
PenaltyMD_FE::getM_Force(const Vector &disp, double fact)
{
 //opserr << "WARNING PenaltyMP_FE::getM_Force() - not yet implemented\n";
 (*resid).Zero();
 return *resid;
}

void  
PenaltyMD_FE::determineTangent(void)
{
    const Matrix &constraint = theMP->getConstraint();
    int noRows = constraint.noRows();
    int noCols = constraint.noCols();
    
    for (int i=0; i<noRows; i++)
	for (int j=0; j<noCols; j++)
	    (*tang)(i,j) = alpha * constraint(i,j);
    
}


