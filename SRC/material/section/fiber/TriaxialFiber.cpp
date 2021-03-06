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
                                                                        
// $Revision: 1.9 $
// $Date: 2007/02/02 01:18:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/fiber/TriaxialFiber.cpp,v $
                                                                        
                                                                        
// File: ~/fiber/TriaxialFiber.h
//
// Written: Neallee
// Created: 2011
// Revision: 
//
// Description: This file contains the implementation for the
// TriaxialFiber class. TriaxialFiber provides the abstraction of a
// uniaxial fiber that forms a fiber section for 3d frame elements.
// The TriaxialFiber is subjected to a stress state with
// only one nonzero axial stress and corresponding axial strain.
//
// What: "@(#) TriaxialFiber.C, revA"

#include <stdlib.h>
#include <stdio.h>

#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <TriaxialFiber.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ID.h>
#include <SectionForceDeformation.h>
#include <Information.h>
#include <FiberResponse.h>

Matrix TriaxialFiber::ks(3,3); 
Vector TriaxialFiber::fs(3); 
ID TriaxialFiber::code(6); 

// constructor:
TriaxialFiber::TriaxialFiber()
:Fiber(0, FIBER_TAG_Triaxial),
 theMaterial(0), area(0.0)
{
	if (code(0) != SECTION_RESPONSE_P) {
       code(0) = SECTION_RESPONSE_P;
       code(1) = SECTION_RESPONSE_MZ;
       code(2) = SECTION_RESPONSE_VY;
       code(3) = SECTION_RESPONSE_MY;
       code(4) = SECTION_RESPONSE_VZ;
       code(5) = SECTION_RESPONSE_T;
	}

   as[0] = 0.0;
   as[1] = 0.0;
}

TriaxialFiber::TriaxialFiber(int tag, 
                                 NDMaterial &theMat,
                                 double Area, const Vector &position)
:Fiber(tag, FIBER_TAG_Triaxial),
 theMaterial(0), area(Area)
{
	theMaterial = theMat.getCopy("BeamFiber");  // get a copy of the MaterialModel

	if (theMaterial == 0) {
	  opserr << "TriaxialFiber::TriaxialFiber -- failed to get copy of NDMaterial\n";
	  exit(-1);
	}
	
	if (code(0) != SECTION_RESPONSE_P) {
		code(0) = SECTION_RESPONSE_P;
		code(1) = SECTION_RESPONSE_MZ;
		code(2) = SECTION_RESPONSE_VY;
		code(3) = SECTION_RESPONSE_MY;
		code(4) = SECTION_RESPONSE_VZ;
		code(5) = SECTION_RESPONSE_T;
	}

	as[0] = -position(0);
	as[1] =  position(1);
}

// destructor:
TriaxialFiber::~TriaxialFiber ()
{
   if (theMaterial != 0)
      delete theMaterial;
}


int   
TriaxialFiber::setTrialFiberStrain(const Vector &vs)
{
  Vector strain;
  strain(0)= vs(0) + as[0]*vs(1) + as[1]*vs(2);

  if (theMaterial != 0)
      return theMaterial->setTrialStrain(strain);
  else {
    opserr << "TriaxialFiber::setTrialFiberStrain() - no material!\n";
    return -1; // in case fatal does not exit
  }
}



// get fiber stress resultants 
Vector &
TriaxialFiber::getFiberStressResultants (void)
{
    Vector df = theMaterial->getStress();
	df = df * area;

    // fs = as^ df;
    fs(0) = df(0);
    fs(1) = as[0]*df(0);
    fs(2) = as[1]*df(0);

    return fs;
}



// get contribution of fiber to section tangent stiffness
Matrix &
TriaxialFiber::getFiberTangentStiffContr(void) 
{
    // ks = (as^as) * area * Et;
    Matrix value = theMaterial->getTangent();
	value = value * area;

    double as1 = as[0];
    double as2 = as[1];
    double vas1 = as1*value(0,0);
    double vas2 = as2*value(1,1);
    double vas1as2 = vas1*as2;

    ks(0,0) = value(0,0);
    ks(0,1) = vas1;
    ks(0,2) = vas2;
    
    ks(1,0) = vas1;
    ks(1,1) = vas1*as1;
    ks(1,2) = vas1as2;
    
    ks(2,0) = vas2;
    ks(2,1) = vas1as2;
    ks(2,2) = vas2*as2;

    return ks;
}

Fiber*
TriaxialFiber::getCopy (void)
{
   // make a copy of the fiber 
   static Vector position(2);

   position(0) = -as[0];
   position(1) =  as[1];

   TriaxialFiber *theCopy = new TriaxialFiber (this->getTag(), 
                                                   *theMaterial, area, 
                                                   position);
   return theCopy;
}  

int
TriaxialFiber::getOrder(void)
{
	return 3;
}

const ID&
TriaxialFiber::getType(void)
{
	return code;
}

int   
TriaxialFiber::commitState(void)
{
   return theMaterial->commitState();
}


int   
TriaxialFiber::revertToLastCommit(void)
{
   return theMaterial->revertToLastCommit();
}

int   
TriaxialFiber::revertToStart(void)
{
   return theMaterial->revertToStart();
}


int   
TriaxialFiber::sendSelf(int commitTag, Channel &theChannel)
{
    // 
    // store tag and material info in an ID and send it
    //

    static ID idData(3);
    int dbTag = this->getDbTag();
    idData(0) = this->getTag();
    idData(1) = theMaterial->getClassTag();
    int matDbTag = theMaterial->getDbTag();
    if (matDbTag == 0) {
	matDbTag = theChannel.getDbTag();
	if (matDbTag != 0)
	    theMaterial->setDbTag(matDbTag);
    }
    idData(2) = matDbTag;
    
    if (theChannel.sendID(dbTag, commitTag, idData) < 0)  {
	opserr << "TriaxialFiber::sendSelf() -  failed to send ID data\n";
	return -1;
    }    
    
    // 
    // store area and position data in a vector and send it
    //
    
    static Vector dData(3);
    dData(0) = area;
    dData(1) = as[0];
    dData(2) = as[1];
    if (theChannel.sendVector(dbTag, commitTag, dData) < 0)  {
      opserr << "TriaxialFiber::sendSelf() -  failed to send Vector data\n";
      return -2;
    }    

    // now invoke sendSelf on the material
    if (theMaterial->sendSelf(commitTag, theChannel) < 0) {
      opserr << "TriaxialFiber::sendSelf() -  the material failed in sendSelf()\n";
      return -3;
    }    	
    
    return 0;
}


int   
TriaxialFiber::recvSelf(int commitTag, Channel &theChannel, 
			  FEM_ObjectBroker &theBroker)
{
    // 
    // get tag and material info from an ID
    //

    static ID idData(3);
    int dbTag = this->getDbTag();
    
    if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {
	opserr << "TriaxialFiber::recvSelf() -  failed to recv ID data\n";
	return -1;
    }    

    this->setTag(idData(0));

    // 
    // get area and position datafrom a vector
    //
    
    static Vector dData(3);
    if (theChannel.recvVector(dbTag, commitTag, dData) < 0)  {
      opserr << "TriaxialFiber::recvSelf() -  failed to recv Vector data\n";
	return -2;
    }        
    area = dData(0);
    as[0] = dData(1);
    as[1] = dData(2);

    //
    // now we do the material stuff
    //
    
    int matClassTag = idData(1);    
    
    // if we have a material, check it is of correct type
    if (theMaterial != 0) {
	if (matClassTag != theMaterial->getClassTag()) {
	    delete theMaterial;
	    theMaterial = 0;
	} 
    }

    // if no material we need to get one,
    // NOTE: not an else if in case deleted in if above
    if (theMaterial == 0) {
	theMaterial = theBroker.getNewNDMaterial(matClassTag);
	if (theMaterial == 0) {
	  opserr << "TriaxialFiber::recvSelf() - " << 
	    "failed to get a UniaxialMaterial of type "<< matClassTag << endln;
	    return -3;
	}
    }

    // set the materials dbTag and invoke recvSelf on the material
    theMaterial->setDbTag(idData(2));

    // now invoke recvSelf on the material
    if (theMaterial->recvSelf(commitTag, theChannel, theBroker) < 0) {
      opserr << "TriaxialFiber::recvSelf() -  the material failed in recvSelf()\n";
	return -4;
    }    	

    return 0;
}


void TriaxialFiber::Print(OPS_Stream &s, int flag)
{
    s << "\nTriaxialFiber, tag: " << this->getTag() << endln;
    s << "\tArea: " << area << endln; 
    s << "\tMatrix as: " << 1.0 << " " << as[0] << " " << as[1] << endln; 
    s << "\tMaterial, tag: " << theMaterial->getTag() << endln;
}

Response*
TriaxialFiber::setResponse(const char **argv, int argc, OPS_Stream &s)
{
	if (argc == 0)
		return 0;

	if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)		
		return new FiberResponse(this, 1, Vector(3));

	else
	  return theMaterial->setResponse(argv, argc, s);
}

int
TriaxialFiber::getResponse(int responseID, Information &fibInfo)
{
	switch(responseID) {
		case 1:
			return fibInfo.setVector(this->getFiberStressResultants());

		default:
			return -1;
	}
}

void 
TriaxialFiber::getFiberLocation(double &yLoc, double &zLoc)
{
	yLoc = -as[0];
	zLoc = as[1];
}


