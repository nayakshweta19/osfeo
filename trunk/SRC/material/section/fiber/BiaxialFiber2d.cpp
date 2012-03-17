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
// $Date: 2011/06/16 01:18:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/fiber/BiaxialFiber2d.h,v $


// File: ~/fiber/BiaxialFiber2d.h
//
// Written: Neallee
// Created: 2011
// Revision: 

// Description: This file contains the class definition for 
// BiaxialFiber2d.h. BiaxialFiber2d provides the abstraction of a
// biaxial fiber whose position is defined with only one coordinate.
// The BiaxialFiber2d is subjected to a stress state with 
// nonzero axial and shear stresses and corresponding strains.
//
// What: "@(#) BiaxialFiber2d.cpp, revA"

#include <stdlib.h>

#include <NDMaterial.h>
#include <BiaxialFiber2d.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ID.h>
#include <SectionForceDeformation.h>
#include <Information.h>
#include <Parameter.h>
#include <FiberResponse.h>

Matrix BiaxialFiber2d::ks(2,2); 
Vector BiaxialFiber2d::fs(2); 
ID BiaxialFiber2d::code(3);

// constructor:
BiaxialFiber2d::BiaxialFiber2d(int tag, 
                                 NDMaterial &theMat,
                                 double Area, double position):
                                 Fiber(tag, FIBER_TAG_Biaxial2d),
                                 theMaterial(0), area(Area), y(-position)
{
  theMaterial = theMat.getCopy("BeamFiber2d");  // get a copy of the MaterialModel
  
  if (theMaterial == 0) {
    opserr <<"BiaxialFiber2d::BiaxialFiber2d  -- failed to get copy of NDMaterial\n";
    exit(-1);
  }
  
  if (code(0) != SECTION_RESPONSE_P) {
    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
	code(2) = SECTION_RESPONSE_VY;
  }
}

// constructor for blank object that recvSelf needs to be invoked upon
BiaxialFiber2d::BiaxialFiber2d(): Fiber(0, FIBER_TAG_Biaxial2d),
                                    theMaterial(0), area(0), y(0.0)
{
  if (code(0) != SECTION_RESPONSE_P) {
    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
	code(2) = SECTION_RESPONSE_VY;
  }
}

// Destructor: 
BiaxialFiber2d::~BiaxialFiber2d ()
{
  if (theMaterial != 0)
    delete theMaterial;
}

int   
BiaxialFiber2d::setTrialFiberStrain(const Vector &vs)
{
  // Use the section kinematic matrix to get the fiber strain
  // eps = as * vs;
  Vector strain;
  strain(0) = vs(0) + y*vs(1); // fiber strain, axial
  strain(1) = vs(1);
  strain(2) = vs(2);
  
  return theMaterial->setTrialStrain(strain);
}

// get fiber stress resultants 
Vector &
BiaxialFiber2d::getFiberStressResultants (void) 
{
  // Use the section kinematic matrix to get the fiber 
  // stress resultant vector
  // fs = as^ * area * sigma;
  Vector df = theMaterial->getStress() * area;
  
  fs(0) = df(0);
  fs(1) = y * df(0);
  fs(2) = df(2);

  return fs;
}

// get contribution of fiber to section tangent stiffness
Matrix &
BiaxialFiber2d::getFiberTangentStiffContr(void) 
{
  // Use the section kinematic matrix to get the fiber 
  // tangent stiffness matrix
  // ks = (as^as) * area * Et;
  Matrix value = theMaterial->getTangent() * area;
  double value_as1 = value(0,0)*y;
  
  ks(0,0) = value(0,0);
  ks(0,1) = value_as1;
  ks(1,0) = value_as1;
  ks(1,1) = value_as1 * y;
  
  return ks;
}

Fiber*
BiaxialFiber2d::getCopy (void)
{
   // make a copy of the fiber 
   BiaxialFiber2d *theCopy = new BiaxialFiber2d (this->getTag(), 
                                                   *theMaterial, area, -y);
   return theCopy;
}  

int
BiaxialFiber2d::getOrder(void)
{
	return 3;
}

const ID&
BiaxialFiber2d::getType(void)
{
	return code;
}

int   
BiaxialFiber2d::commitState(void)
{
   return theMaterial->commitState();
}

int   
BiaxialFiber2d::revertToLastCommit(void)
{
   return theMaterial->revertToLastCommit();
}

int   
BiaxialFiber2d::revertToStart(void)
{
   return theMaterial->revertToStart();
}

int   
BiaxialFiber2d::sendSelf(int commitTag, Channel &theChannel)
{
  // 
  // store tag and material info in an ID and send it
  //
  int res = 0;
  
  int dbTag = this->getDbTag();
  
  static ID idData(3);
  
  idData(0) = this->getTag();
  idData(1) = theMaterial->getClassTag();
  
  int matDbTag = theMaterial->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      theMaterial->setDbTag(matDbTag);
  }
  
  idData(2) = matDbTag;
  
  res += theChannel.sendID(dbTag, commitTag, idData);
  if (res < 0) {
    opserr << "BiaxialFiber2d::sendSelf - failed to send ID data\n";
    return res;
  }    
  
  // 
  // store area and position data in a vector and send it
  //
  static Vector dData(2);
  
  dData(0) = area;
  dData(1) = y;
  
  res += theChannel.sendVector(dbTag, commitTag, dData);
  if (res < 0) {
    opserr << "BiaxialFiber2d::sendSelf - failed to send Vector data\n";
    return res;
  }    

  // now invoke sendSelf on the material
  res += theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr << "BiaxialFiber2d::sendSelf - failed to send NDMaterial\n";
      return res;
  }
    
  return res;
}

int
BiaxialFiber2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  // 
  // get tag and material info from an ID
  //
  
  int res = 0;
  
  int dbTag = this->getDbTag();
  
  static ID idData(3);
    
  res += theChannel.recvID(dbTag, commitTag, idData);
  if (res < 0) {
    opserr << "BiaxialFiber2d::rcvSelf - failed to receive ID data\n";
    return res;
  }    
  
  this->setTag(idData(0));

  // 
  // get area from a vector received from channel
  //
  
  static Vector dData(2);
  
  res += theChannel.recvVector(dbTag, commitTag, dData);
  if (res < 0) {
      opserr << "BiaxialFiber2d::recvSelf - failed to receive Vector data\n";
      return res;
  }
  
  area = dData(0);
  y = dData(1);

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
      opserr << "BiaxialFiber2d::recvSelf() - " <<
	  "failed to get a NDMaterial of type " << matClassTag << endln;
      return -1;
    }
  }
  
    // set the materials dbTag and invoke recvSelf on the material
  theMaterial->setDbTag(idData(2));
  
  // now invoke recvSelf on the material
  res += theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "BiaxialFiber2d::recvSelf() - the material failed in recvSelf()\n";
    return res;
  }    	
  
  return res;
}

void BiaxialFiber2d::Print(OPS_Stream &s, int flag)
{
  s << "\nbiaxialFiber2d, tag: " << this->getTag() << endln;
  s << "\tArea: " << area << endln; 
  s << "\tMatrix as: " << 1.0 << " " << y << endln; 
  s << "\tMaterial, tag: " << theMaterial->getTag() << endln;
}

Response*
BiaxialFiber2d::setResponse(const char **argv, int argc, OPS_Stream &s)
{
  if (argc == 0)
    return 0;
  
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
    return new FiberResponse(this, 1, Vector(2));
  
  else
    return theMaterial->setResponse(argv, argc, s);
}

int
BiaxialFiber2d::getResponse(int responseID, Information &fibInfo)
{
  switch(responseID) {
  case 1:
    return fibInfo.setVector(this->getFiberStressResultants());
    
  default:
    return -1;
  }
}

void 
BiaxialFiber2d::getFiberLocation(double &yLoc, double &zLoc)
{
  yLoc = -y;
  zLoc = 0.0;
}

int
BiaxialFiber2d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"A") == 0)
    return param.addObject(1, this);

  if (strcmp(argv[0],"y") == 0)
    return param.addObject(2, this);

  else
    return theMaterial->setParameter(argv, argc, param);
}

int
BiaxialFiber2d::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case 1:
    area = info.theDouble;
    return 0;
  case 2:
    y = -info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
BiaxialFiber2d::activateParameter(int parameterID)
{
  return -1;
}

const Vector&
BiaxialFiber2d::getFiberSensitivity(int gradNumber, bool cond)
{
  return Fiber::getFiberSensitivity(gradNumber, cond);
}

int 
BiaxialFiber2d::commitSensitivity(const Vector &dedh, int gradNumber,
				   int numGrads)
{
  return -1;
}

