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
// $Date: 2010/02/04 19:10:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/TimoshenkoSection2d.cpp,v $
                                                                        
// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition for 
// TimoshenkoSection2d.h. TimoshenkoSection2d provides the abstraction of a 
// rectangular section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.
// The fiber stresses are the 11, 12, and 13 components of stress, from
// which all six beam stress resultants are obtained.

#include <stdlib.h>
#include <math.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <Fiber.h>
#include <classTags.h>
#include <TimoshenkoSection2d.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <NDMaterial.h>

#include <string.h>

ID TimoshenkoSection2d::code(3);

// constructors:
TimoshenkoSection2d::TimoshenkoSection2d(int tag, int num, Fiber **fibers):
  SectionForceDeformation(tag, SEC_TAG_TimoshenkoSection2d),
	  numFibers(num), theMaterials(0), matData(0),
	  yBar(0.0), zBar(0.0), e(3), eCommit(3), s(0), ks(0)
{
  if (numFibers != 0) {
  
    theMaterials = new NDMaterial*[numFibers];
    if (theMaterials == 0) {
 	   opserr << "TimoshenkoSection2d::TimoshenkoSection2d -- failed to allocate Material pointers\n";
 	   exit(-1);
    }
    
    matData = new double [numFibers*3];
    if (matData == 0) {
 	   opserr << "TimoshenkoSection2d::TimoshenkoSection2d -- failed to allocate double array for material data\n";
 	   exit(-1);
    }
    
    double Qz = 0.0;
    double Qy = 0.0;
    double a  = 0.0;
    
	NDMaterial *theMat;
    for (int i = 0; i < numFibers; i++) {
	  double yLoc, zLoc, Area;
	  fibers[i]->getFiberLocation(yLoc, zLoc);
	  Area = fibers[i]->getArea();
	  theMat = fibers[i]->getNDMaterial();
	  if (theMat == 0) {
	    opserr << "TimoshenkoSection3d::TimoshenkoSection3d -- failed to get fiber information" << endln;
	  }
      
 	  Qz += yLoc*Area;
 	  Qy += zLoc*Area;
 	  a  += Area;
      
 	  matData[i*3] = -yLoc;
 	  matData[i*3+1] = zLoc;
 	  matData[i*3+2] = Area;
      
      theMaterials[i] = theMat->getCopy("BeamFiber2d") ; // theMat.getCopy("TimoshenkoFiber");
      if (theMaterials[i] == 0)
         opserr << "TimoshenkoSection2d::TimoshenkoSection2d -- failed to get copy of beam fiber" << endln;
    }
    yBar = -Qz/a;
    zBar = Qy/a;
  } 
  
  s = new Vector(sData, 3);
  ks = new Matrix(kData, 3, 3);

  for (int i=0; i<9; i++)
	kData[i] = 0.0;
  for (int i=0; i<3; i++)
	sData[i] = 0.0;
  
  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_VY;

}

// constructor for blank object that recvSelf needs to be invoked upon
TimoshenkoSection2d::TimoshenkoSection2d():
  SectionForceDeformation(0, SEC_TAG_TimoshenkoSection2d),
	  numFibers(0), theMaterials(0), matData(0),
	  yBar(0.0), zBar(0.0), e(3), eCommit(3), s(0), ks(0)
{
  s = new Vector(sData, 3);
  ks = new Matrix(kData, 3, 3);
  
  for (int i=0; i<9; i++)
	kData[i] = 0.0;
  for (int i=0; i<3; i++)
	sData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_VY;
}

// destructor:
TimoshenkoSection2d::~TimoshenkoSection2d()
{
  if (theMaterials != 0) {
  	for (int i = 0; i < numFibers; i++)
  	  if (theMaterials[i] != 0)
  	    delete theMaterials[i];
  
  	delete [] theMaterials;
  }
  
  if (matData != 0)
  	delete [] matData;
  
  if (s != 0)
  	delete s;
  
  if (ks != 0)
	delete ks;
}

// Compute fiber strain, eps, from section deformation, e, using the
// linear kinematic relationship, eps = a*e, where
// eps = [eps_11, eps_12 ]'
// e   = [eps_a, kappa_z, gamma_y ]'
// a   = [1 -y  0   
//        0  0  1 ]
int TimoshenkoSection2d::setTrialSectionDeformation (const Vector &deforms)
{
  int res = 0;
  
  e = deforms;
  
  // Material strain vector
  static Vector eps(2);
  
  double y, z;
  
  double root56 = sqrt(5.0/6.0);

  for (int i = 0; i < numFibers; i++) {
    y = matData[i*3] - yBar;
    z = 0.0 - zBar; //matData[i*3+1] - zBar;
      
    eps(0) = e(0) - y*e(1);
    eps(1) = e(2);
    
    res += theMaterials[i]->setTrialStrain(eps);
  }

  return res;
}

const Vector&
TimoshenkoSection2d::getSectionDeformation(void)
{
  return e;
}

// Compute section tangent stiffness, ks, from material tangent, Dt,
// by integrating over section area
// eps = [eps_11, eps_12 ]'
// e   = [eps_a, kappa_z, gamma_y ]'
// ks = int_A a'*Dt*a dA
// a   = [1 -y  0   
//        0  0  1 ]
const Matrix&
TimoshenkoSection2d::getSectionTangent(void)
{
  for (int i=0; i<9; i++)
	kData[i] = 0.0;
 
  double y, z, w;
  double y2, z2, yz;
  
  double d00, d01;
  double d10, d11;
  
  double tmp;
  
  double five6  = 5./6.;
  double root56 = sqrt(five6);
  
  for (int i = 0; i < numFibers; i++) {
    
    y = matData[i*3] - yBar;
    z = 0.0 - zBar; //matData[i*3+1] - zBar;
    w = matData[i*3+2];

    y2 = y*y;
    z2 = z*z;
    yz = y*z;
    
    const Matrix &Dt = theMaterials[i]->getTangent();

	d00 = Dt(0,0)*w; d01 = Dt(0,1)*w; 
	d10 = Dt(1,0)*w; d11 = Dt(1,1)*w; 
    
    // Bending terms
    kData[0] += d00;  //0,0
    kData[4] += y2*d00; //1,1
    kData[8] += d11; //2,2

    tmp = -y*d00;
    kData[1] += tmp; // 0,1
    kData[3] += tmp; //1,0
    // Shear terms
    kData[2] += d01; //0,2
    kData[6] += d10; //2,0
    tmp = -y*d01;
    kData[5] += tmp; //1,2
	tmp = -y*d10;
    kData[7] += tmp; //2,1
    
  }

  return *ks;
}

// Compute section stress resultant, s, from material stress, sig,
// by integrating over section area
// s = int_A a'*sig dA
// s = [P, Mz, Vy ]'
// a   = [1 -y  0   
//        0  0  1 ]
const Vector&
TimoshenkoSection2d::getStressResultant(void)
{
  for (int i=0; i<3; i++)
	sData[i] = 0.0;
  
  double y, z, w;
  double sig0, sig1;
  
  double root56 = sqrt(5./6.);

  for (int i = 0; i < numFibers; i++) {
    
    y = matData[i*3];
    z = 0.0; //matData[i*3+1];
    w = matData[i*3+2];
    
    const Vector &sig = theMaterials[i]->getStress();
    
    sig0 = sig(0)*w;
    sig1 = sig(1)*w;
    
    sData[0] += sig0;
    sData[1] -= y*sig0;
    sData[2] += sig1;

  }

  return *s;
}

SectionForceDeformation*
TimoshenkoSection2d::getCopy(void)
{

  TimoshenkoSection2d *theCopy = new TimoshenkoSection2d();
  theCopy->setTag(this->getTag());
  theCopy->numFibers = numFibers;

  if (numFibers != 0) {
    theCopy->theMaterials = new NDMaterial *[numFibers];

    if (theCopy->theMaterials == 0) {
      opserr << "TimoshenkoSection2d::TimoshenkoSection2d -- failed to allocate Material pointers\n";
      exit(-1);			    
    }

    theCopy->matData = new double [numFibers*3];

    if (theCopy->matData == 0) {
      opserr << "TimoshenkoSection2d::TimoshenkoSection2d -- failed to allocate double array for material data\n";
      exit(-1);
    }
			    
    
    for (int i = 0; i < numFibers; i++) {
      theCopy->matData[i*3] = matData[i*3];
      theCopy->matData[i*3+1] = matData[i*3+1];
      theCopy->matData[i*3+2] = matData[i*3+2];
      theCopy->theMaterials[i] = theMaterials[i]->getCopy("BeamFiber2d");

      if (theCopy->theMaterials[i] == 0) {
	opserr << "TimoshenkoSection2d::getCopy -- failed to get copy of a Material\n";
	exit(-1);
      }
    }    
  }
  theCopy->eCommit = eCommit;
  theCopy->e = e;
  theCopy->yBar = yBar;
  theCopy->zBar = zBar;
  
  for (int i=0; i<9; i++)
	theCopy->kData[i] = kData[i];
  for (int i=0; i<3; i++)
	theCopy->sData[i] = sData[i];

  return theCopy;
}

const ID&
TimoshenkoSection2d::getType(void)
{
  return code;
}

int
TimoshenkoSection2d::getOrder(void) const
{
  return 3;
}

int
TimoshenkoSection2d::commitState(void)
{
  int err = 0;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->commitState();
  
  eCommit = e;

  return err;
}

int
TimoshenkoSection2d::revertToLastCommit(void)
{
  int err = 0;

  // Last committed section deformations
  e = eCommit;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->revertToLastCommit();
  
  return err;
}

int
TimoshenkoSection2d::revertToStart(void)
{
  int err = 0;
  
  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->revertToStart();
  
  return err;
}

int
TimoshenkoSection2d::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
TimoshenkoSection2d::recvSelf(int commitTag, Channel &theChannel,
		     FEM_ObjectBroker &theBroker)
{
  return -1;
}

void
TimoshenkoSection2d::Print(OPS_Stream &s, int flag)
{
  s << "\nTimoshenkoSection2d, tag: " << this->getTag() << endln;
  s << "\tFiber Material, tag: " << theMaterials[0]->getTag() << endln;
  theMaterials[0]->Print(s, flag);
}
