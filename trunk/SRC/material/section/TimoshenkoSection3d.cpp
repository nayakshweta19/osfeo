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
// $Date: 2007/11/30 23:34:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/TimoshenkoSection3d.cpp,v $
                                                                        
// Written: fmk
// Created: 04/04
//
// Description: This file contains the class implementation of TimoshenkoSection3d.

#include <stdlib.h>
#include <math.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <TimoshenkoSection3d.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <NDMaterial.h>

#include <string.h>

ID TimoshenkoSection3d::code(6);

// constructors:
TimoshenkoSection3d::TimoshenkoSection3d(int tag, int num, Fiber **fibers, double gj): 
  SectionForceDeformation(tag, SEC_TAG_TimoshenkoSection3d),
  numFibers(num), theMaterials(0), matData(0), GJ(gj), yh(0.0), zh(0.0),
  yBar(0.0), zBar(0.0), e(6), eCommit(6), s(0), ks(0)
{
  if (numFibers != 0) {
    theMaterials = new NDMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "TimoshenkoSection3d::TimoshenkoSection3d -- failed to allocate Material pointers\n";
      exit(-1);
    }

    matData = new double [numFibers*4];

    if (matData == 0) {
      opserr << "TimoshenkoSection3d::TimoshenkoSection3d -- failed to allocate double array for material data\n";
      exit(-1);
    }

    double Qz = 0.0;
    double Qy = 0.0;
    double a  = 0.0;
    double yHmin=0, zHmin=0, yHmax=0, zHmax=0;
	NDMaterial *theMat;
    for (int i = 0; i < numFibers; i++) {

      double yLoc, zLoc, Area, perpTheta;
	  fibers[i]->getFiberLocation(yLoc, zLoc);
	  Area = fibers[i]->getArea();
	  theMat = fibers[i]->getNDMaterial();
	  perpTheta = fibers[i]->getPerpTheta();

	  if (theMat == 0)
		opserr << "TimoshenkoSection3d::TimoshenkoSection3d -- failed to get fiber information" << endln;

      Qz += yLoc*Area;
      Qy += zLoc*Area;
      a  += Area;

      matData[i*4]   = -yLoc;
      matData[i*4+1] = zLoc;
      matData[i*4+2] = Area;
	  matData[i*4+3] = perpTheta;
	  
	  if (strcmp(theMat->getType(), "BeamFiber2d") == 0 )
        theMaterials[i] = theMat->getCopy("BeamFiber2d");

	  if (strcmp(theMat->getType(), "BeamFiber") == 0 )
		theMaterials[i] = theMat->getCopy("BeamFiber");

      if (theMaterials[i] == 0) 
	    opserr << "TimoshenkoSection3d::TimoshenkoSection3d -- failed to get copy of a Material\n";

	  if (yLoc < yHmin ) yHmin = yLoc;
	  if (zLoc < zHmin ) zHmin = zLoc;
	  if (yLoc > yHmax ) yHmax = yLoc;
	  if (zLoc > zHmax ) zHmax = zLoc;

    }
	zh   = yHmax - yHmin;
	yh   = zHmax - zHmin;
    yBar = -Qz/a;
    zBar = Qy/a;
  }

  s = new Vector(sData, 6);
  ks = new Matrix(kData, 6, 6);

  for (int i=0; i<6; i++)
	sData[i] = 0.0;
  for (int i=0; i<36; i++)
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_VY;
  code(4) = SECTION_RESPONSE_VZ;
  code(5) = SECTION_RESPONSE_T;
}

// constructor for blank object that recvSelf needs to be invoked upon
TimoshenkoSection3d::TimoshenkoSection3d():
  SectionForceDeformation(0, SEC_TAG_TimoshenkoSection3d),
  numFibers(0), theMaterials(0), matData(0), GJ(0), yh(0.0), zh(0.0),
  yBar(0.0), zBar(0.0), e(6), eCommit(6), s(0), ks(0)
{
  s = new Vector(sData, 6);
  ks = new Matrix(kData, 6, 6);

  for (int i=0; i<6; i++)
	sData[i] = 0.0;

  for (int i=0; i<36; i++)
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_VY;
  code(4) = SECTION_RESPONSE_VZ;
  code(5) = SECTION_RESPONSE_T;
}

// destructor:
TimoshenkoSection3d::~TimoshenkoSection3d()
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
// s = [P, Mz, My, Vy, Vz, T]'
// eps = [eps_11, eps_12, eps_13]'
// e   = [eps_a, kappa_z, kappa_y, gamma_y, gamma_z, phi]'
// a   = [1 -y z         0         0  0
//        0  0 0 sqrt(5/6)         0 -z
//        0  0 0         0 sqrt(5/6)  y]
int
TimoshenkoSection3d::setTrialSectionDeformation (const Vector &deforms)
{
  int res = 0;
  e = deforms;

  // Material strain vector
  static Vector eps(3);

  double y, z, A, theta; //

  //double root56 = sqrt(5.0/6.0);

  for (int i = 0; i < numFibers; i++) {

    y = matData[i*4] - yBar;
    z = matData[i*4+1] - zBar;
    A = matData[i*4+2];
	theta = matData[i*4+3];

    // determine material strain and set it

	eps(0) = e(0) - y*e(1) + z*e(2);
	eps(1) = e(3) - z*e(5);       // root56*e(3) - z*e(5);
	eps(2) = e(4) + y*e(5);       // root56*e(4) + y*e(5);

	if (strcmp(theMaterials[i]->getType(), "BeamFiber2d") == 0 )
	  res += theMaterials[i]->setTrialStrain(eps, theta);
	else
      res += theMaterials[i]->setTrialStrain(eps);

  }

  return res;
}

const Vector&
TimoshenkoSection3d::getSectionDeformation(void)
{
  return e;
}

double
TimoshenkoSection3d::getZh(void)
{
  double yHmin=0, zHmin=0, yHmax=0, zHmax=0;
  double yLoc, zLoc;
  for (int i = 0; i < numFibers; i++) {
	yLoc=matData[i*4];
 	zLoc=matData[i*4+1];

	if (yLoc < yHmin ) yHmin=yLoc;
	if (zLoc < zHmin ) zHmin=zLoc;
	if (yLoc > yHmax ) yHmax=yLoc;
	if (zLoc > zHmax ) zHmax=zLoc;
  }

  zh = yHmax - yHmin;
  yh = zHmax - zHmin;

  return zh;
}

double
TimoshenkoSection3d::getYh(void)
{
  double yHmin=0, zHmin=0, yHmax=0, zHmax=0;
  double yLoc, zLoc;
  for (int i = 0; i < numFibers; i++) {
	yLoc=matData[i*4];
 	zLoc=matData[i*4+1];

	if (yLoc < yHmin ) yHmin=yLoc;
	if (zLoc < zHmin ) zHmin=zLoc;
	if (yLoc > yHmax ) yHmax=yLoc;
	if (zLoc > zHmax ) zHmax=zLoc;
  }

  zh = yHmax - yHmin;
  yh = zHmax - zHmin;

  return yh;
}

double
TimoshenkoSection3d::getEIz(void)
{
  double EIz = 0.;
  double y, z, A;
  for (int i = 0; i < numFibers; i++) {
    y = matData[i*4] - yBar;
    z = matData[i*4+1] - zBar;
	A = matData[i*4+2];

	const Matrix &Dt = theMaterials[i]->getTangent();

	EIz += Dt(0,0) * A * pow(y,2.);
  }

  return EIz;
}

double
TimoshenkoSection3d::getEIy(void)
{
  double EIy = 0.;
  double y, z, A;
  for (int i = 0; i < numFibers; i++) {
    y = matData[i*4] - yBar;
    z = matData[i*4+1] - zBar;
	A = matData[i*4+2];

	const Matrix &Dt = theMaterials[i]->getTangent();

	EIy += Dt(0,0) * A * pow(z,2.);
  }

  return EIy;
}

double
TimoshenkoSection3d::getGAy(void)
{
  double GAy = 0.;
  double y, z, A;
  for (int i = 0; i < numFibers; i++) {
	//y =matData[i*4] - yBar;
 	//z =matData[i*4+1] - zBar;
	A =matData[i*4+2];

	const Matrix &Dt = theMaterials[i]->getTangent();

	GAy += Dt(1,1) * A;
  }

  return GAy;
}

double
TimoshenkoSection3d::getGAz(void)
{
  double GAz = 0.;
  double y, z, A;
  for (int i = 0; i < numFibers; i++) {
	//y =matData[i*4] - yBar;
 	//z =matData[i*4+1] - zBar;
	A =matData[i*4+2];

	const Matrix &Dt = theMaterials[i]->getTangent();

	GAz += Dt(2,2) * A;
  }

  return GAz;
}

// Compute section tangent stiffness, ks, from material tangent, Dt,
// by integrating over section area
// ks = int_A a'*Dt*a dA
// a   = [1 -y z         0         0  0
//        0  0 0 sqrt(5/6)         0 -z
//        0  0 0         0 sqrt(5/6)  y]
const Matrix&
TimoshenkoSection3d::getSectionTangent(void)
{
  for (int i=0; i<36; i++)
	kData[i] = 0.0;
 
  double y, z, w, theta; //
  double y2, z2, yz;
  
  double d00, d01, d02;
  double d10, d11, d12;
  double d20, d21, d22;
  
  double tmp;
  
  //double five6  = 5./6.;
  //double root56 = sqrt(five6);
  
  for (int i = 0; i < numFibers; i++) {
    
    y = matData[i*4] - yBar;
    z = matData[i*4+1] - zBar;
    w = matData[i*4+2];

    y2 = y*y;
    z2 = z*z;
    yz = y*z;

	Matrix Dt;  // it is a dangerous declaration
    if (strcmp(theMaterials[i]->getType(), "BeamFiber2d") == 0 ) {
	  theta = matData[i*4+3];
	  Dt = theMaterials[i]->getTangent(theta);
	} else
      Dt = theMaterials[i]->getTangent();

    d00 = Dt(0,0)*w; d01 = Dt(0,1)*w; d02 = Dt(0,2)*w;
    d10 = Dt(1,0)*w; d11 = Dt(1,1)*w; d12 = Dt(1,2)*w;
    d20 = Dt(2,0)*w; d21 = Dt(2,1)*w; d22 = Dt(2,2)*w;
    
    // Bending terms
    kData[0] += d00;      //(0,0)
    kData[7] += y2*d00;   //(1,1)
    kData[14] += z2*d00;  //(2,2)
    tmp = -y*d00;
    kData[1] += tmp;      //(0,1)
    kData[6] += tmp;      //(1,0)
    tmp = z*d00;
    kData[2] += tmp;      //(0,2)
    kData[12] += tmp;     //(2,0)
    tmp = -yz*d00;
    kData[8] += tmp;      //(1,2)
    kData[13] += tmp;     //(2,1)
    
    // Shear terms
    kData[21] += d11; //(3,3)
    kData[22] += d12; //(3,4)
    kData[27] += d21; //(4,3)
    kData[28] += d22; //(4,4)
    
    // Torsion term
	kData[35] += z2*d11 - yz*(d12+d21) + y2*d22; //(5,5)

    // Bending-torsion coupling terms
    tmp = -z*d01 + y*d02;
    kData[5] += tmp;    //(0,5)
    kData[11] -= y*tmp; //(1,5)
    kData[17] += z*tmp; //(2,5)
    tmp = -z*d10 + y*d20;
    kData[30] += tmp;   //(5,0)
    kData[31] -= y*tmp; //(5,1)
    kData[32] += z*tmp; //(5,2)
    
    // Hit tangent terms with root56
    //d01 *= root56; d02 *= root56;
    //d10 *= root56; d11 *= root56; d12 *= root56;
    //d20 *= root56; d21 *= root56; d22 *= root56;
    
    // Bending-shear coupling terms
    kData[3] += d01;    //(0,3)
    kData[4] += d02;    //(0,4)
    kData[9] -= y*d01;  //(1,3)
    kData[10] -= y*d02; //(1,4)
    kData[15] += z*d01; //(2,3)
    kData[16] += z*d02; //(2,4)
    kData[18] += d10;   //(3,0)
    kData[24] += d20;   //(4,0)
    kData[19] -= y*d10; //(3,1)
    kData[25] -= y*d20; //(4,1)
    kData[20] += z*d10; //(3,2)
    kData[26] += z*d20; //(4,2)
    
    // Torsion-shear coupling terms
    y2 =  y*d22;
    z2 = -z*d11;
    kData[33] +=  z2 + y*d21; //(5,3)
    kData[34] += -z*d12 + y2; //(5,4)
    kData[23] +=  z2 + y*d12; //(3,5)
    kData[29] += -z*d21 + y2; //(4,5)
  }

  // invert the stiffness matrix needed, for z-axial symmetric section
  if (kData[14] == 0.0)
    kData[14] += 1.0e-8;

  return *ks;
}

// Compute section stress resultant, s, from material stress, sig,
// by integrating over section area
// s = int_A a'*sig dA
// s = [P, Mz, My, Vy, Vz, T]'
// a = [1 -y z         0         0  0
//      0  0 0 sqrt(5/6)         0 -z
//      0  0 0         0 sqrt(5/6)  y]
const Vector&
TimoshenkoSection3d::getStressResultant(void)
{
  for (int i=0; i<6; i++)
	sData[i] = 0.0;
  
  double y, z, w, theta;
  double sig0, sig1, sig2;
  
  double root56 = sqrt(5./6.);

  for (int i = 0; i < numFibers; i++) {
    
    y = matData[i*4] - yBar;
    z = matData[i*4+1] - zBar;
    w = matData[i*4+2];

	Vector sig;  // it is a dangerous declaration
	if (strcmp(theMaterials[i]->getType(), "BeamFiber2d") == 0 ) {
	  theta = matData[i*4+3];
	  sig = theMaterials[i]->getStress(theta);
	} else
	  sig = theMaterials[i]->getStress();
    
    sig0 = sig(0)*w;
    sig1 = sig(1)*w;
    sig2 = sig(2)*w;
    
    sData[0] += sig0;
    sData[1] -= y*sig0;
    sData[2] += z*sig0;
    sData[3] += sig1;   //sig1;
    sData[4] += sig2;   //sig2;
    sData[5] += -z*sig1 + y*sig2;
  }

  return *s;
}

SectionForceDeformation*
TimoshenkoSection3d::getCopy(void)
{
  TimoshenkoSection3d *theCopy = new TimoshenkoSection3d ();
  theCopy->setTag(this->getTag());

  theCopy->numFibers = numFibers;

  if (numFibers != 0) {
    theCopy->theMaterials = new NDMaterial *[numFibers];

    if (theCopy->theMaterials == 0) {
      opserr << "TimoshenkoSection3d::TimoshenkoSection3d -- failed to allocate Material pointers\n";
      exit(-1);			    
    }

    theCopy->matData = new double [numFibers*4];

    if (theCopy->matData == 0) {
      opserr << "TimoshenkoSection3d::TimoshenkoSection3d -- failed to allocate double array for material data\n";
      exit(-1);
    }
	
    for (int i = 0; i < numFibers; i++) {
      theCopy->matData[i*4]   = matData[i*4];
      theCopy->matData[i*4+1] = matData[i*4+1];
      theCopy->matData[i*4+2] = matData[i*4+2];
	  theCopy->matData[i*4+3] = matData[i*4+3];
      theCopy->theMaterials[i] = theMaterials[i]->getCopy();

      if (theCopy->theMaterials[i] == 0) {
	opserr << "TimoshenkoSection3d::getCopy -- failed to get copy of a Material\n";
	exit(-1);
      }
    }    
  }

  theCopy->eCommit = eCommit;
  theCopy->e = e;
  theCopy->yBar = yBar;
  theCopy->zBar = zBar;

  for (int i=0; i<36; i++)
    theCopy->kData[i] = kData[i];
  for (int i=0; i<6; i++)
	theCopy->sData[i] = sData[i];

  return theCopy;
}

const ID&
TimoshenkoSection3d::getType ()
{
  return code;
}

int
TimoshenkoSection3d::getOrder () const
{
  return 6;
}

int
TimoshenkoSection3d::commitState(void)
{
  int err = 0;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->commitState();

  eCommit = e;

  return err;
}

int
TimoshenkoSection3d::revertToLastCommit(void)
{
  int err = 0;

  // Last committed section deformations
  e = eCommit;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->revertToLastCommit();

  return err;
}

int
TimoshenkoSection3d::revertToStart(void)
{
  // revert the fibers to start    
  int err = 0;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->revertToStart();

  return err;
}

int
TimoshenkoSection3d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // create an id to send objects tag and numFibers, 
  //     size 3 so no conflict with matData below if just 1 fiber
  static ID data(3);
  data(0) = this->getTag();
  data(1) = numFibers;
  int dbTag = this->getDbTag();
  res += theChannel.sendID(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "TimoshenkoSection2d::sendSelf - failed to send ID data\n";
    return res;
  }    

  if (numFibers != 0) {
    
    // create an id containing classTag and dbTag for each material & send it
    ID materialData(2*numFibers);
    for (int i=0; i<numFibers; i++) {
      NDMaterial *theMat = theMaterials[i];
      materialData(2*i) = theMat->getClassTag();
      int matDbTag = theMat->getDbTag();
      if (matDbTag == 0) {
	matDbTag = theChannel.getDbTag();
	if (matDbTag != 0)
	  theMat->setDbTag(matDbTag);
      }
      materialData(2*i+1) = matDbTag;
    }    
    
    res += theChannel.sendID(dbTag, commitTag, materialData);
    if (res < 0) {
     opserr << "TimoshenkoSection2d::sendSelf - failed to send material data\n";
     return res;
    }    

    // send the fiber data, i.e. area and loc
    Vector fiberData(matData, 4*numFibers);
    res += theChannel.sendVector(dbTag, commitTag, fiberData);
    if (res < 0) {
     opserr << "TimoshenkoSection2d::sendSelf - failed to send material data\n";
     return res;
    }    

    // now invoke send(0 on all the materials
    for (int j=0; j<numFibers; j++)
      theMaterials[j]->sendSelf(commitTag, theChannel);
  }

  return res;
}

int
TimoshenkoSection3d::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static ID data(3);
  
  int dbTag = this->getDbTag();
  res += theChannel.recvID(dbTag, commitTag, data);

  if (res < 0) {
   opserr << "TimoshenkoSection2d::sendSelf - failed to recv ID data\n";
   return res;
  } 
   
  this->setTag(data(0));

  // recv data about materials objects, classTag and dbTag
  if (data(1) != 0) {
    ID materialData(2*data(1));
    res += theChannel.recvID(dbTag, commitTag, materialData);
    if (res < 0) {
     opserr << "TimoshenkoSection2d::sendSelf - failed to send material data\n";
     return res;
    }    

    // if current arrays not of correct size, release old and resize
    if (theMaterials == 0 || numFibers != data(1)) {
      // delete old stuff if outa date
      if (theMaterials != 0) {
	for (int i=0; i<numFibers; i++)
	  delete theMaterials[i];
	delete [] theMaterials;
	if (matData != 0)
	  delete [] matData;
	matData = 0;
	theMaterials = 0;
      }

      // create memory to hold material pointers and fiber data
      numFibers = data(1);
      if (numFibers != 0) {

	theMaterials = new NDMaterial *[numFibers];
	
	if (theMaterials == 0) {
	  opserr << "TimoshenkoSection2d::recvSelf -- failed to allocate Material pointers\n";
	  exit(-1);
	}

	for (int j=0; j<numFibers; j++)
	  theMaterials[j] = 0;
	
	matData = new double [numFibers*4];

	if (matData == 0) {
	  opserr << "TimoshenkoSection2d::recvSelf  -- failed to allocate double array for material data\n";
	  exit(-1);
	}
      }
    }

    Vector fiberData(matData, 4*numFibers);
    res += theChannel.recvVector(dbTag, commitTag, fiberData);
    if (res < 0) {
     opserr << "TimoshenkoSection2d::sendSelf - failed to send material data\n";
     return res;
    }    
    
    int i;
    for (i=0; i<numFibers; i++) {
      int classTag = materialData(2*i);
      dbTag = materialData(2*i+1);

      // if material pointed to is blank or not of correct type, 
      // release old and create a new one
      if (theMaterials[i] == 0)
	theMaterials[i] = theBroker.getNewNDMaterial(classTag);
      else if (theMaterials[i]->getClassTag() != classTag) {
	delete theMaterials[i];
	theMaterials[i] = theBroker.getNewNDMaterial(classTag);      
      }

      if (theMaterials[i] == 0) {
	opserr << "TimoshenkoSection2d::recvSelf -- failed to allocate double array for material data\n";
	exit(-1);
      }

      theMaterials[i]->setDbTag(dbTag);
      res += theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
    }

    double Qz = 0.0;
    double Qy = 0.0;
    double A  = 0.0;
    double yLoc, zLoc, Area, perpTheta;

    // Recompute centroid
    for (i = 0; i < numFibers; i++) {
      yLoc = -matData[4*i];
      zLoc =  matData[4*i+1];
      Area =  matData[4*i+2];
	  perpTheta = matData[4*i+3];
      A  += Area;
      Qz += yLoc*Area;
      Qy += zLoc*Area;
    }
    
    yBar = -Qz/A;
    zBar = Qy/A;
  }    

  return res;
}

void
TimoshenkoSection3d::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {
    int loc = 0;
    for (int i = 0; i < numFibers; i++) {
      s << -matData[loc] << " "  << matData[loc+1] << " "  << matData[loc+2] << " " << matData[loc+3] << " ";
      s << theMaterials[i]->getStress() << " "  << theMaterials[i]->getStrain() << endln;
      loc += 4;
    } 
  } else {
    s << "\nTimoshenkoSection3d, tag: " << this->getTag() << endln;
    s << "\tSection code: " << code;
    s << "\tNumber of Fibers: " << numFibers << endln;
    s << "\tCentroid: (" << -yBar << ", " << zBar << ')' << endln;
    
    if (flag == 1) {
      int loc = 0;
      for (int i = 0; i < numFibers; i++) {
	s << "\nLocation (y, z) = (" << -matData[loc] << ", " << matData[loc+1] << ")";
	s << "\nArea = " << matData[loc+2] << "\nperpTheta = " << matData[loc+3] << endln;
	loc+= 4;
	theMaterials[i]->Print(s, flag);
      }
    }
  }
}

Response*
TimoshenkoSection3d::setResponse(const char **argv, int argc,
				 OPS_Stream &output)
{
  const ID &type = this->getType();
  int typeSize = this->getOrder();

  Response *theResponse =0;

  output.tag("SectionOutput");
  output.attr("secType", this->getClassType());
  output.attr("secTag", this->getTag());

  // deformations
  if (strcmp(argv[0],"deformations") == 0 || strcmp(argv[0],"deformation") == 0) {
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","kappaZ");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","eps");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","gammaY");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","kappaY");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","gammaZ");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","theta");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    theResponse =  new MaterialResponse(this, 1, this->getSectionDeformation());
  
  // forces
  } else if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0) {
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","Mz");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","P");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","Vy");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","My");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","Vz");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","T");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    theResponse =  new MaterialResponse(this, 2, this->getStressResultant());
  
  // force and deformation
  } else if (strcmp(argv[0],"forceAndDeformation") == 0) { 
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","kappaZ");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","eps");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","gammaY");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","kappaY");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","gammaZ");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","theta");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    for (int j=0; j<typeSize; j++) {
      int code = type(j);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","Mz");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","P");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","Vy");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","My");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","Vz");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","T");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }

    theResponse =  new MaterialResponse(this, 4, Vector(2*this->getOrder()));
  
  } else if (strcmp(argv[0],"fiberData") == 0) {
    int numData = numFibers*5;
    for (int j = 0; j < numFibers; j++) {
      output.tag("FiberOutput");
      output.attr("yLoc", -matData[3*j]);
      output.attr("zLoc", matData[3*j+1]);
      output.attr("area", matData[3*j+2]);    
      output.tag("ResponseType","yCoord");
      output.tag("ResponseType","zCoord");
      output.tag("ResponseType","area");
      output.tag("ResponseType","stress");
      output.tag("ResponseType","strain");
      output.endTag();
    }
    Vector theResponseData(numData);
    theResponse = new MaterialResponse(this, 5, theResponseData);
  }

  else {
    if (argc > 2 || strcmp(argv[0],"fiber") == 0) {
      
      int key = numFibers;
      int passarg = 2;
      
      if (argc <= 3)	{  // fiber number was input directly
	key = atoi(argv[1]);
      } else if (argc > 4) {         // find fiber closest to coord. with mat tag
	int matTag = atoi(argv[3]);
	double yCoord = atof(argv[1]);
	double zCoord = atof(argv[2]);
	double closestDist = 0.0;
	double ySearch, zSearch, dy, dz;
	double distance;
	int j;
	
	// Find first fiber with specified material tag
	for (j = 0; j < numFibers; j++) {
	  if (matTag == theMaterials[j]->getTag()) {
	    ySearch = -matData[3*j];
	    zSearch =  matData[3*j+1];
	    dy = ySearch-yCoord;
	    dz = zSearch-zCoord;
	    closestDist = sqrt(dy*dy + dz*dz);
	    key = j;
	    break;
	  }
	}
	
	// Search the remaining fibers
	for ( ; j < numFibers; j++) {
	  if (matTag == theMaterials[j]->getTag()) {
	    ySearch = -matData[3*j];
	    zSearch =  matData[3*j+1];
	    dy = ySearch-yCoord;
	    dz = zSearch-zCoord;
	    distance = sqrt(dy*dy + dz*dz);
	    if (distance < closestDist) {
	      closestDist = distance;
	      key = j;
	    }
	  }
	}
	passarg = 4;
      }
      
      else {                  // fiber near-to coordinate specified
	double yCoord = atof(argv[1]);
	double zCoord = atof(argv[2]);
	double closestDist;
	double ySearch, zSearch, dy, dz;
	double distance;
	ySearch = -matData[0];
	zSearch =  matData[1];
	dy = ySearch-yCoord;
	dz = zSearch-zCoord;
	closestDist = sqrt(dy*dy + dz*dz);
	key = 0;
	for (int j = 1; j < numFibers; j++) {
	  ySearch = -matData[3*j];
	  zSearch =  matData[3*j+1];
	  dy = ySearch-yCoord;
	  dz = zSearch-zCoord;
	  distance = sqrt(dy*dy + dz*dz);
	  if (distance < closestDist) {
	    closestDist = distance;
	    key = j;
	  }
	}
	passarg = 3;
      }
      
      if (key < numFibers && key >= 0) {
	output.tag("FiberOutput");
	output.attr("yLoc",-matData[3*key]);
	output.attr("zLoc",matData[3*key+1]);
	output.attr("area",matData[3*key+2]);
	
	theResponse =  theMaterials[key]->setResponse(&argv[passarg], argc-passarg, output);
	
	output.endTag();
      }
    }
  }

  output.endTag();
  return theResponse;
}

int 
TimoshenkoSection3d::getResponse(int responseID, Information &sectInfo)
{
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  if (responseID == 5) {
    double yLoc, zLoc, A, theta;
    Vector stress, strain;
	int numData = 9*numFibers;
    Vector data(numData);
    int count = 0;
    for (int j = 0; j < numFibers; j++) {
	  yLoc  = -matData[4*j];
	  zLoc  = matData[4*j+1];
      A     = matData[4*j+2];
	  if (strcmp(theMaterials[j]->getType(), "BeamFiber2d") == 0 ) {
        theta = matData[4*j+3];
	    stress = theMaterials[j]->getStress(theta);
        strain = theMaterials[j]->getStrain(theta);
	  } else {
		stress = theMaterials[j]->getStress();
		strain = theMaterials[j]->getStrain();
	  }
      data(count) = yLoc;
	  data(count+1) = zLoc;
	  data(count+2) = A;
      data(count+3) = stress(0); 
	  data(count+4) = stress(1); 
	  data(count+5) = stress(2);
	  data(count+6) = strain(0);
	  data(count+7) = strain(1);
	  data(count+8) = strain(2);
      count += 9;
    }
    return sectInfo.setVector(data);	
  } else {
    return SectionForceDeformation::getResponse(responseID, sectInfo);
  }
}

int
TimoshenkoSection3d::setParameter(const char **argv, int argc, 
				  Parameter &param)
{
  if (argc < 1)
    return -1;

  int result = -1;

  // A material parameter
  if (strstr(argv[0],"material") != 0) {

    if (argc < 3)
      return -1;
    
    // Get the tag of the material
    int paramMatTag = atoi(argv[1]);
    
    // Loop over fibers to find the right material(s)
    int ok = 0;
    for (int i = 0; i < numFibers; i++)
      if (paramMatTag == theMaterials[i]->getTag()) {
	ok = theMaterials[i]->setParameter(&argv[2], argc-2, param);
	if (ok != -1)
	  result = ok;
      }
  } 

  return result;
}
