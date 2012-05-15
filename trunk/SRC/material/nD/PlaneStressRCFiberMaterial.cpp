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
// $Date: 2003-02-14 23:01:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlaneStressRCFiberMaterial.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition of PlaneStressRCFiberMaterial.
// The PlaneStressRCFiberMaterial class is a wrapper class that performs static
// condensation on a three-dimensional material model to give the 11 and 12
// stress components which can then be integrated over an area to model a
// shear flexible 2D beam.

#include <PlaneStressRCFiberMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>
#include <float.h>

Vector PlaneStressRCFiberMaterial::stress(2);
Matrix PlaneStressRCFiberMaterial::tangent(2,2);

#define ND_TAG_PlaneStressRCFiberMaterial 265891


PlaneStressRCFiberMaterial::PlaneStressRCFiberMaterial(void)
: NDMaterial(0, ND_TAG_PlaneStressRCFiberMaterial),
Tstrain22(0.0), Cstrain22(0.0), twoDtgLastCommit(3,3), strain(2), theMaterial(0)
{
	// Nothing to do
}

PlaneStressRCFiberMaterial::PlaneStressRCFiberMaterial(int tag, NDMaterial &theMat)
: NDMaterial(tag, ND_TAG_PlaneStressRCFiberMaterial),
Tstrain22(0.0), Cstrain22(0.0), twoDtgLastCommit(3,3), theMaterial(0), strain(2)
{
  // Get a copy of the material
  theMaterial = theMat.getCopy();
  
  if (theMaterial == 0) {
    opserr << "PlaneStressRCFiberMaterial::PlaneStressRCFiberMaterial -- failed to get copy of material\n";
    exit(-1);
  }
}

PlaneStressRCFiberMaterial::~PlaneStressRCFiberMaterial(void) 
{ 
  if (theMaterial != 0)
    delete theMaterial;
} 

NDMaterial*
PlaneStressRCFiberMaterial::getCopy(void) 
{
	PlaneStressRCFiberMaterial *theCopy =
		new PlaneStressRCFiberMaterial(this->getTag(), *theMaterial);

	theCopy->Tstrain22 = this->Tstrain22;
	theCopy->Cstrain22 = this->Cstrain22;
	theCopy->twoDtgLastCommit = this->twoDtgLastCommit;

	return theCopy;
}

NDMaterial* 
PlaneStressRCFiberMaterial::getCopy(const char *type)
{
	if (strcmp(type, "BeamFiber2d") == 0)
		return this->getCopy();
	else
		return 0;
}

int 
PlaneStressRCFiberMaterial::getOrder(void) const
{
	return 2;
}

const char*
PlaneStressRCFiberMaterial::getType(void) const 
{
	return "BeamFiber2d";
}

int 
PlaneStressRCFiberMaterial::commitState(void)
{
  Cstrain22 = Tstrain22;
  if (theMaterial->commitState() != 0)
	  return -1;
  twoDtgLastCommit = theMaterial->getTangent();
  
  return 0;
}

int 
PlaneStressRCFiberMaterial::revertToLastCommit(void)
{
  Tstrain22 = Cstrain22;
  if (theMaterial->revertToLastCommit() != 0)
	  return -1;

  return 0;
}

int
PlaneStressRCFiberMaterial::revertToStart()
{
  this->Tstrain22 = 0.0;
  this->Cstrain22 = 0.0;
  this->twoDtgLastCommit = theMaterial->getInitialTangent();

  return theMaterial->revertToStart();
}

double
PlaneStressRCFiberMaterial::getRho(void)
{
  return theMaterial->getRho();
}

//receive the strain
//2D Plane Stress Material strain order = 11, 22, 12
//PlaneStressRCFiberMaterial strain order = 11, 12, 22
//11, 22, 33, 12, 23, 31
//11, 12, 31, 22, 33, 23
int 
PlaneStressRCFiberMaterial::setTrialStrain(const Vector &strainFromElement)
{
  static const double tolerance = 1.0e-05;

  this->strain(0) = strainFromElement(0);
  this->strain(1) = strainFromElement(1);

  //newton loop to solve for out-of-plane strains
  double norm;
  static double condensedStress;
  static double strainIncrement;
  static Vector twoDstress(3);
  static Vector twoDstrain(3);
  static Matrix twoDtangent(3,3);
  static Vector twoDstressCopy(3); 
  static Matrix twoDtangentCopy(3,3);

  int i, j;
  int ii, jj;

  do 
  {  
    //set two dimensional strain
    twoDstrain(0) = this->strain(0);     //Tstrain11    eps_xx
    twoDstrain(1) = this->Tstrain22;     //Tstrain22;   eps_yy
    twoDstrain(2) = this->strain(1);     //Tstrain12;   gamma_xy -> eps_xy*2.0

	if (theMaterial->setTrialStrain(twoDstrain) < 0) {
      opserr << "PlaneStressRCFiberMaterial::setTrialStrain - setStrain failed in material with strain " << twoDstrain;
      return -1;
    }

	//two dimensional stress
    twoDstress = theMaterial->getStress();

	//two dimensional tangent 
	twoDtangent = theMaterial->getTangent();

    //Plane stress material strain order = 11, 22, 33, 12, 23, 31
    //BeamFiber 2d material strain order = 11, 12, 31, 22, 33, 23
    //swap matrix indices to sort out-of-plane components
    for (i=0; i<6; i++) {
      ii = this->indexMap(i);
      twoDstressCopy(ii) = twoDstress(i);
      for (j=0; j<6; j++) {
	jj = this->indexMap(j);
	twoDtangentCopy(ii,jj) = twoDtangent(i,j);
      }//end for j
    }//end for i

  } while (norm > tolerance);

	//this->tangent(0,0) = twoDtangent(0,0) + (twoDtangent(0,1)+twoDtangent(1,0))/dd22; 
	//this->tangent(0,1) = twoDtangent(0,2) + (twoDtangent(0,1)+twoDtangent(1,2))/dd22;
	//this->tangent(1,0) = twoDtangent(2,0) + (twoDtangent(2,1)+twoDtangent(1,0))/dd22;
	//this->tangent(1,1) = twoDtangent(2,2) + (twoDtangent(2,1)+twoDtangent(1,2))/dd22;
	
	//this->stress(0) = tangent(0,0)*strain(0)+tangent(0,1)*strain(1);
	//this->stress(1) = tangent(1,0)*strain(0)+tangent(1,1)*strain(1);

  return 0;
}

const Vector& 
PlaneStressRCFiberMaterial::getStrain(void)
{
  return this->strain;
}

const Vector&  
PlaneStressRCFiberMaterial::getStress()
{
  return this->stress;
}

const Matrix&  
PlaneStressRCFiberMaterial::getTangent()
{
  return this->tangent;
}

const Matrix&  
PlaneStressRCFiberMaterial::getInitialTangent()
{
  //static Matrix dd11(2,2);
  //static Matrix dd12(2,1);
  //static Matrix dd21(1,2);
  //static Matrix dd22(1,1);
  //static Matrix dd22invdd21(1,2);
  static Matrix twoDtangentCopy(3,3);

  const Matrix &twoDtangent = theMaterial->getInitialTangent();

  double dd22 = twoDtangent(1,1);

  this->tangent(0,0) = twoDtangent(0,0) + (twoDtangent(0,1)+twoDtangent(1,0))/dd22; 
  this->tangent(0,1) = twoDtangent(0,2) + (twoDtangent(0,1)+twoDtangent(1,2))/dd22;
  this->tangent(1,0) = twoDtangent(2,0) + (twoDtangent(2,1)+twoDtangent(1,0))/dd22;
  this->tangent(1,1) = twoDtangent(2,2) + (twoDtangent(2,1)+twoDtangent(1,2))/dd22;

  return this->tangent;
}

//2D material strain order        = 11, 22, 12
//PlaneStressRCFiberMaterial strain order = 11, 12, 22
int 
PlaneStressRCFiberMaterial::indexMap(int i)
{
  int ii;

  if (i == 1) 
	  ii = 2;
  else if (i == 2)
	  ii = 1;
  else 
	  ii = i;

  return ii;
}

void  
PlaneStressRCFiberMaterial::Print(OPS_Stream &s, int flag)
{
	s << "PlaneStressRCFiberMaterial, tag: " << this->getTag() << endln;
	s << "\tWrapped material: "<< theMaterial->getTag() << endln;

	theMaterial->Print(s, flag);
}

int 
PlaneStressRCFiberMaterial::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0;

  // put tag and associated materials class and database tags into an id and send it
  static ID idData(3);
  idData(0) = this->getTag();
  idData(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMaterial->setDbTag(matDbTag);
  }
  idData(2) = matDbTag;

  res = theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "PlaneStressRCFiberMaterial::sendSelf() - failed to send id data\n";
    return res;
  }

  // put the strains in a vector and send it
  static Vector vecData(1);
  vecData(0) = Cstrain22;

  res = theChannel.sendVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "PlaneStressRCFiberMaterial::sendSelf() - failed to send vector data\n";
    return res;
  }

  // now send the materials data
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "PlaneStressRCFiberMaterial::sendSelf() - failed to send vector material\n";

  return res;
}

int 
PlaneStressRCFiberMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  // recv an id contain the tag and associated materials class and db tags
  static ID idData(3);
  res = theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "PlaneStressRCFiberMaterial::sendSelf() - failed to send id data\n";
    return res;
  }

  this->setTag(idData(0));
  int matClassTag = idData(1);
  
  // if the associated material has not yet been created or is of the wrong type
  // create a new material for recvSelf later
  if (theMaterial == 0 || theMaterial->getClassTag() != matClassTag) {
    if (theMaterial != 0)
      delete theMaterial;
    theMaterial = theBroker.getNewNDMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "PlaneStressRCFiberMaterial::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMaterial->setDbTag(idData(2));

  // recv a vector containing strains and set the strains
  static Vector vecData(1);
  res = theChannel.recvVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "PlaneStressRCFiberMaterial::sendSelf() - failed to send vector data\n";
    return res;
  }

  Cstrain22 = vecData(0);

  Tstrain22 = Cstrain22;

  // now receive the materials data
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) 
    opserr << "PlaneStressRCFiberMaterial::sendSelf() - failed to send vector material\n";
  
  return res;
}
