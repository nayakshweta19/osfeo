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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlaneStressFiberMaterial.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition of PlaneStressFiberMaterial.
// The PlaneStressFiberMaterial class is a wrapper class that performs static
// condensation on a three-dimensional material model to give the 11 and 12
// stress components which can then be integrated over an area to model a
// shear flexible 2D beam.

#include <PlaneStressFiberMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>
#include <float.h>

Vector PlaneStressFiberMaterial::stress(2);
Matrix PlaneStressFiberMaterial::tangent(2,2);
//       0   1   2
// PS   11, 22, 12
// FB   11, 12, 22
int PlaneStressFiberMaterial::iMap[] = {0, 2, 1};

PlaneStressFiberMaterial::PlaneStressFiberMaterial(void)
: NDMaterial(0, ND_TAG_PlaneStressFiberMaterial),
Tstrain22(0.0), Cstrain22(0.0), twoDtgLastCommit(3,3), theMaterial(0), strain(2)
{
	// Nothing to do
}

PlaneStressFiberMaterial::PlaneStressFiberMaterial(int tag, NDMaterial &theMat)
: NDMaterial(tag, ND_TAG_PlaneStressFiberMaterial),
Tstrain22(0.0), Cstrain22(0.0), twoDtgLastCommit(3,3), theMaterial(0), strain(2)
{
  // Get a copy of the material
  theMaterial = theMat.getCopy();
  
  if (theMaterial == 0) {
    opserr << "PlaneStressFiberMaterial::PlaneStressFiberMaterial -- failed to get copy of material\n";
    exit(-1);
  }
}

PlaneStressFiberMaterial::~PlaneStressFiberMaterial(void) 
{ 
  if (theMaterial != 0)
    delete theMaterial;
} 

NDMaterial*
PlaneStressFiberMaterial::getCopy(void) 
{
	PlaneStressFiberMaterial *theCopy =
		new PlaneStressFiberMaterial(this->getTag(), *theMaterial);

	theCopy->Tstrain22 = this->Tstrain22;
	theCopy->Cstrain22 = this->Cstrain22;
	theCopy->twoDtgLastCommit = this->twoDtgLastCommit;

	return theCopy;
}

NDMaterial* 
PlaneStressFiberMaterial::getCopy(const char *type)
{
	if (strcmp(type, "BeamFiber2d") == 0)
		return this->getCopy();
	else
		return 0;
}

int 
PlaneStressFiberMaterial::getOrder(void) const
{
	return 2;
}

const char*
PlaneStressFiberMaterial::getType(void) const 
{
	return "BeamFiber2d";
}

int 
PlaneStressFiberMaterial::commitState(void)
{
  Cstrain22 = Tstrain22;
  if (theMaterial->commitState() != 0)
	  return -1;
  twoDtgLastCommit = theMaterial->getTangent();
  
  return 0;
}

int 
PlaneStressFiberMaterial::revertToLastCommit(void)
{
  Tstrain22 = Cstrain22;
  if (theMaterial->revertToLastCommit() != 0)
	  return -1;

  return 0;
}

int
PlaneStressFiberMaterial::revertToStart()
{
  this->Tstrain22 = 0.0;
  this->Cstrain22 = 0.0;
  this->twoDtgLastCommit = theMaterial->getInitialTangent();

  return theMaterial->revertToStart();
}

double
PlaneStressFiberMaterial::getRho(void)
{
  return theMaterial->getRho();
}

//receive the strain
//2D Plane Stress Material strain order = 11, 22, 12
//PlaneStressFiberMaterial strain order = 11, 12, 22
//11, 22, 33, 12, 23, 31
//11, 12, 31, 22, 33, 23
int 
PlaneStressFiberMaterial::setTrialStrain(const Vector &strainFromElement)
{

  this->strain(0) = strainFromElement(0);
  this->strain(1) = strainFromElement(1);

  //newton loop to solve for out-of-plane strains
  int i, j, ii, jj;
  double norm;
  static double factor = 10;
 
  static Vector condensedStress(1);
  static Vector strainIncrement(1);
  static Vector twoDstress(3);
  static Vector twoDstrain(3);
  static Matrix twoDtangent(3,3);
  static Vector twoDstressCopy(3); 
  static Matrix twoDtangentCopy(3,3);
  static Matrix dd22(1,1);

  static Vector twoDstrainTrial(3);
  static Vector twoDstrainToDo(3);

  //double dW;                    // section strain energy (work) norm 
  int maxSubdivisions = 4;
  int numSubdivide    = 1;
  bool converged      = false;
  maxIters            = 10;
  tol                 = 1.0e-5;

  //set two dimensional strain
  
  twoDstrain(0) = this->strain(0); //Tstrain11
  twoDstrain(1) = this->Tstrain22; //Tstrain22;
  twoDstrain(2) = this->strain(1); //Tstrain33;
  //twoDstrain(3) = this->Tgamma12; 
  //twoDstrain(4) = this->Tgamma23;
  //twoDstrain(5) = this->Tgamma31;

  twoDstrainToDo = twoDstrain;
  twoDstrainTrial = twoDstrainToDo;

  while (converged == false && numSubdivide <= maxSubdivisions) {
	
	// try regular newton (if l==0), or
	// initial tangent on first iteration then regular newton (if l==1), or 
	// initial tangent iterations (if l==2)

	for (int l=0; l<3; l++) {

  int numIters = maxIters;

  if (l == 1) 
	numIters = 10*maxIters; // allow 10 times more iterations for initial tangent

  for (int cnt=0; cnt <numIters; cnt++) {

    if (theMaterial->setTrialStrain(twoDstrainTrial) < 0) {
      opserr << "PlaneStressRCFiberMaterial::setTrialStrain - setStrain failed in material with strain " << twoDstrain;
      return -1;
    }

	if (l == 0)	{

  // regular newton
  twoDtangent = theMaterial->getTangent();

	} else if (l == 2) {

  // newton with initial tangent in first iteration
  // otherwise regular newton
  if (cnt == 0) {
	twoDtangent = twoDtgLastCommit;
  } else {
	twoDtangent = theMaterial->getTangent();
  }

	} else {

  //  newton with initial tangent
  twoDtangent = twoDtgLastCommit;

	}

    // two dimensional stress
    twoDstress = theMaterial->getStress();
    
	// PlaneStressRC NDmaterial strain order   = 11, 22, 12;  ##, 33, 23, 31
	// PlaneStressRCFiberMaterial strain order = 11, 12, 22;  ##, 31, 33, 23
	// swap matrix indices to sort out-of-plane components 
    for (i=0; i<3; i++) {
    
      ii = iMap[i];
    
      twoDstressCopy(ii) = twoDstress(i);
    
      for (j=0; j<3; j++) {
    
    jj = iMap[j];
    
    twoDtangentCopy(ii,jj) = twoDtangent(i,j);
    
      }//end for j
    }//end for i
    
    //out of plane stress and tangents
    condensedStress(0) = twoDstressCopy(2);
    dd22(0,0) = twoDtangentCopy(2,2);
    
    //set norm
    norm = condensedStress.Norm();
    
    //condensation 
    dd22.Solve(condensedStress, strainIncrement);
    
    //update out of plane strains
    //this->Tgamma31  -= strainIncrement(0);
    this->Tstrain22 -= strainIncrement(0);
    //this->Tstrain33 -= strainIncrement(2);
    //this->Tgamma23  -= strainIncrement(3);
    
    if (norm <= tol) {
	  converged = true;
	  //twoDstrainToDo -= twoDstrainTrial;

	  // break out of cnt & l loops
	  cnt = numIters+1;
	  l   = 4;
    } else {
	  // for next iteration
	  
	  twoDstrainTrial(1) = this->Tstrain22;

	  // if we have failed to converge for all of our newton schemes
	  // - reduce step size by the factor specified
	  
      if ((cnt == (numIters-1)) && (l == 2)) {
	    twoDstrainTrial /= factor;
	    numSubdivide++;
	  }
    } // test norm <? tol

  }  // for (cnt=0; cnt<numIters; cnt++)

    } // for (int l=0; l<2; l++)

  } // (converged == false && numSubdivide <= maxSubdivisions)
  
    // if fail to converge we return an error flag & print an error message

  if (converged == false) {
    opserr << "WARNING - PlaneStressFiberMaterial::setTrialStrain - failed to get compatible ";
    opserr << "inter-fiber forces & deformations FiberMaterial: " << endln;
    opserr << this->getTag() << "(norm: << " << norm << ")" << endln;
    return -1;
  }

  return 0;
}

const Vector& 
PlaneStressFiberMaterial::getStrain(void)
{
  return this->strain;
}

const Vector&  
PlaneStressFiberMaterial::getStress()
{
  //two dimensional stress
  const Vector &twoDstress = theMaterial->getStress();
  static Vector twoDstressCopy(3);

  //NDmaterial strain order        = 11, 22, 33, 12, 23, 31  
  //PlaneStressFiberMaterial strain order = 11, 12, 31, 22, 33, 23

  //swap matrix indices to sort out-of-plane components 
  int i, ii;
  
  for (i=0; i<3; i++) {

    ii = iMap[i];
    
    twoDstressCopy(ii) = twoDstress(i);

  }//end for i

  for (i=0; i<2; i++) 
    this->stress(i) = twoDstressCopy(i);

  return this->stress;
}

const Matrix&  
PlaneStressFiberMaterial::getTangent()
{
  static Matrix dd11(2,2);
  static Matrix dd12(2,1);
  static Matrix dd21(1,2);
  static Matrix dd22(1,1);
  static Matrix dd22invdd21(1,2);
  static Matrix twoDtangentCopy(3,3);

  const Matrix &twoDtangent = theMaterial->getTangent();

  //swap matrix indices to sort out-of-plane components 
  int i, j, ii, jj;
  for (i=0; i<3; i++) {
    ii = iMap[i];
    for (j=0; j<3; j++) {
      jj = iMap[j];
      twoDtangentCopy(ii,jj) = twoDtangent(i,j);
    }//end for j
  }//end for i

  for (i=0; i<2; i++) 
    for (j=0; j<2; j++) 
      dd11(i,j) = twoDtangentCopy(i,j);

  for (i = 0; i < 2; i++)
    dd12(i,0) = twoDtangentCopy(i,2);

  for (j = 0; j < 2; j++)
    dd21(0,j) = twoDtangentCopy(2,j);

  dd22(0,0)   = twoDtangentCopy(2,2);

  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  dd22.Solve(dd21, dd22invdd21);
  this->tangent   = dd11; 
  this->tangent  -= (dd12*dd22invdd21);

  return this->tangent;
}

const Matrix&  
PlaneStressFiberMaterial::getInitialTangent()
{
  static Matrix dd11(2,2);
  static Matrix dd12(2,1);
  static Matrix dd21(1,2);
  static Matrix dd22(1,1);
  static Matrix dd22invdd21(1,2);
  static Matrix twoDtangentCopy(3,3);

  const Matrix &twoDtangent = theMaterial->getInitialTangent();

  //swap matrix indices to sort out-of-plane components 
  int i, j, ii, jj;
  for (i=0; i<3; i++) {
    ii = iMap[i];
    for (j=0; j<3; j++) {
      jj = iMap[j];
      twoDtangentCopy(ii,jj) = twoDtangent(i,j);
    }//end for j
  }//end for i

  for (i=0; i<2; i++) 
    for (j=0; j<2; j++) 
      dd11(i,j) = twoDtangentCopy(i,j);

  for (i = 0; i < 2; i++)
    dd12(i,0) = twoDtangentCopy(i,2);

  for (j = 0; j < 2; j++)
    dd21(0,j) = twoDtangentCopy(2,j);

  dd22(0,0)   = twoDtangentCopy(2,2);

  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  dd22.Solve(dd21, dd22invdd21);
  this->tangent   = dd11; 
  this->tangent  -= (dd12*dd22invdd21);

  return this->tangent;
}

void  
PlaneStressFiberMaterial::Print(OPS_Stream &s, int flag)
{
	s << "PlaneStressFiberMaterial, tag: " << this->getTag() << endln;
	s << "\tWrapped material: "<< theMaterial->getTag() << endln;

	theMaterial->Print(s, flag);
}

int 
PlaneStressFiberMaterial::sendSelf(int commitTag, Channel &theChannel) 
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
    opserr << "PlaneStressFiberMaterial::sendSelf() - failed to send id data\n";
    return res;
  }

  // put the strains in a vector and send it
  static Vector vecData(1);
  vecData(0) = Cstrain22;

  res = theChannel.sendVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "PlaneStressFiberMaterial::sendSelf() - failed to send vector data\n";
    return res;
  }

  // now send the materials data
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "PlaneStressFiberMaterial::sendSelf() - failed to send vector material\n";

  return res;
}

int 
PlaneStressFiberMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  // recv an id contain the tag and associated materials class and db tags
  static ID idData(3);
  res = theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "PlaneStressFiberMaterial::sendSelf() - failed to send id data\n";
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
      opserr << "PlaneStressFiberMaterial::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMaterial->setDbTag(idData(2));

  // recv a vector containing strains and set the strains
  static Vector vecData(1);
  res = theChannel.recvVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "PlaneStressFiberMaterial::sendSelf() - failed to send vector data\n";
    return res;
  }

  Cstrain22 = vecData(0);

  Tstrain22 = Cstrain22;

  // now receive the materials data
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) 
    opserr << "PlaneStressFiberMaterial::sendSelf() - failed to send vector material\n";
  
  return res;
}
