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

// $Revision: 1.6 $
// $Date: 2006/08/15 00:41:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ConcreteDPM01.cpp,v $

// Written: neallee
// Created: Aug 2014
//
// Description: This file contains the class implementation for 
// ConcreteDPM01 Material.

#include <elementAPI.h>

#include <ConcreteDPM01.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

void *
OPS_NewConcreteDPM01(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData1[1];
  double dData[13];
  int    iData2[1];
  int numData;

  numData = 1;
  if (OPS_GetIntInput(&numData, iData1) != 0) {
    opserr << "WARNING invalid uniaxialMaterial BWBN tag" << endln;
    return 0;
  }

  numData = 13;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid Double Values\n";
    return 0;	
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, iData2) != 0) {
    opserr << "WARNING invalid maxNumIter" << endln;
    return 0;
  }

  //theMaterial = new ConcreteDPM01(iData1[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
	//		 dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12],iData2[0]);       
  theMaterial = new ConcreteDPM01();

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type BWBN\n";
    return 0;
  }

  return theMaterial;
}

ConcreteDPM01::ConcreteDPM01(int tag)
  :UniaxialMaterial(tag,MAT_TAG_ConcreteDPM01),
   trialStrain(0.0), trialStress(0.0), trialTangent(0.0),
   commitStrain(0.0), commitStress(0.0), commitTangent(0.0)
{

}

ConcreteDPM01::ConcreteDPM01()
  :UniaxialMaterial(0,MAT_TAG_ConcreteDPM01),
   trialStrain(0.0), trialStress(0.0), trialTangent(0.0),
   commitStrain(0.0), commitStress(0.0), commitTangent(0.0)
{

}

ConcreteDPM01::~ConcreteDPM01()
{

}

int 
ConcreteDPM01::setTrialStrain(double strain, double strainRate)
{
  // set the trial strain
  trialStrain = strain;

  // determine trial stress and tangent
  trialStress = 0.0;
  trialTangent = 0.0;

  return 0;
}

double 
ConcreteDPM01::getStress(void)
{
  return trialStress;
}

double 
ConcreteDPM01::getTangent(void)
{
  return trialTangent;
}

double 
ConcreteDPM01::getInitialTangent(void)
{
  // return the initial tangent
  return 0.0;
}

double 
ConcreteDPM01::getStrain(void)
{
  return trialStrain;
}

int 
ConcreteDPM01::commitState(void)
{
  commitStrain  = trialStrain;
  commitStress  = trialStress;
  commitTangent = trialTangent;

  return 0;
}

int 
ConcreteDPM01::revertToLastCommit(void)
{
  trialStrain = commitStrain;
  trialStress = commitStress;
  trialTangent = commitTangent;

  return 0;
}

int 
ConcreteDPM01::revertToStart(void)
{
  trialStrain = 0.;
  trialStress = 0.0;
  trialTangent = 0.0;
  commitStrain = 0.;
  commitStress = 0.0;
  commitTangent = 0.0;

  return 0;
}

UniaxialMaterial *
ConcreteDPM01::getCopy(void)
{
  ConcreteDPM01 *theCopy = new ConcreteDPM01(this->getTag());

  return theCopy;
}

int 
ConcreteDPM01::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int 
ConcreteDPM01::recvSelf(int cTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  return -1;
}

void 
ConcreteDPM01::Print(OPS_Stream &s, int flag)
{
  s << "ConcreteDPM01 : " << this->getTag();

  return;
}


