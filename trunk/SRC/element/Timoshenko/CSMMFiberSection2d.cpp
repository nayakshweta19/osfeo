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

// $Source: /usr/local/cvs/OpenSees/SRC/element/Timoshenko/CSMMFiberSection2d.cpp,v $
// $Revision: 1.1 $
// $Date: 2009/01/10 21:22:20 $

// Created: 09/09
// Modified by: Li Ning 
// Description: This file contains the class implementation of CSMMFiberSection2d.Based on FiberSection2d.cpp.


#include <stdlib.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include "CSMMFiberSection2d.h"
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <UniaxialMaterial.h>
#include <NDMaterial.h>

#include <float.h>

#define MyMIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MyMAX(a, b)  (((a) > (b)) ? (a) : (b))
#define sign(a)    (((a) < (0)) ?(-1) : (1))

ID CSMMFiberSection2d::code(3);            

// constructors:
CSMMFiberSection2d::CSMMFiberSection2d(int tag, int num, Fiber **fibers, int strip1, double t1, int strip2, double t2, int strip3, double t3): 
  SectionForceDeformation(tag, SEC_TAG_CSMMFiberSection2d),
  numFibers(num), theMaterials(0), matData(0),
  NStrip(strip1+strip2+strip3), NStrip1(strip1), tavg1(t1), NStrip2(strip2), tavg2(t2), NStrip3(strip3), tavg3(t3), 
  StripCenterLoc(100), StripLoc(100,1000), FiberLoc(1000),
  yBar(0.0), sectionIntegr(0), ymax(0.0), ymin(0.0), e(3), eCommit(3), s(0), ks(0), sigmaY(0), tau(0), alpha(0), alphaCommit(0), 
  iterFile(0),exf(0),e1f(0),e2f(0),eyf(0),sxf(0),s1f(0),s2f(0),syf(0)
{
  
    if (numFibers != 0) {
    theMaterials = new NDMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "CSMMFiberSection2d::CSMMFiberSection2d -- failed to allocate Material pointers";
      exit(-1);
    }

    matData = new double [numFibers*2];

    if (matData == 0) {
      opserr << "CSMMFiberSection2d::CSMMFiberSection2d -- failed to allocate double array for material data\n";
      exit(-1);
    }

    double Qz = 0.0;
    double A  = 0.0;
    for (int i = 0; i < numFibers; i++) {
      Fiber *theFiber = fibers[i];
      double yLoc, zLoc, Area;
      theFiber->getFiberLocation(yLoc, zLoc);
      Area = theFiber->getArea();
      A  += Area;
      Qz += yLoc*Area;
      matData[i*2] = -yLoc;
      matData[i*2+1] = Area;
      NDMaterial *theMat = theFiber->getNDMaterial();        
      theMaterials[i] = theMat->getCopy();
 
      if (theMaterials[i] == 0) {
    opserr << "CSMMFiberSection2d::CSMMFiberSection2d -- failed to get copy of a Material\n";
    exit(-1);
      }

    if (-yLoc>ymax) ymax=-yLoc;
    if (-yLoc<ymin) ymin=-yLoc;
    }    
    yBar = -Qz/A;  
  }

  double YLoc[100];
  
  int ycount = 0;
  int loc = 0;
  for (int i = 0; i < numFibers; i++) {
    double y = matData[loc++];
    double A = matData[loc++];
    if (i==0) {
      YLoc[0] = y;
      ycount += 1; 
    } else {
      if (fabs(YLoc[ycount-1] - y) >= DBL_EPSILON) {
    YLoc[ycount] = y;
    ycount += 1;
      }
    }
    FiberLoc(i)=ycount-1;
  }
  
  if (ycount != NStrip) {
    opserr <<  "\n Failed - Not consistent number of fibers \n";
    exit(-1);
  }
    
  for (int j = 0; j < NStrip; j++)
	StripCenterLoc(j) = +(YLoc[j] - yBar);
  
  for (int k = 0; k < NStrip; k++) {
    int count=0;
    double Ac=0.0;
    for (int i = 0; i < numFibers; i++) {
      if (FiberLoc(i)==k){
    count++;    
    StripLoc(k,count+1)=i;
    Ac += matData[2*i+1];
      }
    }
    StripLoc(k,0)=count;    //num fibers in strip
    StripLoc(k,1)=Ac;        //total concrete area in strip
  }
  
  s = new Vector(sData, 3);
  ks = new Matrix(kData, 3, 3);

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;

  kData[0] = 0.0;
  kData[1] = 0.0;
  kData[2] = 0.0;
  kData[3] = 0.0;
  kData[4] = 0.0;
  kData[5] = 0.0;
  kData[6] = 0.0;
  kData[7] = 0.0;
  kData[8] = 0.0;

  code(0) = SECTION_RESPONSE_P;                            
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_VY;

// AddingSensitivity:BEGIN ////////////////////////////////////
  parameterID = 0;
// AddingSensitivity:END //////////////////////////////////////

}

// constructor for blank object that recvSelf needs to be invoked upon
CSMMFiberSection2d::CSMMFiberSection2d():
  SectionForceDeformation(0, SEC_TAG_CSMMFiberSection2d),
  numFibers(0), theMaterials(0), matData(0),
  NStrip1(0), tavg1(0.0), NStrip2(0), tavg2(0.0), NStrip3(0), tavg3(0.0), StripCenterLoc(100), StripLoc(100,1000), FiberLoc(1000),
  yBar(0.0), ymax(0.0), ymin(0.0), e(3), eCommit(3), s(0), ks(0), sigmaY(0), tau(0), alpha(0), alphaCommit(0), iterFile(0), exf(0), e1f(0),e2f(0),eyf(0), sxf(0), s1f(0),s2f(0),syf(0)
{
  s = new Vector(sData, 3);
  ks = new Matrix(kData, 3, 3);

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;

  kData[0] = 0.0;
  kData[1] = 0.0;
  kData[2] = 0.0;
  kData[3] = 0.0;
  kData[4] = 0.0;
  kData[5] = 0.0;
  kData[6] = 0.0;
  kData[7] = 0.0;
  kData[8] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_VY;

// AddingSensitivity:BEGIN ////////////////////////////////////
  parameterID = 0;
// AddingSensitivity:END //////////////////////////////////////
}


// destructor:
CSMMFiberSection2d::~CSMMFiberSection2d()
{
  if (theMaterials != 0) {
    for (int i = 0; i < numFibers; i++)
      if (theMaterials[i] != 0) {
    delete theMaterials[i];
      }
    delete [] theMaterials;
  }
  
  if (matData != 0)
    delete [] matData;
  
  if (s != 0)
    delete s;
  
  if (ks != 0)
    delete ks;
  
  if (sigmaY != 0)
    delete sigmaY;
  
  if (tau != 0)
    delete tau;
  
  if (alpha != 0)
    delete alpha;
  
  if (alphaCommit != 0)
    delete alphaCommit;
  
  if (iterFile != 0)
    delete iterFile;
  
  if (exf != 0)
    delete exf;
  
  if (e1f != 0)
    delete e1f;
  
  if (e2f != 0)
    delete e2f;
  
  if (eyf != 0)
    delete eyf;
  
  if (sxf != 0)
    delete sxf;
  
  if (s1f != 0)
    delete s1f;
  
  if (s2f != 0)
    delete s2f;
  
  if (syf != 0)
    delete syf;
}

// for dispBeamColumn and nonlinearBeamColumn
int
CSMMFiberSection2d::setTrialSectionDeformation (const Vector &deforms)    
{
  int res = 0;
  
  e = deforms;    // axial strain, curvature and shear trail strain
  
  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0;
  kData[3] = 0.0; kData[4] = 0.0; kData[5] = 0.0;
  kData[6] = 0.0; kData[7] = 0.0; kData[8] = 0.0;
  
  sData[0] = 0.0; sData[1] = 0.0; sData[2] = 0.0;
  
  double d0 = e(0);  // axial strain ; N
  double d1 = e(1);  // curvature ;    M
  double d2 = e(2);  // shear strain ; V
  
  double Fx=0.0, Fy=0.0, Fxy=0.0;
  double stiff11=0.0, stiff12=0.0, stiff22=0.0, stiff21=0.0;
  Vector strain(3);
  Matrix D; // nDMat Tangent stiffness matrix
  Vector d; // nDMat Stress Vector

  for (int jj = 0; jj < NStrip; jj++) {

    int iterNum = 0;
	int enlarge = 0;

    double tavg; //depth for different strip
    if (jj<NStrip1)
      tavg=tavg1;
    else {
      if (jj<NStrip2+NStrip1)
    tavg=tavg2;
      else
    tavg=tavg3;
    }

    double y     = StripCenterLoc(jj);
    double ex    = d0 + y*d1;  //axial strain
    double gamma = d2;         //shear strain
    
    if ((fabs(gamma) <= DBL_EPSILON) && (fabs(ex) <= DBL_EPSILON) ) {
      *ks = this->getInitialTangent();
  
      sData[0] = 0.0;
      sData[1] = 0.0;
      sData[2] = 0.0;
  
      continue;
    }

  	//for (int i = 0; i < numFibers; i++) {  

	int fibNum = jj;
    NDMaterial *theMat = theMaterials[fibNum];  
    //int tag=theMat->getTag();
    //double y = matData[2*fibNum] - yBar;
    double A = matData[2*fibNum+1];

    int status = 0; // status to check if iteration for ey

	// get the nDMat tangent matrix
    D = theMat->getInitialTangent();

	//set trial strain of nDMaterial, transverse strain is tentative value
    double ey = -(ex*D(1,0)+gamma*D(1,2))/D(1,1);
	strain(0) = ex;
	strain(1) = ey;   //tentative value for transverse strain
    strain(2) = gamma;
	
	opserr << "strip No.: " << jj+1 << endln;
	opserr << "strain: " << strain(0) << "\t" << strain(1) << "\t" << strain(2) << endln;
	
	// for transverse stress equal to zero (sigma_y == 0)
	// iterative procedure for section state determination
	
	status = 1;
	int iterMax = 200;
	double tol = 1e-4;
	int count = 0;

	double trialE0 = -0.2, trialE1 = ey, trialE2 = 0.1; // min, middle, max divided
	double trialS0 = 0.0, trialS1 = 0.0, trialS2 = 0.0;

	do {
      count++;
	  // iteration for transverse strain (1)
	  // minimum
	  strain(1) = trialE0;
	  
	  if ( theMat->setTrialStrain(strain) ) {
	  	opserr << "CSMMFiberSection2d::setTrialSectionDeformation - section failed in nDMat->setTrialStrain\n";
	  	return -1;
	  }
	  
	  trialS0 = (theMat->getStress())(1);
	  
	  // middle
	  strain(1) = trialE1;
	  
	  if ( theMat->setTrialStrain(strain) ) {
	  	opserr << "CSMMFiberSection2d::setTrialSectionDeformation - section failed in nDMat->setTrialStrain\n";
	  	return -1;
	  }
	  
	  trialS1 = (theMat->getStress())(1);
	  
	  //maximum
	  strain(1) = trialE2;
	  
	  if ( theMat->setTrialStrain(strain) ) {
	  	opserr << "CSMMFiberSection2d::setTrialSectionDeformation - section failed in nDMat->setTrialStrain\n";
	  	return -1;
	  }
	  
	  trialS2 = (theMat->getStress())(1);
	  
	  // divide the region into more  zero point
	  if ( sign(trialS0) != sign(trialS1) || sign(trialS1) != sign(trialS2)) {
		opserr << "find ONE zero point" << endln;
		
		if (sign(trialS0) != sign(trialS1)) {
	      trialE2 = trialE1;
		  trialE1 = (trialE1+trialE0)/2.0;
	    }
		if ( sign(trialS1) != sign(trialS2)) {
		  trialE0 = trialE1;
		  trialE1 = (trialE1+trialE2)/2.0;
		}
		if ( fabs(trialE2-trialE0) <= 10.*DBL_EPSILON ) {
		  //for nearly equivalent trial strain
			status = 0;
		}
	  } else if ( sign(trialS0) != sign(trialS1) && sign(trialS1) != sign(trialS2)) {
		opserr << "find TWO zero points" << endln;
	    status = 0;
	  } else {
		opserr << "connot find zero points" << endln;
	  }
	  
	  if ( (trialS0 < 0.0 && trialS1 < 0.0 && trialS2 < 0.0) ||
		   (trialS0 > 0.0 && trialS1 > 0.0 && trialS2 > 0.0) ) {

	    if ( trialS0 < 0.0 )      opserr << "all negative stress points" << endln;
		else if ( trialS0 > 0.0 ) opserr << "all positive stress points" << endln;

        trialE2 = trialE2*2.0;
		trialE0 = trialE0*2.0;

		enlarge ++;

		if ( enlarge >3) status = 0;

	  }

	  opserr << "trial transverse strain: " << trialE1 << endln;
	  opserr << "transverse stress:       " << trialS1 << endln;

	  // determine the final sectional strain state
	  if ( fabs(trialS1) <= tol ) {
		strain(1) =trialE1;
		status = 0;
	  }

	} while ((status != 0 ) && (count <= iterMax));
	
	// Assemble the stiffness matrix for section level
    if ( theMat->setTrialStrain(strain) ) {
	  opserr << "CSMMFiberSection2d::setTrialSectionDeformation - section failed in nDMat->setTrialStrain\n";
	  return -1;
	}

	D = theMat->getTangent();

	d = theMat->getStress();

    stiff11 += D(0,0)*A;
    stiff12 += D(0,2)*A;
    stiff21 += D(2,0)*A;
    stiff22 += D(2,2)*A;

    Fx += d(0)*A;
    Fy += d(1)*A;
    Fxy+= d(2)*A;

    //}

    kData[0] += stiff11;
    kData[1] += stiff11*y;
    kData[2] += stiff12;
    kData[3] += stiff11*y;
    kData[4] += stiff11*y*y;
    kData[5] += stiff12*y;
    kData[6] += stiff21;
    kData[7] += stiff21*y;
    kData[8] += stiff22;
    
    sData[0] += Fx;
    sData[1] += Fx*y;
    sData[2] += Fxy;
  }

  return res;
}

// for Timoshenko and dispBeamColumnInt
int 
CSMMFiberSection2d::setTrialSectionDeformation (const Vector &deforms, double L)    
{
  opserr << "error in CSMMFiberSection2d::setTrialSectionDeformation (const Vector &deforms, double L) ";
  return 0;

}

int
CSMMFiberSection2d::revertToLastCommit (void)    
{
  int res = 0;
 
  e = eCommit;    // axial strain, curvature and shear strain

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0;
  kData[3] = 0.0; kData[4] = 0.0; kData[5] = 0.0;
  kData[6] = 0.0; kData[7] = 0.0; kData[8] = 0.0;
  
  sData[0] = 0.0; sData[1] = 0.0; sData[2] = 0.0;
  
  double d0 = e(0);  // axial strain ; N
  double d1 = e(1);  // curvature ;    M
  double d2 = e(2);  // shear strain ; V
  
  int iterMax = 200;
  double tol = 1e-4;
  
  double Fx=0.0, Fy=0.0, Fxy=0.0;
  double stiff11=0.0, stiff12=0.0, stiff22=0.0, stiff21=0.0;
  Vector strain(3);
  Matrix D; // nDMat Tangent stiffness matrix
  Vector d; // nDMat Stress Vector

  for (int jj = 0; jj < NStrip; jj++) {
    int iterNum = 0;
	int enlarge = 0;
    double tavg; //depth for different strip
    if (jj<NStrip1)
      tavg=tavg1;
    else {
      if (jj<NStrip2+NStrip1)
    tavg=tavg2;
      else
    tavg=tavg3;
    }
    
	double y     = StripCenterLoc(jj);
    double ex    = d0 + y*d1;  //axial strain
    double gamma = d2;         //shear strain
    
    if ((fabs(gamma) <= DBL_EPSILON) && (fabs(ex) <= DBL_EPSILON) ) {
      *ks = this->getInitialTangent();
  
      sData[0] = 0.0;
      sData[1] = 0.0;
      sData[2] = 0.0;
  
      continue;
    }

  	//for (int i = 0; i < numFibers; i++) {  

	int fibNum = jj;
    NDMaterial *theMat = theMaterials[fibNum];  
    //int tag=theMat->getTag();
    //double y = matData[2*fibNum] - yBar;
    double A = matData[2*fibNum+1];

    int status = 0; // status to check if iteration for ey

	// get the nDMat tangent matrix
    D = theMat->getInitialTangent();

	//set trial strain of nDMaterial, transverse strain is tentative value
    double ey = -(ex*D(1,0)+gamma*D(1,2))/D(1,1);
	strain(0) = ex;
	strain(1) = ey;   //tentative value for transverse strain
    strain(2) = gamma;
	
	opserr << "strip No.: " << jj+1 << endln;
	opserr << "strain: " << strain(0) << "\t" << strain(1) << "\t" << strain(2) << endln;
	
	// for transverse stress equal to zero (sigma_y == 0)
	// iterative procedure for section state determination
	
	status = 1;
	double trialE0 = -0.2, trialE1 = ey, trialE2 = 0.1; // min, middle, max divided
	double trialS0 = 0.0, trialS1 = 0.0, trialS2 = 0.0;
	while(status) {
	  // iteration for transverse strain (1)
	  // minimum
	  strain(1) = trialE0;
	  
	  if ( theMat->setTrialStrain(strain) ) {
	  	opserr << "CSMMFiberSection2d::setTrialSectionDeformation - section failed in nDMat->setTrialStrain\n";
	  	return -1;
	  }
	  
	  trialS0 = (theMat->getStress())(1);
	  
	  // middle
	  strain(1) = trialE1;
	  
	  if ( theMat->setTrialStrain(strain) ) {
	  	opserr << "CSMMFiberSection2d::setTrialSectionDeformation - section failed in nDMat->setTrialStrain\n";
	  	return -1;
	  }
	  
	  trialS1 = (theMat->getStress())(1);
	  
	  //maximum
	  strain(1) = trialE2;
	  
	  if ( theMat->setTrialStrain(strain) ) {
	  	opserr << "CSMMFiberSection2d::setTrialSectionDeformation - section failed in nDMat->setTrialStrain\n";
	  	return -1;
	  }
	  
	  trialS2 = (theMat->getStress())(1);
	  
	  // divide the region into more  zero point
	  if ( sign(trialS0) != sign(trialS1) || sign(trialS1) != sign(trialS2)) {
		opserr << "find ONE zero point" << endln;
		if (sign(trialS0) != sign(trialS1)) {
	      trialE2 = trialE1;
		  trialE1 = (trialE1+trialE0)/2.0;
	    }
		if ( sign(trialS1) != sign(trialS2)) {
		  trialE0 = trialE1;
		  trialE1 = (trialE1+trialE2)/2.0;
		}
	  } else if ( sign(trialS0) != sign(trialS1) && sign(trialS1) != sign(trialS2)) {
		opserr << "find TWO zero points" << endln;
	    status = 0;
	  } else {

	  }
	  
	  if ( (trialS0 < 0.0 && trialS1 < 0.0 && trialS2 < 0.0) ||
		   (trialS0 > 0.0 && trialS1 > 0.0 && trialS2 > 0.0) ) {
	    if ( trialS0 < 0.0 ) opserr << "all negative stress points" << endln;
		else                 opserr << "all positive stress points" << endln;
        trialE2 = trialE2*2.0;
		trialE0 = trialE0*2.0;
		enlarge ++;
		if ( enlarge >3) status = 0;
	  }
	  opserr << "trial transverse strain: " << trialE1 << endln;
	  opserr << "transverse stress:       " << trialS1 << endln;
	  // determine the final sectional strain state
	  if ( fabs(trialS1) <= tol ) {
		strain(1) =trialE1;
		status = 0;
	  }

	  //iterNum ++; 
	  //if ( iterNum > iterMax ) {
		//
	  //}
	
	  

	}
	
	// Assemble the stiffness matrix for section level
    if ( theMat->setTrialStrain(strain) ) {
	  opserr << "CSMMFiberSection2d::setTrialSectionDeformation - section failed in nDMat->setTrialStrain\n";
	  return -1;
	}

	D = theMat->getTangent();

	d = theMat->getStress();

    stiff11 += D(0,0)*A;
    stiff12 += D(0,2)*A;
    stiff21 += D(2,0)*A;
    stiff22 += D(2,2)*A;

    Fx += d(0)*A;
    Fy += d(1)*A;
    Fxy+= d(2)*A;

    //}
	
    kData[0] += stiff11;
    kData[1] += stiff11*y;
    kData[2] += stiff12;
    kData[3] += stiff11*y;
    kData[4] += stiff11*y*y;
    kData[5] += stiff12*y;
    kData[6] += stiff21;
    kData[7] += stiff21*y;
    kData[8] += stiff22;
    
    sData[0] += Fx;
    sData[1] += Fx*y;
    sData[2] += Fxy;
  }

  return res;
}

int
CSMMFiberSection2d::commitState(void)            
{
  int err = 0;
  
  for (int i = 0; i < numFibers; i++){
    err += theMaterials[i]->commitState();
  }

  eCommit = e;

  return err;
}


const Vector&
CSMMFiberSection2d::getSigmaY(void)    
{
  return *sigmaY;
}

const Vector&
CSMMFiberSection2d::getTau(void)    
{
  return *tau;
}

const Vector&
CSMMFiberSection2d::getAlpha(void)    
{
  return *alpha;
}

const Vector&
CSMMFiberSection2d::getIter(void)    
{
  return *iterFile;
}

const Vector&
CSMMFiberSection2d::getEX(void)    
{
  return *exf;
}

const Vector&
CSMMFiberSection2d::getEY(void)    
{
  return *eyf;
}

const Vector&
CSMMFiberSection2d::getE1(void)    
{
  return *e1f;
}

const Vector&
CSMMFiberSection2d::getE2(void)    
{
  return *e2f;
}

const Vector&
CSMMFiberSection2d::getSX(void)    
{
  return *sxf;
}

const Vector&
CSMMFiberSection2d::getSY(void)    
{
  return *syf;
}

const Vector&
CSMMFiberSection2d::getS1(void)
{
  return *s1f;
}

const Vector&
CSMMFiberSection2d::getS2(void)    
{
  return *s2f;
  
}


const Vector &
CSMMFiberSection2d::getSectionDeformation (void)
{
    return e;
}

const Matrix&
CSMMFiberSection2d::getInitialTangent(void)    
{
  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0;
  kData[3] = 0.0; kData[4] = 0.0; kData[5] = 0.0;
  kData[6] = 0.0; kData[7] = 0.0; kData[8] = 0.0;

  for (int i = 0; i < numFibers; i++) {
    
    NDMaterial *theMat = theMaterials[i];  
    
    //int tag=theMat->getTag();
    double y = StripCenterLoc(FiberLoc(i));
    double A = matData[2*i+1];
    
    Matrix tangent = theMat->getInitialTangent();

    double d11 = tangent(0,0);
    //double d12 = tangent(0,1);
    double d13 = tangent(0,2);
    //double d21 = tangent(1,0);
    //double d22 = tangent(1,1);
    //double d23 = tangent(1,2);
    double d31 = tangent(2,0);
    //double d32 = tangent(2,1);
    double d33 = tangent(2,2);

    kData[0] += d11*A;
    kData[1] += d11*A*y;
    kData[2] += d13*A;
    kData[3] += d11*A*y;
    kData[4] += d11*A*y*y;
    kData[5] += d13*A*y;
    kData[6] += d31*A;
    kData[7] += d31*A*y;
    kData[8] += d33*A;
  }
  
  return *ks;
}

const Matrix&
CSMMFiberSection2d::getSectionTangent(void)    
{
  return *ks;
}

const Vector&
CSMMFiberSection2d::getStressResultant(void)    
{
  return *s;
}

SectionForceDeformation*
CSMMFiberSection2d::getCopy(void)
{
  CSMMFiberSection2d *theCopy = new CSMMFiberSection2d ();
  theCopy->setTag(this->getTag());

  theCopy->numFibers = numFibers;

  if (numFibers != 0) {
    theCopy->theMaterials = new NDMaterial *[numFibers];

    if (theCopy->theMaterials == 0) {
      opserr <<"CSMMFiberSection2d::getCopy -- failed to allocate Material pointers\n";
      exit(-1);
    }
  
    theCopy->matData = new double [numFibers*2];

    if (theCopy->matData == 0) {
      opserr << "CSMMFiberSection2d::getCopy -- failed to allocate double array for material data\n";
      exit(-1);
    }
                
    for (int i = 0; i < numFibers; i++) {
      theCopy->matData[i*2] = matData[i*2];
      theCopy->matData[i*2+1] = matData[i*2+1];
      theCopy->theMaterials[i] = theMaterials[i]->getCopy();
        
      if (theCopy->theMaterials[i] == 0) {
    opserr <<"CSMMFiberSection2d::getCopy -- failed to get copy of a Material";
    exit(-1);
      }
    }  
  }

  theCopy->NStrip = NStrip;
  theCopy->NStrip1 = NStrip1;
  theCopy->NStrip2 = NStrip2;
  theCopy->NStrip3 = NStrip3;
  theCopy->tavg1 = tavg1;
  theCopy->tavg2 = tavg2;
  theCopy->tavg3 = tavg3;
  
  for (int j = 0; j < NStrip; j++) {
    theCopy->sy[j] = sy[j];
    theCopy->txy[j] = txy[j];
    theCopy->alfa[j] = alfa[j]; 
    theCopy->alfaCommit[j] = alfaCommit[j];
    theCopy->iterOut[j] = iterOut[j];
    theCopy->iterCommit[j] = iterCommit[j];
    theCopy->exOut[j] = exOut[j];
    theCopy->exCommit[j] = exCommit[j];
    theCopy->eyCommit[j] = eyCommit[j];
    theCopy->e1Commit[j] = e1Commit[j];
    theCopy->e2Commit[j] = e2Commit[j];
    theCopy->sxCommit[j] = sxCommit[j];
    theCopy->syCommit[j] = syCommit[j];
    theCopy->s1Commit[j] = s1Commit[j];
    theCopy->s2Commit[j] = s2Commit[j];
  }

  theCopy->StripCenterLoc = StripCenterLoc;
  theCopy->StripLoc = StripLoc;
  theCopy->FiberLoc = FiberLoc;

  theCopy->eCommit = eCommit;
  theCopy->e = e;
  theCopy->yBar = yBar;
  theCopy->ymin = ymin;
  theCopy->ymax = ymax;

  theCopy->kData[0] = kData[0];
  theCopy->kData[1] = kData[1];
  theCopy->kData[2] = kData[2];
  theCopy->kData[3] = kData[3];
  theCopy->kData[4] = kData[4];
  theCopy->kData[5] = kData[5];
  theCopy->kData[6] = kData[6];
  theCopy->kData[7] = kData[7];
  theCopy->kData[8] = kData[8];

  theCopy->sData[0] = sData[0];
  theCopy->sData[1] = sData[1];
  theCopy->sData[2] = sData[2];

  return theCopy;
}

const ID&
CSMMFiberSection2d::getType ()            
{
  return code;
}

int
CSMMFiberSection2d::getOrder () const    
{
  return 3;                                
}


int
CSMMFiberSection2d::revertToStart(void)        
{
  // revert the fibers to start    
  int err = 0;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0;
  kData[3] = 0.0; kData[4] = 0.0; kData[5] = 0.0;
  kData[6] = 0.0; kData[7] = 0.0; kData[8] = 0.0;

  sData[0] = 0.0; sData[1] = 0.0;    sData[2] = 0.0;
  
    
  for (int i = 0; i < numFibers; i++) {
    
    NDMaterial *theMat = theMaterials[i];  
    
    //int tag=theMat->getTag();
    double y = StripCenterLoc(FiberLoc(i));
    double A = matData[2*i+1];
    
    Matrix tangent = theMat->getInitialTangent();

    double d11 = tangent(0,0);
    //double d12 = tangent(0,1);
    double d13 = tangent(0,2);
    //double d21 = tangent(1,0);
    //double d22 = tangent(1,1);
    //double d23 = tangent(1,2);
    double d31 = tangent(2,0);
    //double d32 = tangent(2,1);
    double d33 = tangent(2,2);

    kData[0] += d11*A;
    kData[1] += d11*A*y;
    kData[2] += d13*A;
    kData[3] += d11*A*y;
    kData[4] += d11*A*y*y;
    kData[5] += d13*A*y;
    kData[6] += d31*A;
    kData[7] += d31*A*y;
    kData[8] += d33*A;

    sData[0] += 0.0; 
    sData[1] += 0.0; 
    sData[2] += 0.0; 
  }

  return err;
}

int
CSMMFiberSection2d::sendSelf(int commitTag, Channel &theChannel)    
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
    opserr <<  "CSMMFiberSection2d::sendSelf - failed to send ID data\n";
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
      opserr <<  "CSMMFiberSection2d::sendSelf - failed to send material data\n";
      return res;
    }    

    // send the fiber data, i.e. area and loc
    Vector fiberData(matData, 2*numFibers);
    res += theChannel.sendVector(dbTag, commitTag, fiberData);
    if (res < 0) {
      opserr <<  "CSMMFiberSection2d::sendSelf - failed to send material data\n";
      return res;
    }    

    // now invoke send(0 on all the materials
    for (int j=0; j<numFibers; j++){
      theMaterials[j]->sendSelf(commitTag, theChannel);
    }
  }
  return res;
}

int
CSMMFiberSection2d::recvSelf(int commitTag, Channel &theChannel,
             FEM_ObjectBroker &theBroker)                        
{
  int res = 0;
  static ID data(3);
  
  int dbTag = this->getDbTag();
  res += theChannel.recvID(dbTag, commitTag, data);
  if (res < 0) {
    opserr <<  "CSMMFiberSection2d::recvSelf - failed to recv ID data\n";
    return res;
  }    
  this->setTag(data(0));

  // recv data about materials objects, classTag and dbTag
  if (data(1) != 0) {
    ID materialData(2*data(1));
    res += theChannel.recvID(dbTag, commitTag, materialData);
    if (res < 0) {
      opserr <<  "CSMMFiberSection2d::recvSelf - failed to recv material data\n";
      return res;
    }    

    // if current arrays not of correct size, release old and resize
    if (theMaterials == 0 || numFibers != data(1)) {
      // delete old stuff if outa date
      if (theMaterials != 0) {
    for (int i=0; i<numFibers; i++){
      delete theMaterials[i];
    }
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
      opserr <<"CSMMFiberSection2d::recvSelf -- failed to allocate Material pointers\n";
      exit(-1);
    }
    
    for (int j=0; j<numFibers; j++){
      theMaterials[j] = 0;
    }
    matData = new double [numFibers*2];

    if (matData == 0) {
      opserr <<"CSMMFiberSection2d::recvSelf  -- failed to allocate double array for material data\n";
      exit(-1);
    }
      }
    }

    Vector fiberData(matData, 2*numFibers);
    res += theChannel.recvVector(dbTag, commitTag, fiberData);
    if (res < 0) {
      opserr <<  "CSMMFiberSection2d::recvSelf - failed to recv material data\n";
      return res;
    }    

    int i;
    for (i=0; i<numFibers; i++) {
      int classTag = materialData(2*i);
      int dbTag = materialData(2*i+1);

      // if material pointed to is blank or not of corrcet type, 
      // release old and create a new one
      if (theMaterials[i] == 0){
    theMaterials[i] = theBroker.getNewNDMaterial(classTag);

      }
      else if (theMaterials[i]->getClassTag() != classTag) {
    delete theMaterials[i];
    theMaterials[i] = theBroker.getNewNDMaterial(classTag);       
    
      }

      if (theMaterials[i] == 0) {
    opserr <<"CSMMFiberSection2d::recvSelf -- failed to allocate double array for material data\n";
    exit(-1);
      }

      theMaterials[i]->setDbTag(dbTag);
      res += theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
    }

    double Qz = 0.0;
    double A  = 0.0;
    double yLoc, Area;

    // Recompute centroid
    for (i = 0; i < numFibers; i++) {
      yLoc = -matData[2*i];
      Area = matData[2*i+1];
      A  += Area;
      Qz += yLoc*Area;
    }
    
    yBar = -Qz/A;
  }    

  return res;
}

void
CSMMFiberSection2d::Print(OPS_Stream &s, int flag)
{
  s << "\nFiberSection2d, tag: " << this->getTag() << endln;
  s << "\tSection code: " << code;
  s << "\tNumber of Fibers: " << numFibers << endln;
  s << "\tCentroid: " << -yBar << endln;

  if (flag == 1) {
    int loc = 0;
    for (int i = 0; i < numFibers; i++) {
      s << "\nLocation (y) = (" << -matData[loc++] << ")";
      s << "\nArea = " << matData[loc++] << endln;
      theMaterials[i]->Print(s, flag);
    }
  }
}

Response*
CSMMFiberSection2d::setResponse(const char **argv, int argc, OPS_Stream &s)    
{
  // See if the response is one of the defaults
  Response *res = SectionForceDeformation::setResponse(argv, argc, s);
  if (res != 0)
    return res;
  
  // Check if fiber response is requested
  else if ((strcmp(argv[0],"fiber") == 0) || (strcmp(argv[0],"fiber1") == 0)) {    
    int key = numFibers;
    int passarg = 2;
    
    if (argc <= 2)          
      return 0;

    if (argc <= 3) {          
      key = atoi(argv[1]);
      if (key < numFibers)
         return theMaterials[key]->setResponse(&argv[passarg],argc-passarg,s);
      else
         return 0;
    }

    if (argc > 4) {         // find fiber closest to coord. with mat tag
      int matTag = atoi(argv[3]);
      double yCoord = atof(argv[1]);
      double closestDist;
      double ySearch, dy;
      double distance;
      int j;
      // Find first fiber with specified material tag
      for (j = 0; j < numFibers; j++) {
    if (matTag == theMaterials[j]->getTag()) {
      ySearch = -matData[2*j];
      dy = ySearch-yCoord;
      closestDist = fabs(dy);
      key = j;
      break;
    }
      }
      // Search the remaining fibers
      for ( ; j < numFibers; j++) {
    if (matTag == theMaterials[j]->getTag()) {
      ySearch = -matData[2*j];
      dy = ySearch-yCoord;
      distance = fabs(dy);
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
      double closestDist;
      double ySearch, dy;
      double distance;
      ySearch = -matData[0];
      dy = ySearch-yCoord;
      closestDist = fabs(dy);
      key = 0;
      for (int j = 1; j < numFibers; j++) {
    ySearch = -matData[2*j];
    dy = ySearch-yCoord;
    distance = fabs(dy);
    if (distance < closestDist) {
      closestDist = distance;
      key = j;
    }
      }
      passarg = 3;
    }
    
    if (key < numFibers)
      return theMaterials[key]->setResponse(&argv[passarg],argc-passarg,s);
    else
      return 0;
  }


  else if (strcmp(argv[0],"fiber2") == 0) {
    int key = numFibers;
    int passarg = 2;
    
    if (argc <= 2)          
      return 0;

    if (argc <= 3) {          
      key = atoi(argv[1]);
      if (key < numFibers)
         return theMaterials[key]->setResponse(&argv[passarg],argc-passarg,s);
      else
         return 0;
    }

    if (argc > 4) {         // find fiber closest to coord. with mat tag
      int matTag = atoi(argv[3]);
      double yCoord = atof(argv[1]);
      double closestDist;
      double ySearch, dy;
      double distance;
      int j;
      // Find first fiber with specified material tag
      for (j = 0; j < numFibers; j++) {
    if (matTag == theMaterials[j]->getTag()) {
      ySearch = -matData[2*j];
      dy = ySearch-yCoord;
      closestDist = fabs(dy);
      key = j;
      break;
    }
      }
      // Search the remaining fibers
      for ( ; j < numFibers; j++) {
    if (matTag == theMaterials[j]->getTag()) {
      ySearch = -matData[2*j];
      dy = ySearch-yCoord;
      distance = fabs(dy);
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
      double closestDist;
      double ySearch, dy;
      double distance;
      ySearch = -matData[0];
      dy = ySearch-yCoord;
      closestDist = fabs(dy);
      key = 0;
      for (int j = 1; j < numFibers; j++) {
    ySearch = -matData[2*j];
    dy = ySearch-yCoord;
    distance = fabs(dy);
    if (distance < closestDist) {
      closestDist = distance;
      key = j;
    }
      }
      passarg = 3;
    }
    
    if (key < numFibers)
      return theMaterials[key]->setResponse(&argv[passarg],argc-passarg,s);
    else
      return 0;
  }
 
  else
    return 0;
}


int 
CSMMFiberSection2d::getResponse(int responseID, Information &sectInfo)
{
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  return SectionForceDeformation::getResponse(responseID, sectInfo);
}



// AddingSensitivity:BEGIN ////////////////////////////////////
const Vector &
CSMMFiberSection2d::getSectionDeformationSensitivity(int gradNumber)
{
    static Vector dummy(2);
    return dummy;
}

const Vector &
CSMMFiberSection2d::getStressResultantSensitivity(int gradNumber, bool conditional)
{
    static Vector dummy(2);    
    return dummy;
}

const Matrix &
CSMMFiberSection2d::getSectionTangentSensitivity(int gradNumber)
{
    static Matrix something(2,2);
    something.Zero();
    return something;
}

int
CSMMFiberSection2d::commitSensitivity(const Vector& defSens, int gradNumber, int numGrads)
{
    return 0;
}
// AddingSensitivity:END ///////////////////////////////////
