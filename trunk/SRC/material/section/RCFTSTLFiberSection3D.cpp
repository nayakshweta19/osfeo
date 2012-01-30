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
                                                                        
// $Revision: 1.3 $
// $Date: 2008-11-06 20:52:20 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/OpenSEESComp/OpenSees/SRC/material/section/RCFTSTLFiberSection3D.cpp,v $
                                                                        
// Written: fmk
// Created: 04/04
//
// Description: This file contains the class implementation of FiberSection2d.

#include <stdlib.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <RCFTSTLFiberSection3D.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <UniaxialMaterial.h>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
using namespace::std;

//using std::ofstream;
//using std::ios;
//using std::endl;

ID RCFTSTLFiberSection3D::code(4);
Vector RCFTSTLFiberSection3D::s(4);
Matrix RCFTSTLFiberSection3D::ks(4,4);

// constructors:
RCFTSTLFiberSection3D::RCFTSTLFiberSection3D(int tag, int num, Fiber **fibers, double gj): 
  SectionForceDeformation(tag, SEC_TAG_RCFTSTLFiberSection3D),
  numFibers(num), theMaterials(0), matData(0),
  yBar(0.0), zBar(0.0), e(4), eCommit(4), GJ(gj)
{
  ofstream stlfib;
  stlfib.open("stlfib.dat",ios::app);
  double EA = 0.0;
  double EQz = 0.0;
  double EQy = 0.0;
  double EIz = 0.0;
  double EIy = 0.0;
  double EIyz = 0.0;
  
  double Es = 0.0;    
  if (numFibers != 0) {
    theMaterials = new UniaxialMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "RCFTSTLFiberSection3D::RCFTSTLFiberSection3D -- failed to allocate Material pointers\n";
      exit(-1);
    }

    matData = new double [numFibers*3];

    if (matData == 0) {
      opserr << "RCFTSTLFiberSection3D::RCFTSTLFiberSection3D -- failed to allocate double array for material data\n";
      exit(-1);
    }
    double Qz = 0.0;
    double Qy = 0.0;
    double A  = 0.0;
   
    for (int i = 0; i < numFibers; i++) {
      Fiber *theFiber = fibers[i];
      double yLoc, zLoc, Area;
      theFiber->getFiberLocation(yLoc, zLoc);
      Area = theFiber->getArea();
      UniaxialMaterial *theMat = theFiber->getMaterial();
      
      stlfib<<i<<"  "<<yLoc<<"  "<<zLoc<<"  "<<Area<<endl;
      Es = theMat->getInitialTangent();
      Qz += yLoc*Area;
      Qy += zLoc*Area;
      A  += Area;
      EQz += yLoc*Area*Es;
      EQy += zLoc*Area*Es;
      EA  += Area*Es;
      EIz += yLoc*yLoc*Area*Es;
      EIy += zLoc*zLoc*Area*Es;
      EIyz += yLoc*zLoc*Area*Es;
      
      matData[i*3] = -yLoc;
      matData[i*3+1] = zLoc;
      matData[i*3+2] = Area;
      theMaterials[i] = theMat->getCopy();

      if (theMaterials[i] == 0) {
	opserr << "RCFTSTLFiberSection3D::RCFTFiberSection3D -- failed to get copy of a Material\n";
	exit(-1);
      }
    }

    yBar = -Qz/A;
    zBar = Qy/A;
  }

  ks.Zero();

  ks(0,0) = EA;
  ks(0,1) = ks(1,0) = -EQz;
  ks(0,2) = ks(2,0) = -EQy;
  ks(1,1) = EIz;
  ks(2,2) = EIy;
  ks(2,1) = ks(1,2) = EIyz;
  ks(3,3) = GJ;

  //ks(0,0) = 449602.6;
  //ks(0,1) = ks(1,0) = 0.0;
  //ks(0,2) = ks(2,0) = 0.0;
  //ks(1,1) = 6968845.86;
  //ks(2,2) = 6968845.86;
  //ks(2,1) = ks(1,2) = 0.0;
  //ks(3,3) = GJ;

  //ks(0,0) = 1000.0;
  //ks(0,1) = ks(1,0) = 0.0;
  //ks(0,2) = ks(2,0) = 0.0;
  //ks(1,1) = 1.0;
  //ks(2,2) = 1.0;
  //ks(2,1) = ks(1,2) = 0.0;
  //ks(3,3) = GJ;

  //ks(0,0) = 300000.0;
  //ks(0,1) = ks(1,0) = 0.0;
  //ks(0,2) = ks(2,0) = 0.0;
  //ks(1,1) = 3000000.0;
  //ks(2,2) = 3000000.0;
  //ks(2,1) = ks(1,2) = 0.0;
  //ks(3,3) = GJ;

  //ks(0,0) = 43200000.0;
  //ks(0,1) = ks(1,0) = 0.0;
  //ks(0,2) = ks(2,0) = 0.0;
  //ks(1,1) = 14400000.0;
  //ks(2,2) = 14400000.0;
  //ks(2,1) = ks(1,2) = 0.0;
  //ks(3,3) = GJ;
  
  //ks(0,0) = 449509.3;
  //ks(0,1) = ks(1,0) = 0.0;
  //ks(0,2) = ks(2,0) = 0.0;
  //ks(1,1) = 6967684.16;
  //ks(2,2) = 6967684.16;
  //ks(2,1) = ks(1,2) = 0.0;
  //ks(3,3) = GJ;

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;

  kData[0] = EA;
  kData[1] = -EQz;
  kData[2] = EIyz;
  kData[3] = 0.0;
  kData[4] = -EQz;
  kData[5] = EIz;
  kData[6] = -EQy;
  kData[7] = 0.0;
  kData[8] = EIyz;
  kData[9] = -EQy;
  kData[10] = EIy;
  kData[11] = 0.0;
  kData[12] = 0.0;
  kData[13] = 0.0;
  kData[14] = 0.0;
  kData[15] = GJ;

  //kData[0] = 43200000.0;
  //kData[1] = 0.0;
  //kData[2] = 0.0;
  //kData[3] = 0.0;
  //kData[4] = 0.0;
  //kData[5] = 14400000.0;
  //kData[6] = 0.0;
  //kData[7] = 0.0;
  //kData[8] = 0.0;
  //kData[9] = 0.0;
  //kData[10] = 14400000.0;
  //kData[11] = 0.0;
  //kData[12] = 0.0;
  //kData[13] = 0.0;
  //kData[14] = 0.0;
  //kData[15] = GJ;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_T;
}

// constructor for blank object that recvSelf needs to be invoked upon
RCFTSTLFiberSection3D::RCFTSTLFiberSection3D():
  SectionForceDeformation(0, SEC_TAG_RCFTSTLFiberSection3D),
  numFibers(0), theMaterials(0), matData(0),
  yBar(0.0), zBar(0.0), e(4), eCommit(4), GJ(1.0)
{
  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;

  for (int i=0; i<6; i++)
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_T;
}

int
RCFTSTLFiberSection3D::addFiber(Fiber &newFiber)
{
  // need to create a larger array
  int newSize = numFibers+1;

  UniaxialMaterial **newArray = new UniaxialMaterial *[newSize]; 
  double *newMatData = new double [3 * newSize];
  
  if (newArray == 0 || newMatData == 0) {
    opserr << "RCFTSTLFiberSection3D::addFiber -- failed to allocate Fiber pointers\n";
    return -1;
  }

  // copy the old pointers
  int i;
  for (i = 0; i < numFibers; i++) {
    newArray[i] = theMaterials[i];
    newMatData[3*i] = matData[3*i];
    newMatData[3*i+1] = matData[3*i+1];
    newMatData[3*i+2] = matData[3*i+2];
  }
  // set the new pointers
  double yLoc, zLoc, Area;
  newFiber.getFiberLocation(yLoc, zLoc);
  Area = newFiber.getArea();
  newMatData[numFibers*3] = -yLoc;
  newMatData[numFibers*3+1] = zLoc;
  newMatData[numFibers*3+2] = Area;
  UniaxialMaterial *theMat = newFiber.getMaterial();
  newArray[numFibers] = theMat->getCopy();

  if (newArray[numFibers] == 0) {
    opserr << "RCFTSTLFiberSection3D::addFiber -- failed to get copy of a Material\n";
			  

    delete [] newArray;
    delete [] newMatData;
    return -1;
  }

  numFibers++;
  
  if (theMaterials != 0) {
    delete [] theMaterials;
    delete [] matData;
  }

  theMaterials = newArray;
  matData = newMatData;

  double Qz = 0.0;
  double Qy = 0.0;
  double A  = 0.0;

  // Recompute centroid
  for (i = 0; i < numFibers; i++) {
    yLoc = -matData[2*i];
    zLoc = matData[2*i+1];
    Area = matData[2*i+2];
    A  += Area;
    Qz += yLoc*Area;
    Qy += zLoc*Area;
  }

  yBar = -Qz/A;
  zBar = Qy/A;

  return 0;
}

// destructor:
RCFTSTLFiberSection3D::~RCFTSTLFiberSection3D()
{
  if (theMaterials != 0) {
    for (int i = 0; i < numFibers; i++)
      if (theMaterials[i] != 0)
	delete theMaterials[i];
      
    delete [] theMaterials;
  }

  if (matData != 0)
    delete [] matData;
}

int
RCFTSTLFiberSection3D::setTrialSectionDeformation (const Vector &deforms)
{
  int res = 0;
  e = deforms;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0;
  kData[3] = 0.0; kData[4] = 0.0; kData[5] = 0.0; 
  kData[6] = 0.0; kData[7] = 0.0; kData[8] = 0.0;
  kData[9] = 0.0; kData[10] = 0.0; kData[11] = 0.0;
  kData[12] = 0.0; kData[13] = 0.0; kData[14] = 0.0; kData[15] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0; sData[2] = 0.0; 

  int loc = 0;

  double d0 = deforms(0);
  double d1 = deforms(1);
  double d2 = deforms(2);

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++];// - yBar;
    double z = matData[loc++];// - zBar;
    double A = matData[loc++];

    // determine material strain and set it
    double strain = d0 - y*d1 - z*d2;
    double tangent, stress;
    res = theMat->setTrial(strain, stress, tangent);

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;
    double vas2as2 = vas2*z;
    double vas1as1 = vas1*y;

    kData[0] += value;
    kData[1] -= vas1;
    kData[2] -= vas2;
    kData[3]  = 0.0;
    kData[4] -= vas1;
    kData[5] += vas1as1;
    kData[6] += vas1as2;
    kData[7]  = 0.0;
    kData[8] -= vas2;
    kData[9] += vas1as2;
    kData[10]+= vas2as2;
    kData[11] = 0.0;
    kData[12] = 0.0;
    kData[13] = 0.0;
    kData[14] = 0.0;
    kData[15] = GJ;

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] -= fs0 * y;
    sData[2] -= fs0 * z;
  }

  ks(0,0) = kData[0];
  ks(0,1) = ks(1,0) = kData[1];
  ks(0,2) = ks(2,0) = kData[2];
  ks(1,1) = kData[5];
  ks(1,2) = ks(2,1) = kData[6];
  ks(2,2) = kData[10];
  ks(3,3) = GJ;

  //kData[0] = 449602.6;
  //kData[5] = 6968845.86;
  //kData[10] = 6968845.86;

  //ks(0,0) = 449602.6;
  //ks(0,1) = ks(1,0) = 0.0;
  //ks(0,2) = ks(2,0) = 0.0;
  //ks(1,1) = 6968845.86;
  //ks(1,2) = ks(2,1) = 0.0;
  //ks(2,2) = 6968845.86;
  //ks(3,3) = GJ;

  //kData[0] = 1000.0;
  //kData[5] = 1.0;
  //kData[10] = 1.0;

  //ks(0,0) = 1000.0;
  //ks(0,1) = ks(1,0) = 0.0;
  //ks(0,2) = ks(2,0) = 0.0;
  //ks(1,1) = 1.0;
  //ks(1,2) = ks(2,1) = 0.0;
  //ks(2,2) = 1.0;
  //ks(3,3) = GJ;

  //kData[0] = 43200000.0;
  //kData[5] = 14400000.0;
  //kData[10] = 14400000.0;

  //ks(0,0) = 43200000.0;
  //ks(0,1) = ks(1,0) = 0.0;
  //ks(0,2) = ks(2,0) = 0.0;
  //ks(1,1) = 14400000.0;
  //ks(1,2) = ks(2,1) = 0.0;
  //ks(2,2) = 14400000.0;
  //ks(3,3) = GJ;
  
  //kData[0] = 449509.3;
  //kData[5] = 6967684.16;
  //kData[10] = 6967684.16;
  
  //ks(0,0) = 449509.3;
  //ks(0,1) = ks(1,0) = 0.0;
  //ks(0,2) = ks(2,0) = 0.0;
  //ks(1,1) = 6967684.16;
  //ks(1,2) = ks(2,1) = 0.0;
  //ks(2,2) = 6967684.16;
  //ks(3,3) = GJ;

  //kData[0] = 300000.0;
  //kData[5] = 3000000.0;
  //kData[10]= 3000000.0;

  //ks(0,0) = 300000.0;
  //ks(0,1) = ks(1,0) = 0.0;
  //ks(0,2) = ks(2,0) = 0.0;
  //ks(1,1) = 3000000.0;
  //ks(1,2) = ks(2,1) = 0.0;
  //ks(2,2) = 3000000.0;
  //ks(3,3) = GJ;

  //kData[0] = 1884.9;
  //kData[5] = 9.27;
  //kData[10]= 9.27;

  //ks(0,0) = 1884.9;
  //ks(0,1) = ks(1,0) = 0.0;
  //ks(0,2) = ks(2,0) = 0.0;
  //ks(1,1) = 9.27;
  //ks(1,2) = ks(2,1) = 0.0;
  //ks(2,2) = 9.27;
  //ks(3,3) = GJ;

  //kData[0] = 43200000;
  //kData[5] = 14400000;
  //kData[10]= 14400000;

  //ks(0,0) = 43200000;
  //ks(0,1) = ks(1,0) = 0.0;
  //ks(0,2) = ks(2,0) = 0.0;
  //ks(1,1) = 14400000;
  //ks(1,2) = ks(2,1) = 0.0;
  //ks(2,2) = 14400000;
  //ks(3,3) = GJ;

  return res;
}

const Matrix&
RCFTSTLFiberSection3D::getInitialTangent(void)
{
  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0;
  kData[3] = 0.0; kData[4] = 0.0; kData[5] = 0.0;
  kData[6] = 0.0; kData[7] = 0.0; kData[8] = 0.0;
  kData[9] = 0.0; kData[10] = 0.0; kData[11] = 0.0;
  kData[12] = 0.0; kData[13] = 0.0; kData[14] = 0.0; kData[15] = 0.0;

  sData[0] = 0.0; sData[1] = 0.0; sData[2] = 0.0;

  int loc = 0;

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++]; //- yBar;
    double z = matData[loc++]; //- zBar;
    double A = matData[loc++];

    double tangent = theMat->getInitialTangent();
    
    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;
    double vas2as2 = vas2*z;
    double vas1as1 = vas1*y;

    kData[0] += value;
    kData[1] -= vas1;
    kData[1]  = 0.0;
    kData[2] -= vas2;
    kData[2]  = 0.0;
    kData[3]  = 0.0;
    kData[4] -= vas1;
    kData[4]  = 0.0;
    kData[5] += vas1as1;
    kData[6] += vas1as2;
    kData[6]  = 0.0;
    kData[7]  = 0.0;
    kData[8] -= vas2;
    kData[8]  = 0.0; 
    kData[9] += vas1as2;
    kData[9]  = 0.0;
    kData[10]+= vas2as2;
    kData[11] = 0.0;
    kData[12] = 0.0;
    kData[13] = 0.0;
    kData[14] = 0.0;
    kData[15] = GJ;
  }

  ks(0,0) = kData[0];
  ks(0,1) = ks(1,0) = kData[1];
  ks(0,2) = ks(2,0) = kData[2];
  ks(1,1) = kData[5];
  ks(1,2) = ks(2,1) = kData[6];
  ks(2,2) = kData[10];
  ks(3,3) = GJ;

  return ks;
}

const Vector&
RCFTSTLFiberSection3D::getSectionDeformation(void)
{
  return e;
}

const Matrix&
RCFTSTLFiberSection3D::getSectionTangent(void)
{
  ks(0,0) = kData[0];
  ks(0,1) = ks(1,0) = kData[1];
  ks(0,2) = ks(2,0) = kData[2];
  ks(1,1) = kData[5];
  ks(1,2) = ks(2,1) = kData[6];
  ks(2,2) = kData[10];
  ks(3,3) = GJ;
  return ks;
}

const Vector&
RCFTSTLFiberSection3D::getStressResultant(void)
{
  s(0) = sData[0];
  s(1) = sData[1];
  s(2) = sData[2];

  s(3) = GJ*e(3);

  return s;
}

SectionForceDeformation*
RCFTSTLFiberSection3D::getCopy(void)
{
  RCFTSTLFiberSection3D *theCopy = new RCFTSTLFiberSection3D();
  theCopy->setTag(this->getTag());

  theCopy->numFibers = numFibers;

  if (numFibers != 0) {
    theCopy->theMaterials = new UniaxialMaterial *[numFibers];

    if (theCopy->theMaterials == 0) {
      opserr << "RCFTSTLFiberSection3D::RCFTSTLFiberSection3D -- failed to allocate Material pointers\n";
      exit(-1);
    }

    theCopy->matData = new double [numFibers*3];

    if (theCopy->matData == 0) {
      opserr << "RCFTSTLFiberSection3D::RCFTSTLFiberSection3D -- failed to allocate double array for material data\n";
      exit(-1);
    }    
    for (int i = 0; i < numFibers; i++) {
      theCopy->matData[i*3] = matData[i*3];
      theCopy->matData[i*3+1] = matData[i*3+1];
      theCopy->matData[i*3+2] = matData[i*3+2];
      theCopy->theMaterials[i] = theMaterials[i]->getCopy();

      if (theCopy->theMaterials[i] == 0) {
	opserr << "RCFTSTLFiberSection3D::getCopy -- failed to get copy of a Material\n";
	exit(-1);
      }
    }    
  }

  theCopy->eCommit = eCommit;
  theCopy->e = e;
  theCopy->yBar = yBar;
  theCopy->zBar = zBar;

  for (int i=0; i<16; i++)
    theCopy->kData[i] = kData[i];

  theCopy->sData[0] = sData[0];
  theCopy->sData[1] = sData[1];
  theCopy->sData[2] = sData[2];

  theCopy->GJ = GJ;

  return theCopy;
}

const ID&
RCFTSTLFiberSection3D::getType ()
{
  return code;
}

int
RCFTSTLFiberSection3D::getOrder () const
{
  return 4;
}

int
RCFTSTLFiberSection3D::commitState(void)
{
  int err = 0;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->commitState();

  //ofstream stlb4;
  //stlb4.open("stlb4.dat",ios::app); 

  //ofstream stlb12;
  //stlb12.open("stlb12.dat",ios::app);

  //ofstream stlb6;
  //stlb6.open("stlb6.dat",ios::app);

//  for( int i = 0; i < 16; i++ ){
//    stringstream number;
//    number << i;
//    string conc = "stlb";
//    string type = ".dat";
//    string name = conc + number.str() + type;
//    ofstream output;
//    output.open(name.c_str(),ios::app);
//    output<<theMaterials[i]->getStrain()<<" "<<theMaterials[i]->getStress()<<" "<<theMaterials[i]->getTangent()<<endl;
//  }

  //stlb4<<theMaterials[4]->getStrain()<<" "<<theMaterials[4]->getStress()<<" "<<theMaterials[4]->getTangent()<<endl;
  //stlb12<<theMaterials[12]->getStrain()<<" "<<theMaterials[12]->getStress()<<" "<<theMaterials[12]->getTangent()<<endl;
  //stlb6<<theMaterials[6]->getStrain()<<" "<<theMaterials[6]->getStress()<<" "<<theMaterials[6]->getTangent()<<endl;

  eCommit = e;

  return err;
}

int
RCFTSTLFiberSection3D::revertToLastCommit(void)
{
  int err = 0;

  // Last committed section deformations
  e = eCommit;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0;
  kData[3] = 0.0; kData[4] = 0.0; kData[5] = 0.0;
  kData[6] = 0.0; kData[7] = 0.0; kData[8] = 0.0;
  kData[9] = 0.0; kData[10] = 0.0; kData[11] = 0.0;
  kData[12] = 0.0; kData[13] = 0.0; kData[14] = 0.0;
  kData[15] = 0.0;

  sData[0] = 0.0; sData[1] = 0.0; sData[2] = 0.0; 

  int loc = 0;

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];

    // invoke revertToLast on the material
    err += theMat->revertToLastCommit();

    double tangent = theMat->getTangent();
    double stress = theMat->getStress();

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;
    double vas2as2 = vas2*z;
    double vas1as1 = vas1*y;

    kData[0] += value;
    kData[1] -= vas1;
    kData[2] -= vas2;
    kData[3]  = 0.0;
    kData[4] -= vas1;
    kData[5] += vas1as1;
    kData[6] += vas1as2;
    kData[7]  = 0.0;
    kData[8] -= vas2;
    kData[9] += vas1as2;
    kData[10]+= vas2as2;
    kData[11] = 0.0;
    kData[12] = 0.0;
    kData[13] = 0.0;
    kData[14] = 0.0;
    kData[15] = GJ;

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] -= fs0 * y;
    sData[2] -= fs0 * z;
  }

  ks(0,0) = kData[0];
  ks(0,1) = ks(1,0) = kData[1];
  ks(0,2) = ks(2,0) = kData[2];
  ks(1,1) = kData[5];
  ks(1,2) = ks(2,1) = kData[6];
  ks(2,2) = kData[10];
  ks(3,3) = GJ;

  return err;
}

int
RCFTSTLFiberSection3D::revertToStart(void)
{
  // revert the fibers to start    
  int err = 0;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0;
  kData[3] = 0.0; kData[4] = 0.0; kData[5] = 0.0;
  kData[6] = 0.0; kData[7] = 0.0; kData[8] = 0.0;
  kData[9] = 0.0; kData[10] = 0.0; kData[11] = 0.0;
  kData[12] = 0.0; kData[13] = 0.0; kData[14] = 0.0;
  kData[15] = 0.0;

  sData[0] = 0.0; sData[1] = 0.0; sData[2] = 0.0; 

  int loc = 0;

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];

    // invoke revertToStart on the material
    err += theMat->revertToStart();

    double tangent = theMat->getTangent();
    double stress = theMat->getStress();

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;
    double vas2as2 = vas2*z;
    double vas1as1 = vas1*y;

    kData[0] += value;
    kData[1] -= vas1;
    kData[1]  = 0.0;
    kData[2] -= vas2;
    kData[2]  = 0.0;
    kData[3]  = 0.0;
    kData[4] -= vas1;
    kData[4]  = 0.0;
    kData[5] += vas1as1;
    kData[6] += vas1as2;
    kData[6]  = 0.0;
    kData[7]  = 0.0;
    kData[8] -= vas2;
    kData[8]  = 0.0;
    kData[9] += vas1as2;
    kData[9]  = 0.0;
    kData[10]+= vas2as2;
    kData[11] = 0.0;
    kData[12] = 0.0;
    kData[13] = 0.0;
    kData[14] = 0.0;
    kData[15] = GJ;

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += fs0 * y;
    sData[2] += fs0 * z;
  }

  return err;
}

int
RCFTSTLFiberSection3D::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  return res;
}

int
RCFTSTLFiberSection3D::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  return res;
}

void
RCFTSTLFiberSection3D::Print(OPS_Stream &s, int flag)
{
  s << "\nRCFTFiberSection3D, tag: " << this->getTag() << endln;
  s << "\tSection code: " << code;
  s << "\tNumber of Fibers: " << numFibers << endln;
  s << "\tCentroid: (" << -yBar << ", " << zBar << ')' << endln;
  s << "\tTorsional Stiffness: " << GJ << endln;

  if (flag == 1) {
    int loc = 0;
    for (int i = 0; i < numFibers; i++) {
      s << "\nLocation (y, z) = (" << -matData[loc++] << ", " << matData[loc++] << ")";
      s << "\nArea = " << matData[loc++] << endln;
      theMaterials[i]->Print(s, flag);
    }
  }
}

Response*
RCFTSTLFiberSection3D::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  // See if the response is one of the defaults
  Response *res = SectionForceDeformation::setResponse(argv, argc, output);
  if (res != 0)
    return res;
  
  // Check if fiber response is requested
  else if (strcmp(argv[0],"fiber") == 0) {
    int key = numFibers;
    int passarg = 2;
    
    if (argc <= 2)          // not enough data input
      return 0;
    
    if (argc <= 3)		  // fiber number was input directly
      key = atoi(argv[1]);
    
    if (argc > 4) {         // find fiber closest to coord. with mat tag
      int matTag = atoi(argv[3]);
      double yCoord = atof(argv[1]);
      double zCoord = atof(argv[2]);
      double closestDist;
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
    
    if (key < numFibers)
      return theMaterials[key]->setResponse(&argv[passarg],argc-passarg,output);
    else
      return 0;
  }
  
  // otherwise response quantity is unknown for the FiberSection class
  else
    return 0;
}


int 
RCFTSTLFiberSection3D::getResponse(int responseID, Information &sectInfo)
{
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  return SectionForceDeformation::getResponse(responseID, sectInfo);
}

int
RCFTSTLFiberSection3D::setParameter (const char **argv, int argc, Parameter &param)
{
	// Initial declarations
	int ok = -1;

	// A material parameter
	if (strcmp(argv[0],"material") == 0) {

		// Get the tag of the material
		int paramMatTag = atoi(argv[1]);

		// Loop over fibers to find the right material(s)
		for (int i=0; i<numFibers; i++) {
			if (paramMatTag == theMaterials[i]->getTag()) {
				ok = theMaterials[i]->setParameter(&argv[2], argc-2, param);
			}
		}
		if (ok<0) {
		  opserr << "RCFTSTLFiberSection3D::setParameter() - could not set parameter. " << endln;
		  return -1;
		}
		else {
			return ok + 100;
		}
	} 
	else
		return -1;
}

int
RCFTSTLFiberSection3D::updateParameter (int parameterID, Information &info)
{
	int ok = -1;

	switch (parameterID) {
	case 1:
		return -1;
	default:
		if (parameterID >= 100) {
			ID *paramIDPtr;
			paramIDPtr = info.theID;
			ID paramID = (*paramIDPtr);
			int paramMatrTag = paramID(1);

			for (int i=0; i<numFibers; i++) {
				if (paramMatrTag == theMaterials[i]->getTag()) {
					ok =theMaterials[i]->updateParameter(parameterID-100, info);
				}
			}
			if (ok < 0) {
				opserr << "RCFTSTLFiberSection3D::updateParameter() - could not update parameter. " << endln;
				return ok;
			}
			else {
				return ok;
			}
		}
		else
			return -1;
	}
}
